"""36-calendar-month estimation; monthly buy-and-hold backtest, 2020–2025.
Run alongside returns.json with NumPy/SciPy: python rolling_study.py.
Daily percent units. Ten industries plus zero-return cash; alpha=.975.
Oracle: exact next-month terminal-return optimum under that month's empirical ES.
For <40 daily observations, 97.5% empirical ES equals maximum loss, making the
buy-and-hold oracle a linear program in opening holdings. No transaction costs.
"""
import json, hashlib
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from scipy import sparse
from scipy.optimize import linprog, minimize
from scipy.stats import norm
ROOT=Path(__file__).parent
ALPHA=.975
K=float(norm.pdf(norm.ppf(ALPHA))/(1-ALPHA))
BUDGETS=np.linspace(.5,4.5,17)

def es(r):
    losses=np.sort(-np.asarray(r))[::-1];mass=len(losses)*(1-ALPHA)
    return float(losses@np.clip(mass-np.arange(len(losses)),0,1)/mass)

def held_returns(x,w):
    g=np.cumprod(1+x/100,axis=0)@w
    return 100*(g/np.r_[1,g[:-1]]-1)

def summary(x,w):
    r=x@w;return dict(mean=float(r.mean()),es=es(r),normal_es=float(-r.mean()+K*r.std(ddof=1)),sd=float(r.std(ddof=1)))

def month_job(task):
    month,dates,train,train_dates,future=task
    train=np.array(train);future=np.array(future);n,p=train.shape
    mu=train.mean(0);cov=np.cov(train,rowvar=False);cash=np.eye(p)[-1]
    # Mean maximization under empirical ES: Rockafellar–Uryasev LP.
    A=sparse.hstack([-sparse.csr_matrix(train),-np.ones((n,1)),-sparse.eye(n)],format='csr')
    risk=np.r_[np.zeros(p),1,np.full(n,1/(n*(1-ALPHA)))]
    aub=sparse.vstack([A,sparse.csr_matrix(risk.reshape(1,-1))])
    eq=sparse.csr_matrix(np.r_[np.ones(p),np.zeros(n+1)].reshape(1,-1))
    bounds=[(0,1)]*p+[(None,None)]+[(0,None)]*n
    G=np.cumprod(1+future/100,axis=0);Gprev=np.vstack([np.ones(p),G[:-1]])
    assert len(future)*(1-ALPHA)<1
    out=[]
    for b in BUDGETS:
        lp=linprog(np.r_[-mu,np.zeros(n+1)],A_ub=aub,b_ub=np.r_[np.zeros(n),b],A_eq=eq,b_eq=[1],bounds=bounds,method='highs')
        assert lp.success,lp.message
        weights={'historical':lp.x[:p]}
        for key in ['normal','markowitz']:
            v=b/K
            def riskfn(w):
                sd=np.sqrt(max(0,w@cov@w));return -mu@w+K*sd if key=='normal' else sd
            def riskjac(w):
                sd=np.sqrt(max(1e-24,w@cov@w));return -mu+K*(cov@w)/sd if key=='normal' else (cov@w)/sd
            limit=b if key=='normal' else v
            init=np.full(p,.1/p);init[-1]+=0.9
            cons=[{'type':'eq','fun':lambda w:w.sum()-1,'jac':lambda w:np.ones(p)}, {'type':'ineq','fun':lambda w:limit-riskfn(w),'jac':lambda w:-riskjac(w)}]
            best=None
            for scale in [1.,10.]:
                fit=minimize(lambda w:-scale*mu@w,init,jac=lambda w:-scale*mu,bounds=[(0,1)]*p,constraints=cons,method='SLSQP',options={'ftol':1e-11,'maxiter':1500})
                w=np.maximum(0,fit.x);w=w/w.sum()
                # Exact homogeneous scaling into cash removes numerical infeasibility.
                rr=riskfn(w)
                if rr>limit:
                    w*=limit/rr;w[-1]+=1-w.sum()
                g=riskjac(w)
                upper=linprog(-mu,A_ub=g.reshape(1,-1),b_ub=[limit-riskfn(w)+g@w],A_eq=np.ones((1,p)),b_eq=[1],bounds=[(0,1)]*p,method='highs')
                # Convex risk tangent defines a relaxation: its return bounds the optimum.
                if upper.success and -upper.fun-mu@w<1e-6:
                    best=w;break
                init=w*.99;init[-1]+=.01
            assert best is not None,(month,key,b,fit.message)
            assert riskfn(best)<=limit+1e-7
            weights[key]=best
        # Exact buy-and-hold ES constraint: every daily gross return >= 1-b/100.
        oracle=linprog(-G[-1],A_ub=(1-b/100)*Gprev-G,b_ub=np.zeros(len(G)),A_eq=np.ones((1,p)),b_eq=[1],bounds=[(0,1)]*p,method='highs')
        assert oracle.success,oracle.message
        weights['oracle']=oracle.x
        row={}
        for key,w in weights.items():
            w=np.clip(w,0,1);w=w/w.sum();r=held_returns(future,w)
            row[key]=dict(weights=w.tolist(),forecast=summary(train,w),returns=r.tolist())
            assert abs(w.sum()-1)<1e-8
            if key=='historical':assert es(train@w)<=b+1e-6
            if key=='oracle':
                assert es(r)<=b+1e-5
                # Every feasible comparator has no greater month terminal wealth.
                for wk in weights.values():
                    if es(held_returns(future,wk))<=b+1e-6:assert G[-1]@w>=G[-1]@wk-1e-6
        out.append(row)
    return dict(month=month,dates=dates,train_start=train_dates[0],train_end=train_dates[-1],train_n=n,records=out)

def build():
    raw=(ROOT/'returns.json').read_bytes();D=json.loads(raw);rows=D['rows']
    months=sorted({r[0][:7] for r in rows if '2020-01-01'<=r[0]<'2026-01-01'})
    tasks=[]
    for m in months:
        y,mo=map(int,m.split('-'));start=f'{y-3:04d}-{mo:02d}-01';end=m+'-01'
        tr=[r for r in rows if start<=r[0]<end];te=[r for r in rows if r[0][:7]==m]
        assert tr[-1][0]<te[0][0]
        tasks.append((m,[r[0] for r in te],[r[1:]+[0] for r in tr],[r[0] for r in tr],[r[1:]+[0] for r in te]))
    results=[]
    with ProcessPoolExecutor(max_workers=4) as pool:
        for i,r in enumerate(pool.map(month_job,tasks)):
            results.append(r)
            print(f'Calculated {i+1}/72: {r["month"]}',flush=True)
    curves=[]
    for i,b in enumerate(BUDGETS):
        record=dict(budget=float(b),volatility_budget=float(b/K))
        for key in ['historical','normal','markowitz','oracle']:
            rs=np.concatenate([m['records'][i][key]['returns'] for m in results]);wealth=np.r_[1,np.cumprod(1+rs/100)];draw=1-wealth/np.maximum.accumulate(wealth)
            record[key]=dict(forecast={k:float(np.mean([m['records'][i][key]['forecast'][k] for m in results])) for k in ['mean','es','normal_es','sd']},mean=float(rs.mean()),es=es(rs),sd=float(rs.std(ddof=1)),terminal=float(wealth[-1]),cagr=float(wealth[-1]**(1/6)-1),max_drawdown=float(draw.max()),returns=rs.tolist(),wealth=wealth.tolist(),drawdown=draw.tolist())
        curves.append(record)
    out=dict(alpha=ALPHA,k=K,window_months=36,rebalance='monthly',cash_return=0,cost_bps=0,source_sha256=hashlib.sha256(raw).hexdigest(),columns=D['meta']['columns']+['Cash'],names={**D['meta']['names'],'Cash':'Cash (0% return)'},dates=['2019-12-31']+[d for m in results for d in m['dates']],months=results,curves=curves)
    s=json.dumps(out,separators=(',',':'))
    (ROOT/'rolling-results.json').write_text(s+'\n');(ROOT/'rolling-data.js').write_text('window.ROLLING_STUDY='+s+';\n')
    print('Saved verified rolling study',len(s),'bytes',flush=True)
if __name__=='__main__':build()
