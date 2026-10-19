"""Matched-target ES/Markowitz and full-evaluation-period hindsight benchmark.
Run alongside returns.json: python hindsight_study.py (NumPy, SciPy).
All numbers use daily percent returns. No future data enter trained weights.
"""
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.optimize import linprog, minimize
root=Path(__file__).parent
D=json.loads((root/'returns.json').read_text())
train=np.array([r[1:] for r in D['rows'] if '2010-01-01'<=r[0]<'2020-01-01'])
test=np.array([r[1:] for r in D['rows'] if r[0]>='2020-01-01'])
p=train.shape[1];mu=train.mean(0);cov=np.cov(train,rowvar=False);eq=np.ones(p)/p

def es(r,a):
    l=np.sort(-r)[::-1];f=np.clip(len(l)*(1-a)-np.arange(len(l)),0,1)
    return float(l@f/f.sum())

def summary(w,a):
    tr=train@w;te=test@w
    return dict(weights=w.tolist(),mean=float(tr.mean()),sd=float(tr.std(ddof=1)),es=es(tr,a),test_mean=float(te.mean()),test_sd=float(te.std(ddof=1)),test_es=es(te,a))

def es_solver(x,a):
    n=len(x);m=x.mean(0)
    obj=np.r_[np.zeros(p),1,np.full(n,1/(n*(1-a)))]
    A=sparse.hstack([-sparse.csr_matrix(x),-np.ones((n,1)),-sparse.eye(n)],format='csr')
    E=sparse.csr_matrix(np.r_[np.ones(p),np.zeros(n+1)].reshape(1,-1))
    bounds=[(0,1)]*p+[(None,None)]+[(0,None)]*n
    def solve(target=None):
        aub=A if target is None else sparse.vstack([A,sparse.csr_matrix(np.r_[-m,np.zeros(n+1)].reshape(1,-1))])
        bub=np.zeros(n) if target is None else np.r_[np.zeros(n),-target]
        r=linprog(obj,A_ub=aub,b_ub=bub,A_eq=E,b_eq=[1],bounds=bounds,method='highs')
        assert r.success,r.message
        w=r.x[:p];assert abs(w.sum()-1)<1e-7
        if target is not None:assert m@w>=target-1e-7
        assert abs(es(x@w,a)-r.fun)<1e-6
        return w
    return solve

def markowitz(target=None):
    constraints=[{'type':'eq','fun':lambda w:w.sum()-1,'jac':lambda w:np.ones(p)}]
    if target is not None:constraints.append({'type':'ineq','fun':lambda w:mu@w-target,'jac':lambda w:mu})
    r=minimize(lambda w:w@cov@w,eq,jac=lambda w:2*cov@w,bounds=[(0,1)]*p,constraints=constraints,method='SLSQP',options={'ftol':1e-11,'maxiter':1000})
    assert r.success,r.message
    if target is not None:assert mu@r.x>=target-1e-7
    return r.x

minvar=markowitz()
out=dict(training='2010–2019',evaluation='2020–2025',columns=D['meta']['columns'],names=D['meta']['names'],units='daily percent',levels={})
for a in [.95,.975,.99]:
    fit=es_solver(train,a);oracle=es_solver(test,a)
    # Common target range starts above both global risk-minimum means.
    lo=max(float(mu@minvar),float(mu@fit()))
    qs=np.linspace(lo,float(mu.max()),25)
    records=[]
    for q in qs:
        h=summary(fit(q),a);m=summary(markowitz(q),a)
        assert h['es']<=m['es']+1e-6
        assert m['sd']<=h['sd']+1e-6
        for s in [h,m]:
            o=summary(oracle(s['test_mean']),a)
            s['oracle']=o
            s['hindsight_gap']=s['test_es']-o['test_es']
            assert s['hindsight_gap']>=-1e-6
            assert o['test_mean']>=s['test_mean']-1e-7
        records.append(dict(target=float(q),historical=h,markowitz=m))
    low=float((test@oracle()).mean());high=float(test.mean(0).max())
    curve=[summary(oracle(q),a) for q in np.linspace(low,high,41)]
    # Include exact matched oracle points so comparisons don't rely on interpolation.
    curve.extend(s[k]['oracle'] for s in records for k in ['historical','markowitz'])
    curve.sort(key=lambda s:s['test_mean'])
    out['levels'][str(a)]=dict(records=records,oracle=curve,equal=summary(eq,a))
    print('Verified matched targets and hindsight bound',a,flush=True)
(root/'hindsight-results.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
(root/'hindsight-data.js').write_text('window.HINDSIGHT_STUDY='+json.dumps(out,separators=(',',':'))+';\n')
