"""Ten-industry mean–ES study. Run: python -m pip install numpy scipy; python portfolio_study.py.
Requires returns.json beside this script. Returns and risks are in daily percent.
Published 10Ind universe: DeMiguel et al. (2009), Management Science, Table 1.
This new daily-data exercise uses 2010–19 training and untouched 2020–25 evaluation.
"""
import json
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.optimize import linprog, minimize
from scipy.stats import norm

root = Path(__file__).parent
data = json.loads((root/'returns.json').read_text())
train = np.array([r[1:] for r in data['rows'] if '2010-01-01' <= r[0] < '2020-01-01'])
test = np.array([r[1:] for r in data['rows'] if r[0] >= '2020-01-01'])
n, p = train.shape
mu, cov = train.mean(axis=0), np.cov(train, rowvar=False)
eq = np.ones(p)/p

def contributions(x, w, alpha):
    losses = -x@w
    order = np.argsort(-losses, kind='stable')
    mass = len(x)*(1-alpha)
    f = np.clip(mass-np.arange(len(x)), 0, 1)
    return ((-x[order]*w)*f[:, None]).sum(axis=0)/mass

def stats(w, alpha):
    cs = contributions(train, w, alpha)
    return dict(weights=w.tolist(), mean=float(mu@w), es=float(cs.sum()),
                normal_es=float(-mu@w+norm.pdf(norm.ppf(alpha))/(1-alpha)*np.sqrt(w@cov@w)),
                contributions=cs.tolist(), test_mean=float((test@w).mean()),
                test_es=float(contributions(test, w, alpha).sum()))

def study(alpha):
    # Variables: ten weights, unrestricted VaR threshold, n nonnegative excesses.
    c = np.r_[np.zeros(p), 1, np.full(n, 1/(n*(1-alpha)))]
    A = sparse.hstack([-sparse.csr_matrix(train), -np.ones((n,1)), -sparse.eye(n)], format='csr')
    E = sparse.csr_matrix(np.r_[np.ones(p),np.zeros(n+1)].reshape(1,-1))
    bounds = [(0,1)]*p+[(None,None)]+[(0,None)]*n
    def lp(budget=None):
        obj = c if budget is None else np.r_[-mu,np.zeros(n+1)]
        aub = A if budget is None else sparse.vstack([A,sparse.csr_matrix(c.reshape(1,-1))])
        bub = np.zeros(n) if budget is None else np.r_[np.zeros(n),budget]
        res = linprog(obj,A_ub=aub,b_ub=bub,A_eq=E,b_eq=[1],bounds=bounds,method='highs')
        if res.status == 2: return None
        assert res.success,res.message
        return res.x[:p]
    k = norm.pdf(norm.ppf(alpha))/(1-alpha)
    risk = lambda w: -mu@w+k*np.sqrt(w@cov@w)
    grad = lambda w: -mu+k*(cov@w)/np.sqrt(w@cov@w)
    constraints = [{'type':'eq','fun':lambda w:w.sum()-1,'jac':lambda w:np.ones(p)}]
    mn = minimize(risk,eq,jac=grad,method='SLSQP',bounds=[(0,1)]*p,constraints=constraints,options={'ftol':1e-12,'maxiter':500})
    assert mn.success,mn.message
    min_hist = stats(lp(),alpha)['es']
    maxret = np.eye(p)[np.argmax(mu)]
    top = max(stats(maxret,alpha)['es'],risk(maxret))
    budgets = np.linspace(min(min_hist,risk(mn.x))*.95, top, 51)
    historical, gaussian = [], []
    for b in budgets:
        w = lp(b)
        historical.append(None if w is None else stats(w,alpha))
        if b < risk(mn.x)-1e-9:
            gaussian.append(None)
            continue
        result = minimize(lambda w:-100*mu@w,mn.x,jac=lambda w:-100*mu,method='SLSQP',bounds=[(0,1)]*p,
                          constraints=constraints+[{'type':'ineq','fun':lambda w:b-risk(w),'jac':lambda w:-grad(w)}],
                          options={'ftol':1e-7,'maxiter':1000})
        assert result.success,(alpha,b,result.message,risk(result.x)-b)
        assert risk(result.x)<=b+1e-6
        gaussian.append(stats(result.x,alpha))
    for points,key in [(historical,'es'),(gaussian,'normal_es')]:
        for b,s in zip(budgets,points):
            if s is not None:
                assert abs(sum(s['weights'])-1)<1e-7
                assert min(s['weights'])>=-1e-7
                assert s[key]<=b+1e-6
                assert abs(sum(s['contributions'])-s['es'])<1e-10
    return dict(budgets=budgets.tolist(), historical=historical, gaussian=gaussian,
                equal=stats(eq,alpha), minimum_historical_es=min_hist,
                minimum_gaussian_es=float(risk(mn.x)))

out = dict(columns=data['meta']['columns'], names=data['meta']['names'], training='2010–2019',testing='2020–2025',
           training_n=n,test_n=len(test),units='daily percent',source_sha256=data['meta']['zip_sha256'],levels={})
for alpha in [.95,.975,.99]:
    out['levels'][str(alpha)]=study(alpha)
    print('Solved',alpha,flush=True)
(root/'portfolio-results.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
(root/'portfolio-data.js').write_text('window.PORTFOLIO_STUDY='+json.dumps(out,separators=(',',':'))+';\n')
