"""Compare mean-ES allocations on a fixed downloaded data snapshot."""
import json, math, statistics
from pathlib import Path
from prepare_data import es
root=Path(__file__).parent
data=json.loads((root/'returns.json').read_text())
a=data['meta']['columns'].index('HiTec')+1
b=data['meta']['columns'].index('Utils')+1
train=[r for r in data['rows'] if '2010-01-01'<=r[0]<'2020-01-01']
test=[r for r in data['rows'] if r[0]>='2020-01-01']
alpha=.975; budget=2.5
z=statistics.NormalDist().inv_cdf(alpha)
def normal_es(r):return -statistics.mean(r)+statistics.stdev(r)*math.exp(-z*z/2)/math.sqrt(2*math.pi)/(1-alpha)
def returns(rows,w):return [w*r[a]+(1-w)*r[b] for r in rows]
result={'training':'2010–2019','testing':'2020–2025','training_n':len(train),'test_n':len(test),'alpha':alpha,'budget_percent':budget,'methods':[]}
for name,risk in [('Historical ES',lambda r:es(r,alpha)),('Normal ES',normal_es)]:
    candidates=[]
    for j in range(101):
        w=j/100;r=returns(train,w)
        if risk(r)<=budget:candidates.append((statistics.mean(r),w))
    mu,w=max(candidates)
    tr,te=returns(train,w),returns(test,w)
    result['methods'].append({'method':name,'technology_weight':w,'training_mean':mu,'model_training_es':risk(tr),'historical_training_es':es(tr,alpha),'test_mean':statistics.mean(te),'test_es':es(te,alpha)})
(root/'allocation-results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
