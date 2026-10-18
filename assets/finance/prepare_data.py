"""Download French industry returns and reproduce the notebook's two-asset study.
Uses only the Python standard library. Run from any directory.
"""
import csv, hashlib, io, json, math, statistics, urllib.request, zipfile
from pathlib import Path
ROOT=Path(__file__).parent
URL='https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/10_Industry_Portfolios_daily_CSV.zip'
NAMES={'NoDur':'Consumer non-durables','Durbl':'Consumer durables','Manuf':'Manufacturing','Enrgy':'Energy','HiTec':'Technology','Telcm':'Telecom','Shops':'Shops','Hlth':'Healthcare','Utils':'Utilities','Other':'Other'}

def es(values, alpha=.975):
    losses=sorted((-r for r in values),reverse=True)
    mass=len(losses)*(1-alpha)
    k=int(math.floor(mass)); f=mass-k
    return (sum(losses[:k])+(f*losses[k] if k<len(losses) else 0))/mass

def main():
    raw=urllib.request.urlopen(URL).read()
    z=zipfile.ZipFile(io.BytesIO(raw))
    text=z.read(z.namelist()[0]).decode()
    section=text.split('Average Value Weighted Returns -- Daily')[1].split('Average Equal Weighted Returns -- Daily')[0]
    reader=csv.reader(io.StringIO(section.strip()))
    header=[x.strip() for x in next(reader)][1:]
    rows=[]
    for r in reader:
        if not r or not r[0].strip().isdigit(): continue
        day=r[0].strip()
        if not '20000101'<=day<='20251231': continue
        v=[float(x) for x in r[1:]]
        if len(v)!=10 or any(x in (-99.99,-999) for x in v): continue
        rows.append([f'{day[:4]}-{day[4:6]}-{day[6:]}',*v])
    assert len(rows)>6000 and len(set(r[0] for r in rows))==len(rows)
    assert rows==sorted(rows)
    meta={'source':URL,'source_header':text.splitlines()[0], 'download_date':'2026-10-18','zip_sha256':hashlib.sha256(raw).hexdigest(),'units':'daily simple returns, percent','weighting':'value-weighted industry portfolios','columns':header,'names':NAMES,'start':rows[0][0],'end':rows[-1][0],'n':len(rows)}
    data={'meta':meta,'rows':rows}
    (ROOT/'returns.json').write_text(json.dumps(data,separators=(',',':'))+'\n')
    (ROOT/'data.js').write_text('window.FINANCE_DATA='+json.dumps(data,separators=(',',':'))+';\n')
    with (ROOT/'returns.csv').open('w',newline='') as f:
        w=csv.writer(f); w.writerow(['Date',*header]); w.writerows(rows)
    # All statistics in percentage points. Select weights on training dates only.
    a,b=header.index('HiTec')+1,header.index('Utils')+1
    train=[r for r in rows if r[0]<'2020-01-01']; test=[r for r in rows if r[0]>='2020-01-01']
    def series(rs,w):return [w*r[a]+(1-w)*r[b] for r in rs]
    grid=[i/100 for i in range(101)]
    wvar=min(grid,key=lambda w:statistics.variance(series(train,w)))
    wes=min(grid,key=lambda w:es(series(train,w)))
    report=[]
    for name,w in [('50/50',.5),('Minimum variance',wvar),('Minimum ES',wes)]:
        tr,te=series(train,w),series(test,w)
        report.append({'method':name,'technology_weight':w,'train_sd':statistics.stdev(tr),'train_es':es(tr),'test_sd':statistics.stdev(te),'test_es':es(te),'test_mean':statistics.mean(te)})
    whole=series(rows,.5); mu=statistics.mean(whole); sd=statistics.stdev(whole)
    n3=sum(abs(x-mu)>3*sd for x in whole)
    out={'train_n':len(train),'test_n':len(test),'alpha':.975,'pair':['HiTec','Utils'],'allocations':report,'default':{'n':len(whole),'sd':sd,'over_3sd':n3,'normal_expected':len(whole)*math.erfc(3/math.sqrt(2))}}
    (ROOT/'analysis.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(meta,indent=2));print(json.dumps(out,indent=2))

if __name__=='__main__': main()
