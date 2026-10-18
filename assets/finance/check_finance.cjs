// Compare browser-side formulas with the independently calculated Python report.
const assert=require('node:assert/strict');
const fs=require('node:fs');
const path=require('node:path');
const F=require('./panel.js');
const read=p=>JSON.parse(fs.readFileSync(path.join(__dirname,p)));
const {rows,meta}=read('returns.json'),truth=read('analysis.json');
const ia=meta.columns.indexOf('HiTec')+1,ib=meta.columns.indexOf('Utils')+1;
for(const t of truth.allocations){
  for(const [tag,select] of [['train',r=>r[0]<'2020-01-01'],['test',r=>r[0]>='2020-01-01']]){
    const r=rows.filter(select).map(r=>t.technology_weight*r[ia]+(1-t.technology_weight)*r[ib]);
    assert(Math.abs(F.sd(r)-t[tag+'_sd'])<1e-10);
    assert(Math.abs(F.es(r,.975)-t[tag+'_es'])<1e-10);
  }
}
assert(Math.abs(F.invnorm(.975)-1.95996398454)<1e-7);
assert(Math.abs(F.normalTail(3)-.002699796063)<2e-7);
assert.equal(F.es([-10,-5,0,5],.5),7.5);
assert.equal(F.es([-10,-5,0,5],.625),12.5/1.5);
const full=rows.map(r=>.5*r[ia]+.5*r[ib]);
const m=F.mean(full),s=F.sd(full);
assert.equal(full.filter(x=>Math.abs(x-m)>3*s).length,truth.default.over_3sd);
for(const year of ['2008','2017','2020','2022']){
  const rs=rows.filter(r=>r[0].startsWith(year));
  for(let j=1;j<=10;j++)for(const alpha of [.95,.975,.99]){
    const a=rs.map(r=>r[j]);assert(Number.isFinite(F.es(a,alpha)));
    assert(Number.isFinite(F.gaussianES(a,alpha)));
  }
}
console.log('Pass: Python/JS agreement, fractional ES, normal reference, tail counts, and all industry/year edge cases.');
for(const weight of [0,.23,.5,1]){
 const a=rows.map(r=>r[ia]),b=rows.map(r=>r[ib]);
 const tc=F.tailContributions(a,b,weight,.975);
 assert(Math.abs(tc.a+tc.b-F.es(a.map((v,i)=>weight*v+(1-weight)*b[i]),.975))<1e-10);
}
assert.equal(F.bestAllocation([3,4],[.1,.2],2),-1);
assert.equal(F.bestAllocation([3,4],[.1,.2],3.5),0);
assert.equal(F.bestAllocation([3,4],[.1,.2],4),1);
const allocation=read('allocation-results.json');
const training=rows.filter(r=>r[0]>='2010-01-01'&&r[0]<'2020-01-01');
for(const model of allocation.methods){
 const risks=[],means=[];
 for(let w=0;w<=100;w++){
  const rs=training.map(r=>w/100*r[ia]+(1-w/100)*r[ib]);
  means.push(F.mean(rs));risks.push(model.method==='Historical ES'?F.es(rs,.975):F.gaussianES(rs,.975));
 }
 assert.equal(F.bestAllocation(risks,means,2.5),Math.round(model.technology_weight*100));
}
console.log('Pass: tail contributions sum to ES, infeasible budgets, and independent constrained allocation results.');
const multi=read('portfolio-results.json');
const frontier=require('./portfolio.js').frontier;
const train10=rows.filter(r=>r[0]>='2010-01-01'&&r[0]<'2020-01-01'),test10=rows.filter(r=>r[0]>='2020-01-01');
for(const [alphaKey,study] of Object.entries(multi.levels)){
 const alpha=+alphaKey;
 for(const method of ['historical','gaussian']){
  let last=-Infinity;
  for(const [i,s] of study[method].entries()){
   if(!s)continue;
   const r=train10.map(row=>s.weights.reduce((sum,w,j)=>sum+w*row[j+1],0));
   const te=test10.map(row=>s.weights.reduce((sum,w,j)=>sum+w*row[j+1],0));
   assert(Math.abs(s.weights.reduce((a,b)=>a+b,0)-1)<1e-7);
   assert(s.weights.every(w=>w>=-1e-7));
   assert(Math.abs(F.mean(r)-s.mean)<1e-10);
   assert(Math.abs(F.es(r,alpha)-s.es)<1e-8);
   assert(Math.abs(F.gaussianES(r,alpha)-s.normal_es)<1e-7);
   assert(Math.abs(F.es(te,alpha)-s.test_es)<1e-8);
   assert(Math.abs(F.mean(te)-s.test_mean)<1e-10);
   assert(Math.abs(s.contributions.reduce((a,b)=>a+b,0)-s.es)<1e-9);
   assert(s[method==='historical'?'es':'normal_es']<=study.budgets[i]+1e-6);
   assert(s.mean>=last-1e-7);last=s.mean;
  }
  const f=frontier(study[method],method==='historical'?'es':'normal_es');
  assert(f.length>2);
  for(let i=1;i<f.length;i++)assert(f[i].mean>f[i-1].mean);
 }
 assert.equal(study.historical[0],null);
 assert.equal(study.gaussian[0],null);
}
console.log('Pass: all ten-industry allocations independently recomputed, feasibility, later-sample results, ES additivity, and efficient boundary.');
