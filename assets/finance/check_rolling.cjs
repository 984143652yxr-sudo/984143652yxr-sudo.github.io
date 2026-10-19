'use strict';
const fs=require('fs'),assert=require('assert'),crypto=require('crypto');
const root=__dirname,D=JSON.parse(fs.readFileSync(root+'/rolling-results.json')),raw=fs.readFileSync(root+'/returns.json'),input=JSON.parse(raw);
assert.equal(D.source_sha256,crypto.createHash('sha256').update(raw).digest('hex'));
const close=(a,b,tol=1e-7)=>assert(Math.abs(a-b)<tol,`${a} != ${b}`);
const mean=a=>a.reduce((s,x)=>s+x,0)/a.length;
function ES(r){const loss=r.map(x=>-x).sort((a,b)=>b-a),mass=r.length*.025,k=Math.floor(mass);return (loss.slice(0,k).reduce((s,x)=>s+x,0)+(mass-k)*(loss[k]||0))/mass;}
const all={};assert.equal(D.months.length,72);assert.equal(D.months[0].train_start,'2017-01-03');
for(const month of D.months){const year=+month.month.slice(0,4),start=String(year-3)+month.month.slice(4)+'-01',cutoff=month.month+'-01';
 const train=input.rows.filter(r=>r[0]>=start&&r[0]<cutoff),next=input.rows.filter(r=>r[0].startsWith(month.month));assert.equal(train.length,month.train_n);assert.equal(train.at(-1)[0],month.train_end);assert(month.train_end<next[0][0]);assert.equal(next.length,month.dates.length);
 month.records.forEach((record,i)=>{for(const [key,p] of Object.entries(record)){
  const w=p.weights;close(w.reduce((s,x)=>s+x,0),1);assert(w.every(x=>x>=-1e-9));
  const r=train.map(r=>w.slice(0,10).reduce((s,x,j)=>s+x*r[j+1],0)),mu=mean(r),sd=Math.sqrt(r.reduce((s,x)=>s+(x-mu)**2,0)/(r.length-1));close(mu,p.forecast.mean);close(sd,p.forecast.sd);close(ES(r),p.forecast.es);close(-mu+D.k*sd,p.forecast.normal_es);
  const b=D.curves[i].budget;if(key==='historical')assert(p.forecast.es<=b+1e-6);if(key==='normal')assert(p.forecast.normal_es<=b+1e-6);if(key==='markowitz')assert(sd<=D.curves[i].volatility_budget+1e-6);
  // Independently evolve dollar holdings; no within-month weight reset.
  let holdings=w.slice(),prev=1;const realized=[];
  for(const row of next){holdings=holdings.map((h,j)=>h*(1+(j<10?row[j+1]:0)/100));const value=holdings.reduce((s,x)=>s+x,0);realized.push(100*(value/prev-1));prev=value;}
  realized.forEach((x,j)=>close(x,p.returns[j]));if(key==='oracle')assert(ES(realized)<=b+1e-5);
  (all[i+key]??=[]).push(...realized);
 }});
}
for(const [i,c] of D.curves.entries())for(const key of ['historical','normal','markowitz','oracle']){const p=c[key],r=all[i+key];close(mean(r),p.mean);close(ES(r),p.es);assert.equal(r.length,D.dates.length-1);let v=1,peak=1,mdd=0;close(p.wealth[0],1);
 r.forEach((x,j)=>{v*=1+x/100;peak=Math.max(peak,v);const dd=1-v/peak;mdd=Math.max(mdd,dd);close(v,p.wealth[j+1]);close(dd,p.drawdown[j+1]);});close(v,p.terminal);close(mdd,p.max_drawdown);close(v**(1/6)-1,p.cagr);
}
// Hand-check monthly hold vs daily rebalancing and fractional-tail ES.
close(ES([1,-3,-2,-4]),4);
assert(Math.abs(.5*1.1*.9+.5*1 - (1.05*.95))>1e-4);
console.log('Pass: all 72 rolling windows, no training overlap with next month, budgets, drifted holdings, realised ES, wealth, drawdown, and oracle constraints.');
