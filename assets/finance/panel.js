/* Historical diagnostics, using percent returns throughout. No analytics or API calls. */
(function(){
'use strict';
const mean=a=>a.reduce((x,y)=>x+y,0)/a.length;
const sd=a=>{const m=mean(a);return Math.sqrt(a.reduce((s,x)=>s+(x-m)**2,0)/(a.length-1));};
function quantile(a,p){const s=[...a].sort((x,y)=>x-y),h=(s.length-1)*p,i=Math.floor(h);return s[i]+(s[Math.min(i+1,s.length-1)]-s[i])*(h-i);}
function es(a,alpha){const l=a.map(x=>-x).sort((x,y)=>y-x),mass=l.length*(1-alpha),k=Math.floor(mass);return (l.slice(0,k).reduce((s,x)=>s+x,0)+(mass-k)*(l[k]||0))/mass;}
// Acklam inverse-normal approximation, adequate for plotting positions in this sample.
function invnorm(p){const a=[-39.6968302866538,220.946098424521,-275.928510446969,138.357751867269,-30.6647980661472,2.50662827745924],b=[-54.4760987982241,161.585836858041,-155.698979859887,66.8013118877197,-13.2806815528857],c=[-.00778489400243029,-.322396458041136,-2.40075827716184,-2.54973253934373,4.37466414146497,2.93816398269878],d=[.00778469570904146,.32246712907004,2.445134137143,3.75440866190742];let q,r;
if(p<.02425){q=Math.sqrt(-2*Math.log(p));return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);}
if(p>1-.02425)return -invnorm(1-p);q=p-.5;r=q*q;return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q/(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);}
function normalTail(x){const z=x/Math.sqrt(2),t=1/(1+.3275911*z);return (((((1.061405429*t-1.453152027)*t)+1.421413741)*t-.284496736)*t+.254829592)*t*Math.exp(-z*z);}
const gaussianES=(a,alpha)=>-mean(a)+sd(a)*Math.exp(-(invnorm(alpha)**2)/2)/Math.sqrt(2*Math.PI)/(1-alpha);
// Expose pure functions for reproducible cross-language checks.
if(typeof module!=='undefined')module.exports={mean,sd,quantile,es,invnorm,normalTail,gaussianES};
if(typeof document==='undefined')return;
const $=id=>document.getElementById(id),D=window.FINANCE_DATA;
if(!D||!window.Plotly){$('error').textContent='The charts could not load. Please check your connection to the chart library; the CSV remains downloadable.';return;}
const blue='#237c9c',orange='#be6535',gray='#667684';
for(const key of ['asset-a','asset-b'])for(const c of D.meta.columns){const o=document.createElement('option');o.value=c;o.textContent=D.meta.names[c];$(key).append(o);}
$('asset-a').value='HiTec';$('asset-b').value='Utils';
$('provenance').textContent=`Kenneth R. French Data Library · ${D.meta.n.toLocaleString()} trading days · ${D.meta.start} to ${D.meta.end} · value-weighted daily returns (%) · ${D.meta.source_header} Snapshot: ${D.meta.download_date}.`;
const config={responsive:true,displaylogo:false,modeBarButtonsToRemove:['select2d','lasso2d','autoScale2d'],toImageButtonOptions:{format:'svg',filename:'finance-tail-study'}};
const layout=(x,y)=>({margin:{l:60,r:14,t:16,b:55},font:{family:'system-ui, sans-serif',size:12,color:'#314252'},paper_bgcolor:'#fff',plot_bgcolor:'#fff',xaxis:{title:{text:x,font:{size:12}},gridcolor:'#e8edf1',zerolinecolor:'#b7c4ce',automargin:true},yaxis:{title:{text:y,font:{size:12}},gridcolor:'#e8edf1',zerolinecolor:'#b7c4ce',automargin:true},legend:{orientation:'h',x:0,y:1.14,font:{size:11}},hovermode:'closest'});
const line=(x,y,name,color,dash)=>({x,y,name,type:'scatter',mode:'lines',line:{color,width:2,dash:dash||'solid'},hovertemplate:'%{x:.2f}, %{y:.4g}<extra>%{fullData.name}</extra>'});
let curveCache=null;
function render(){
const A=$('asset-a').value,B=$('asset-b').value,ia=D.meta.columns.indexOf(A)+1,ib=D.meta.columns.indexOf(B)+1,p=$('period').value,w=+$('weight').value/100,alpha=+$('alpha').value,h=+$('threshold').value;
const rows=D.rows.filter(r=>p==='all'||(p==='train'?r[0]<'2020-01-01':p==='test'?r[0]>='2020-01-01':r[0].startsWith(p)));
const a=rows.map(r=>r[ia]),b=rows.map(r=>r[ib]),r=a.map((v,i)=>w*v+(1-w)*b[i]),mu=mean(r),s=sd(r),n=r.length,z=r.map(x=>(x-mu)/s),hits=z.filter(x=>Math.abs(x)>h).length,expected=n*normalTail(h),e=es(r,alpha),g=gaussianES(r,alpha);
$('weight-label').textContent=`${D.meta.names[A]}: ${Math.round(w*100)}% · B: ${Math.round((1-w)*100)}%`;
$('summary').textContent=`${n.toLocaleString()} days · daily SD ${s.toFixed(2)}% · |standardized return| > ${h}: ${hits} days (normal reference ${expected.toFixed(1)}) · ${(alpha*100).toFixed(1)}% historical ES ${e.toFixed(2)}% loss vs normal ${g.toFixed(2)}%.`;
const ordered=z.map((v,i)=>({v,date:rows[i][0]})).sort((x,y)=>x.v-y.v),qx=ordered.map((_,i)=>invnorm((i+.5)/n));
Plotly.react('qq',[line([-4,4],[-4,4],'Normal reference',gray,'dash'),{x:qx,y:ordered.map(o=>o.v),text:ordered.map(o=>o.date),type:'scatter',mode:'markers',name:'Observed',marker:{color:blue,size:3},hovertemplate:'%{text}<br>Normal quantile: %{x:.2f}<br>Observed z: %{y:.2f}<extra></extra>'}],layout('Normal quantile','Observed standardized return'),config);
const maxz=Math.max(...z.map(Math.abs)),tx=Array.from({length:81},(_,i)=>maxz*i/80),ty=tx.map(t=>z.filter(v=>Math.abs(v)>t).length/n),nx=Array.from({length:101},(_,i)=>i/20),ny=nx.map(normalTail);
const tl=layout('Magnitude threshold (SD units)','P(|standardized return| > x)');tl.yaxis.type='log';tl.yaxis.range=[Math.log10(1/(n*3)),0];tl.yaxis.tickvals=[1,.1,.01,.001,.0001,.00001];tl.yaxis.ticktext=['100%','10%','1%','0.1%','0.01%','0.001%'];tl.shapes=[{type:'line',x0:h,x1:h,y0:1/(n*3),y1:1,xref:'x',yref:'y',line:{color:orange,width:1,dash:'dot'}}];
Plotly.react('tail',[line(tx,ty.map(v=>v>0?v:null),'Observed',blue),line(nx,ny,'Normal',gray,'dash')],tl,config);
$('tail-comment').textContent=`Logarithmic y-axis. At ${h} SD: observed ${(hits/n*100).toFixed(2)}%, normal ${(normalTail(h)*100).toFixed(3)}%. Zero observed exceedances are omitted, rather than plotted as proof of zero risk.`;
const qa=quantile(a,.05),qb=quantile(b,.05),joint=[],other=[];
for(let i=0;i<n;i++)(a[i]<=qa&&b[i]<=qb?joint:other).push(i);
const ma=mean(a),mb=mean(b),corr=a.reduce((v,x,i)=>v+(x-ma)*(b[i]-mb),0)/((n-1)*sd(a)*sd(b));
const sc=(idx,name,color,size)=>({x:idx.map(i=>a[i]),y:idx.map(i=>b[i]),text:idx.map(i=>rows[i][0]),type:'scatter',mode:'markers',name,marker:{color,size,opacity:.65},hovertemplate:'%{text}<br>A: %{x:.2f}%<br>B: %{y:.2f}%<extra>%{fullData.name}</extra>'});
const jl=layout(D.meta.names[A]+' return (%)',D.meta.names[B]+' return (%)');jl.shapes=[{type:'line',x0:qa,x1:qa,y0:0,y1:1,xref:'x',yref:'paper',line:{color:orange,width:1,dash:'dot'}},{type:'line',y0:qb,y1:qb,x0:0,x1:1,xref:'paper',yref:'y',line:{color:orange,width:1,dash:'dot'}}];
Plotly.react('joint',[sc(other,'Other days',blue,3),sc(joint,'Both ≤ 5th percentile',orange,5)],jl,config);
$('joint-comment').textContent=`Correlation ${corr.toFixed(2)}. Both industries are in their own bottom 5% on ${(100*joint.length/n).toFixed(2)}% of days (independence benchmark ≈ 0.25%). This finite-threshold frequency is not an asymptotic tail-dependence estimate.`;
const key=[A,B,p,alpha].join('|');
if(!curveCache||curveCache.key!==key){const ws=Array.from({length:101},(_,i)=>i),hist=[],norm=[];for(const v of ws){const vrs=a.map((x,i)=>v/100*x+(1-v/100)*b[i]);hist.push(es(vrs,alpha));norm.push(gaussianES(vrs,alpha));}curveCache={key,ws,hist,norm};}
const {ws,hist,norm}=curveCache,rl=layout('Weight in industry A (%)',`${(100*alpha).toFixed(1)}% ES · daily loss (%)`);rl.shapes=[{type:'line',x0:w*100,x1:w*100,y0:0,y1:1,xref:'x',yref:'paper',line:{color:orange,width:1,dash:'dot'}}];
Plotly.react('risk',[line(ws,hist,'Historical ES',blue),line(ws,norm,'Normal ES',gray,'dash'),{x:[w*100],y:[e],type:'scatter',mode:'markers',name:'Your allocation',marker:{color:orange,size:9}}],rl,config);
const min=Math.min(...hist),wi=hist.indexOf(min);
$('risk-comment').textContent=`Selected-period minimum on a 1% weight grid: ${wi}% A, ES ${min.toFixed(2)}%. Tail mass: ${(n*(1-alpha)).toFixed(1)} observations. This is an in-sample comparison, not a future-optimal allocation.`;
const worst=r.map((v,i)=>({v,date:rows[i][0]})).sort((x,y)=>x.v-y.v).slice(0,8);
$('extremes').replaceChildren(...worst.map(o=>{const tr=document.createElement('tr');for(const v of [o.date,o.v.toFixed(2)+'%',(-o.v).toFixed(2)+'%']){const td=document.createElement('td');td.textContent=v;tr.append(td);}return tr;}));
}
let timer;
for(const id of ['asset-a','asset-b','period','weight','threshold','alpha'])$(id).addEventListener('input',()=>{clearTimeout(timer);timer=setTimeout(render,60);});
render();
new ResizeObserver(()=>parent.postMessage({type:'financeLabResize',height:Math.ceil(document.body.getBoundingClientRect().height)+32},location.origin)).observe(document.body);
})();
