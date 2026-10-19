(function(){
'use strict';
const K={'0.95':2.0627128075074253,'0.975':2.3378027922014133,'0.99':2.665214220345805};
function normal(p,a){return -p.mean+K[String(a)]*p.sd;}
function contributions(rows,w,a){const sorted=rows.map(r=>({r,loss:-r.reduce((s,v,j)=>s+v*w[j],0)})).sort((a,b)=>b.loss-a.loss),mass=rows.length*(1-a);return w.map((v,j)=>sorted.reduce((s,o,i)=>s-Math.max(0,Math.min(1,mass-i))*v*o.r[j],0)/mass);}
function diagnostics(e,m,a){return {normalGap:m.es-normal(m,a),trainingBenefit:m.es-e.es,calibration:e.test_es-e.es,hindsight:e.test_es-e.oracle.test_es};}
if(typeof module!=='undefined'&&module.exports){module.exports={normal,diagnostics,contributions};return;}
const D=window.HINDSIGHT_STUDY,$=id=>document.getElementById(id),blue='#237c9c',orange='#b96937',green='#36816b',gray='#83939c';
if(!D||!window.Plotly){$('error').textContent='Unable to load the chart.';return;}
const f=(x,n=3)=>x.toFixed(n),cfg={responsive:true,displaylogo:false};
function line(ps,x,y,name,color,dash='solid'){return {x:ps.map(x),y:ps.map(y),name,type:'scatter',mode:'lines',line:{color,dash,width:2.5},hovertemplate:'ES: %{x:.3f}%<br>Mean: %{y:.4f}%<extra>%{fullData.name}</extra>'};}
function point(x,y,label,color,open=false){return {x:[x],y:[y],text:[label],name:label,type:'scatter',mode:'markers+text',textposition:label==='M'?'bottom center':'top center',marker:{color,size:10,symbol:open?'circle-open':'circle'},hovertemplate:'ES: %{x:.3f}%<br>Mean: %{y:.4f}%<extra>%{fullData.name}</extra>'};}
function layout(period){return {showlegend:false,margin:{l:65,r:20,t:30,b:65},font:{family:'system-ui',size:12},xaxis:{title:period+' · ES (%)',gridcolor:'#e8edf1'},yaxis:{title:period+' · mean return (%)',gridcolor:'#e8edf1'}};}
function render(){
const a=+$('alpha').value,S=D.levels[String(a)],r=S.records[+$('target').value],e=r.historical,m=r.markowitz,d=diagnostics(e,m,a),es=S.records.map(r=>r.historical),ms=S.records.map(r=>r.markowitz);
$('target-label').textContent='Training target q = '+f(r.target,4)+'% / day';
const train=[line(ms,p=>normal(p,a),p=>p.mean,'Normal estimate of M',gray,'dash'),line(ms,p=>p.es,p=>p.mean,'Markowitz · historical ES',orange,'dash'),line(es,p=>p.es,p=>p.mean,'CVaR · historical ES',blue),point(normal(m,a),m.mean,'A',gray,true),point(m.es,m.mean,'M',orange),point(e.es,e.mean,'E',blue)];
const l=layout('2010–2019');l.annotations=[{x:m.es,y:m.mean,ax:normal(m,a),ay:m.mean,xref:'x',yref:'y',axref:'x',ayref:'y',text:'',showarrow:true,arrowhead:2,arrowcolor:gray}];Plotly.react('frontier',train,l,cfg);
const later=[line(S.oracle,p=>p.test_es,p=>p.test_mean,'Hindsight',green,'dot'),line(ms,p=>p.test_es,p=>p.test_mean,'Markowitz',orange,'dash'),line(es,p=>p.test_es,p=>p.test_mean,'CVaR',blue)];
for(const [p,n,c] of [[e,'E',blue],[m,'M',orange]]){later.push(line([p,p.oracle],p=>p.test_es,p=>p.test_mean,'Matched hindsight',c,'dot'),point(p.test_es,p.test_mean,n,c),point(p.oracle.test_es,p.oracle.test_mean,'O'+n,green,true));}
const ll=layout('2020–2025'),ps=[...es,...ms,...es.map(p=>p.oracle),...ms.map(p=>p.oracle)];for(const [axis,key] of [['xaxis','test_es'],['yaxis','test_mean']]){const lo=Math.min(...ps.map(p=>p[key])),hi=Math.max(...ps.map(p=>p[key]));ll[axis].range=[lo-(hi-lo)*.12,hi+(hi-lo)*.12];}Plotly.react('evaluation',later,ll,cfg);
const selected=r[$('strategy').value];Plotly.react('weights',[{x:D.columns,y:selected.weights.map(w=>w*100),type:'bar',marker:{color:$('strategy').value==='historical'?blue:orange},hovertemplate:'%{x}: %{y:.2f}%<extra></extra>'}],{margin:{l:55,r:10,t:15,b:45},yaxis:{title:'Capital weight (%)'}},cfg);
const trainRows=window.FINANCE_DATA.rows.filter(r=>r[0]>='2010-01-01'&&r[0]<'2020-01-01').map(r=>r.slice(1));Plotly.react('contributions',[{x:D.columns,y:contributions(trainRows,selected.weights,a),type:'bar',marker:{color:$('strategy').value==='historical'?blue:orange},hovertemplate:'%{x}: %{y:.3f} pp<extra></extra>'}],{margin:{l:55,r:10,t:15,b:45},yaxis:{title:'Tail contribution (pp)'}},cfg);
$('normal-note').textContent=`1 · Normal risk: M reports ${f(normal(m,a))}% but its training tail average is ${f(m.es)}%. The ${f(d.normalGap)} pp rightward gap measures understated risk for identical weights. Gray's apparent advantage uses a different risk estimate.`;
$('benefit-note').textContent=`2 · CVaR versus variance: at q = ${f(r.target,4)}%, training ES falls from ${f(m.es)}% (M) to ${f(e.es)}% (E), a ${f(d.trainingBenefit)} pp reduction. Later ES is ${f(m.test_es)}% (M) versus ${f(e.test_es)}% (E); later means are ${f(m.test_mean,4)}% and ${f(e.test_mean,4)}%. Compare both risk and return.`;
$('calibration-note').textContent=`3 · Estimate versus later outcomes: E's training ES is ${f(e.es)}%; later ES is ${f(e.test_es)}% (change ${f(d.calibration)} pp). Matched hindsight achieves ${f(e.oracle.test_es)}%, leaving ${f(d.hindsight)} pp allocation opportunity. The forecast change and hindsight gap answer different questions.`;
const rows=[['Training normal ES N',p=>normal(p,a)],['Training historical ES H',p=>p.es],['Later realised ES H′',p=>p.test_es],['Later mean / day',p=>p.test_mean],['Hindsight ES at ≥ own later mean',p=>p.oracle.test_es]];$('results').replaceChildren(...rows.map(([n,fn])=>{const tr=document.createElement('tr');for(const t of [n,f(fn(e),4)+'%',f(fn(m),4)+'%']){const td=document.createElement('td');td.textContent=t;tr.append(td);}return tr;}));
}
for(const id of ['alpha','target','strategy'])$(id).addEventListener('input',render);render();new ResizeObserver(()=>parent.postMessage({type:'financeLabResize',height:Math.ceil(document.body.getBoundingClientRect().height)+32},location.origin)).observe(document.body);
})();