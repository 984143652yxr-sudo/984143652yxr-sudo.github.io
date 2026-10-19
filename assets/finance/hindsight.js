(function(){
'use strict';
const D=window.HINDSIGHT_STUDY,$=id=>document.getElementById(id);
if(!D||!window.Plotly){$('error').textContent='The chart could not load. Calculations remain available in the downloads.';return;}
const blue='#237c9c',purple='#8064a2',green='#36816b';
const cfg={responsive:true,displaylogo:false,toImageButtonOptions:{format:'svg',filename:'markowitz-es-hindsight'}};
function render(){
 const alpha=+$('alpha').value,S=D.levels[String(alpha)],i=+$('target').value,r=S.records[i],isTest=$('view').value==='test',x=isTest?'test_es':'es',y=isTest?'test_mean':'mean',period=isTest?'2020–2025':'2010–2019';
 $('target-label').textContent=`Training mean target q = ${r.target.toFixed(4)}% / day`;
 const traces=[];
 for(const [key,name,color,dash] of [['historical','ES-trained',blue,'solid'],['markowitz','Markowitz',purple,'dash']]){
  const points=S.records.map(r=>r[key]);traces.push({x:points.map(p=>p[x]),y:points.map(p=>p[y]),customdata:S.records.map(r=>r.target),mode:'lines',type:'scatter',name,line:{color,dash,width:2},hovertemplate:'Training target q: %{customdata:.4f}%<br>ES: %{x:.3f}%<br>Mean: %{y:.4f}%<extra>%{fullData.name}</extra>'});
  const p=r[key];traces.push({x:[p[x]],y:[p[y]],mode:'markers+text',type:'scatter',text:[key==='historical'?'E':'M'],textposition:key==='historical'?'top center':'bottom center',textfont:{color},marker:{color,size:11},name:`Selected ${name}`,hovertemplate:'ES: %{x:.3f}%<br>Mean: %{y:.4f}%<extra>%{fullData.name}</extra>'});
 }
 if(isTest){
  traces.unshift({x:S.oracle.map(p=>p.test_es),y:S.oracle.map(p=>p.test_mean),type:'scatter',mode:'lines',line:{color:green,dash:'dot',width:2.5},name:'Hindsight · uses future data',hovertemplate:'ES: %{x:.3f}%<br>Mean: %{y:.4f}%<extra>Hindsight</extra>'});
  for(const key of ['historical','markowitz']){const p=r[key],o=p.oracle;traces.push({x:[p.test_es,o.test_es],y:[p.test_mean,o.test_mean],type:'scatter',mode:'lines+markers',line:{color:key==='historical'?blue:purple,dash:'dot',width:1},marker:{symbol:'circle-open',size:8,color:green},name:'Matched hindsight comparison',hovertemplate:'ES: %{x:.3f}%<br>Mean: %{y:.4f}%<extra>Hindsight at ≥ realised mean</extra>'});}
 }
 traces.push({x:[S.equal[x]],y:[S.equal[y]],mode:'markers',type:'scatter',name:'Equal weight',marker:{symbol:'diamond',color:'#83939c',size:9}});
 const l={showlegend:false,margin:{l:65,r:20,t:20,b:65},font:{family:'system-ui',size:12,color:'#314252'},xaxis:{title:`${period} · ${(100*alpha).toFixed(1)}% historical ES (%)`,gridcolor:'#e8edf1',automargin:true},yaxis:{title:`${period} · mean daily return (%)`,gridcolor:'#e8edf1',automargin:true}};
 if($('range').value==='compare'){const ps=S.records.flatMap(r=>[r.historical,r.markowitz]);if(isTest)ps.push(...S.records.flatMap(r=>[r.historical.oracle,r.markowitz.oracle]));ps.push(S.equal);const xs=ps.map(p=>p[x]),ys=ps.map(p=>p[y]);const xmin=Math.min(...xs),xmax=Math.max(...xs),ymin=Math.min(...ys),ymax=Math.max(...ys);l.xaxis.range=[xmin-(xmax-xmin)*.06,xmax+(xmax-xmin)*.08];l.yaxis.range=[ymin-(ymax-ymin)*.1,ymax+(ymax-ymin)*.1];}
 Plotly.react('performance',traces,l,cfg);
 const h=r.historical,m=r.markowitz,d=m.test_es-h.test_es,dm=m.test_mean-h.test_mean;
 $('reading').textContent=isTest?`At the same training target, Markowitz later ES is ${m.test_es.toFixed(2)}%, versus ES-trained ${h.test_es.toFixed(2)}%: ${Math.abs(d)<.005?'essentially equal':`${Math.abs(d).toFixed(2)} pp ${d>0?'higher':'lower'}`}. Markowitz later mean is ${Math.abs(dm)<.00005?'essentially equal':`${Math.abs(dm).toFixed(4)} pp ${dm>0?'higher':'lower'}`}. Open green points are hindsight optima matching at least each selected portfolio's realised mean.`:`Training check: ES-trained ES ${h.es.toFixed(3)}% ≤ Markowitz ES ${m.es.toFixed(3)}%; Markowitz SD ${m.sd.toFixed(3)}% ≤ ES-trained SD ${h.sd.toFixed(3)}%. Each criterion wins on its own training objective. Green is reserved for the evaluation-period hindsight comparison.`;
 const rows=[['Training mean',p=>p.mean],['Training SD',p=>p.sd],['Training ES',p=>p.es],['Evaluation mean',p=>p.test_mean],['Evaluation SD',p=>p.test_sd],['Evaluation ES',p=>p.test_es],['Hindsight ES at ≥ this evaluation mean',p=>p.oracle.test_es],['ES above matched hindsight (pp)',p=>Math.max(0,p.hindsight_gap)]];
 $('comparison').replaceChildren(...rows.map(([name,f],j)=>{const tr=document.createElement('tr');for(const v of [name,f(h).toFixed(3)+(j===7?' pp':'%'),f(m).toFixed(3)+(j===7?' pp':'%')]){const td=document.createElement('td');td.textContent=v;tr.append(td);}return tr;}));
 $('oracle-note').textContent=`Hindsight uses all evaluation days: for E, mean ${h.oracle.test_mean.toFixed(4)}% and ES ${h.oracle.test_es.toFixed(3)}%; for M, mean ${m.oracle.test_mean.toFixed(4)}% and ES ${m.oracle.test_es.toFixed(3)}%. Positive gaps measure room for hindsight improvement, not a guaranteed future gain.`;
 $('weights-table').replaceChildren(...D.columns.map((c,j)=>{const tr=document.createElement('tr');for(const v of [D.names[c],... [h,m,h.oracle,m.oracle].map(p=>(p.weights[j]*100).toFixed(2)+'%')]){const td=document.createElement('td');td.textContent=v;tr.append(td);}return tr;}));
}
for(const id of ['alpha','view','target','range'])$(id).addEventListener('input',render);
render();new ResizeObserver(()=>parent.postMessage({type:'financeLabResize',height:Math.ceil(document.body.getBoundingClientRect().height)+32},location.origin)).observe(document.body);
})();
