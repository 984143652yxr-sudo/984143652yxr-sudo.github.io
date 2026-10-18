(function(){
'use strict';
function frontier(points,key){let best=-Infinity;return points.filter(Boolean).sort((a,b)=>a[key]-b[key]||b.mean-a.mean).filter(p=>{if(p.mean<=best+1e-9)return false;best=p.mean;return true;});}
function gaps(p,b){return {underestimate:p.es-p.normal_es,budgetExcess:p.es-b};}
if(typeof module!=='undefined'){module.exports={frontier,gaps};return;}
const D=window.PORTFOLIO_STUDY,$=id=>document.getElementById(id);
if(!D||!window.Plotly){$('error').textContent='Charts could not load. Downloadable calculations remain available below.';return;}
const blue='#237c9c',gray='#83939c',orange='#b96937';
const config={responsive:true,displaylogo:false,toImageButtonOptions:{format:'svg',filename:'ten-industry-portfolio'}};
const layout=(x,y)=>({margin:{l:60,r:15,t:35,b:55},font:{family:'system-ui',size:12,color:'#314252'},xaxis:{title:x,gridcolor:'#e8edf1',automargin:true},yaxis:{title:y,gridcolor:'#e8edf1',automargin:true},legend:{orientation:'h',x:0,y:1.18,font:{size:11}},hovermode:'closest'});
const labels=['Non-durables','Durables','Manufacturing','Energy','Technology','Telecom','Shops','Healthcare','Utilities','Other'];
function render(){
 const alpha=+$('alpha').value,S=D.levels[String(alpha)],index=+$('budget').value,b=S.budgets[index],strategy=$('strategy').value;
 const selected=strategy==='equal'?S.equal:S[strategy][index];
 $('budget-label').textContent=`Training ES budget: ${b.toFixed(2)}% loss`;
 $('budget').disabled=strategy==='equal';
 $('summary').textContent=selected?`${strategy==='equal'?'Equal weight':'Optimized weights'} · estimated mean ${selected.mean.toFixed(3)}% / day · historical ES ${selected.es.toFixed(2)}% · normal ES ${selected.normal_es.toFixed(2)}%.`:`No feasible ${strategy==='historical'?'historical':'normal'}-ES allocation at this budget. Increase the budget. Minimum estimated ES: ${(strategy==='historical'?S.minimum_historical_es:S.minimum_gaussian_es).toFixed(2)}%.`;
 for(const [id,key,mult,color,title] of [['weights','weights',100,blue,'Capital weight (%)'],['contributions','contributions',1,orange,'Contribution to ES (percentage points)']]){
   const l=layout(title,'');l.margin.l=105;l.margin.t=15;l.yaxis.autorange='reversed';
   Plotly.react(id,selected?[{x:selected[key].map(v=>v*mult),y:labels,type:'bar',orientation:'h',marker:{color},hovertemplate:'%{y}: %{x:.3f}<extra></extra>'}]:[],l,config);
 }
 $('contribution-note').textContent=selected?`Left: where capital is invested. Right: losses from each industry on this portfolio's worst ${((1-alpha)*100).toFixed(1)}% of training days (fractional boundary). Contributions sum to ${selected.es.toFixed(2)}% ES; a negative bar would offset losses.`:'No weights or loss contributions exist for this infeasible choice.';
 const label=strategy==='equal'?'w<sup>*</sup> = w<sup>EW</sup> = (0.1, …, 0.1)':strategy==='historical'?'w<sup>*</sup> = w<sub>H</sub>(b)':'w<sup>*</sup> = w<sub>N</sub>(b)';
 $('selection-equation').innerHTML=`Selected weights: ${label}; α = ${(100*alpha).toFixed(1)}%; b = ${b.toFixed(2)}% daily tail loss.`;
 const l=layout(`${(alpha*100).toFixed(1)}% ES · daily loss (%)`,'m(w) · estimated mean daily return (%)');
 l.margin.t=50;l.margin.r=28;l.showlegend=false;
 const traces=[];
 for(const [model,key,name,color,dash] of [['historical','es','ΓH · historical frontier',blue,'solid'],['gaussian','normal_es','ΓN · normal frontier',gray,'dash'],['gaussian','es','ΓN→H · normal weights, historical ES',orange,'dot']]){
  const ps=name.startsWith('ΓN→H')?S[model].filter(Boolean):frontier(S[model],key);
  traces.push({x:ps.map(p=>p[key]),y:ps.map(p=>p.mean),type:'scatter',mode:'lines',name,line:{color,dash,width:2},hovertemplate:'ES: %{x:.2f}%<br>Mean: %{y:.3f}%<extra>%{fullData.name}</extra>'});
 }
 traces.push({x:[S.equal.es],y:[S.equal.mean],type:'scatter',mode:'markers',name:'Equal weight · Hα(wᴱᵂ)',marker:{color:'#6a7884',symbol:'diamond',size:9},hovertemplate:'Historical ES: %{x:.2f}%<br>Mean: %{y:.3f}%<extra>10% in each industry</extra>'});
 l.annotations=[{xref:'paper',yref:'paper',x:0.02,y:1.06,xanchor:'left',yanchor:'bottom',text:'← Lower ES · ↑ Higher mean',showarrow:false,font:{color:'#386d56',size:11},bgcolor:'rgba(255,255,255,0.9)'}];
 if(selected){
   for(const [key,letter,symbol,desc] of [['normal_es','A','circle-open','A: Nα(w*) · normal estimate'],['es','B','circle','B: Hα(w*) · historical estimate']]){
     traces.push({x:[selected[key]],y:[selected.mean],type:'scatter',mode:'markers+text',text:[letter],textposition:'top center',textfont:{color:orange,size:14},name:desc,marker:{color:orange,symbol,size:12,line:{width:2}},hovertemplate:'ES: %{x:.2f}%<br>Mean: %{y:.3f}%<extra>%{fullData.name}</extra>'});
   }
   const g=gaps(selected,b),sign=g.underestimate>=0?'+':'−';
   l.annotations.push({x:selected.es,y:selected.mean,ax:selected.normal_es,ay:selected.mean,xref:'x',yref:'y',axref:'x',ayref:'y',text:'',showarrow:true,arrowcolor:orange,arrowwidth:1.5,arrowhead:3,standoff:8,startstandoff:8});
   l.annotations.push({x:(selected.es+selected.normal_es)/2,y:selected.mean,xref:'x',yref:'y',text:`H − N = ${sign}${Math.abs(g.underestimate).toFixed(2)} pp`,showarrow:false,yshift:-22,font:{color:orange,size:11},bgcolor:'rgba(255,255,255,0.9)'});
   $('frontier-note').textContent=`A → B: normal ES ${selected.normal_es.toFixed(2)}% → historical ES ${selected.es.toFixed(2)}%, with identical weights and mean return. ${g.underestimate>=0?'Underestimation':'Overestimation'}: ${Math.abs(g.underestimate).toFixed(2)} percentage points. ${g.budgetExcess>1e-6?`B lies RIGHT of b: historical tail risk exceeds the budget by ${g.budgetExcess.toFixed(2)} pp.`:`B lies on or LEFT of b: historical tail risk meets this training budget.`}`;
 }else{
   $('frontier-note').textContent='No selected A/B points: this risk model has no feasible allocation at the chosen budget. Raise b or choose equal weights.';
 }
 l.shapes=[{type:'line',x0:b,x1:b,y0:0,y1:1,yref:'paper',line:{color:'#596775',dash:'dot',width:1}}];
 l.annotations.push({x:b,y:0,xref:'x',yref:'paper',text:'b',showarrow:false,yshift:-18,font:{color:'#596775',size:13}});
 Plotly.react('frontier',traces,l,config);
 const normal=S.gaussian[index],historical=S.historical[index];
 const signed=v=>Math.abs(v)<.005?'0.00 pp':`${v>=0?'+':'−'}${Math.abs(v).toFixed(2)} pp`;
 const entries=[
  ['TRAIN · estimated mean m(w)',p=>p.mean.toFixed(3)+'%'],
  ['Normal ES Nα(w)',p=>p.normal_es.toFixed(2)+'%'],
  ['Historical ES Hα(w)',p=>p.es.toFixed(2)+'%'],
  ['Risk underestimate Hα(w) − Nα(w)',p=>signed(gaps(p,b).underestimate)],
  ['Budget excess Hα(w) − b',p=>signed(gaps(p,b).budgetExcess)],
  ['LATER · historical ES (2020–2025)',p=>p.test_es.toFixed(2)+'%'],
  ['Later mean daily return',p=>p.test_mean.toFixed(3)+'%']
 ];
 $('results').replaceChildren(...entries.map(([name,format],i)=>{const tr=document.createElement('tr');if(i===5)tr.className='separator';const th=document.createElement('th');th.scope='row';th.textContent=name;tr.append(th);for(const p of [normal,historical]){const td=document.createElement('td');td.textContent=p?format(p):'Infeasible';tr.append(td);}return tr;}));
 if(normal&&historical){
   const d=normal.test_es-historical.test_es,dm=normal.test_mean-historical.test_mean;
   $('comparison-note').textContent=`Same budget b = ${b.toFixed(2)}%: the normal-chosen allocation has historical training ES ${normal.es.toFixed(2)}%, versus ${historical.es.toFixed(2)}% for the historical-chosen allocation. In 2020–2025, normal-chosen ES is ${normal.test_es.toFixed(2)}% versus ${historical.test_es.toFixed(2)}%: ${Math.abs(d).toFixed(2)} pp ${d>=0?'more':'less'} average tail loss. Its mean daily return is ${Math.abs(dm).toFixed(3)} pp ${dm>=0?'higher':'lower'} (${normal.test_mean.toFixed(3)}% versus ${historical.test_mean.toFixed(3)}%).`;
 }else $('comparison-note').textContent='A same-budget comparison requires both models to have feasible allocations. Raise b to compare the two columns.';

}
for(const id of ['alpha','strategy','budget'])$(id).addEventListener('input',render);
render();new ResizeObserver(()=>parent.postMessage({type:'financeLabResize',height:Math.ceil(document.body.getBoundingClientRect().height)+32},location.origin)).observe(document.body);
})();
