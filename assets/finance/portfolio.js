(function(){
'use strict';
function frontier(points,key){let best=-Infinity;return points.filter(Boolean).sort((a,b)=>a[key]-b[key]||b.mean-a.mean).filter(p=>{if(p.mean<=best+1e-9)return false;best=p.mean;return true;});}
if(typeof module!=='undefined'){module.exports={frontier};return;}
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
 const l=layout(`${(alpha*100).toFixed(1)}% ES · daily loss (%)`,'Estimated mean daily return (%)');l.margin.t=window.innerWidth<620?110:60;l.legend.y=1.02;l.legend.yanchor='bottom';
 const traces=[];
 for(const [model,key,name,color,dash] of [['historical','es','Historical frontier',blue,'solid'],['gaussian','normal_es','Normal frontier',gray,'dash']]){
  const ps=frontier(S[model],key);traces.push({x:ps.map(p=>p[key]),y:ps.map(p=>p.mean),type:'scatter',mode:'lines',name,line:{color,dash,width:2},hovertemplate:'ES: %{x:.2f}%<br>Mean: %{y:.3f}%<extra>%{fullData.name}</extra>'});
 }
 traces.push({x:[S.equal.es],y:[S.equal.mean],type:'scatter',mode:'markers',name:'Equal weight',marker:{color:'#6a7884',symbol:'diamond',size:9}});
 if(selected){
   traces.push({x:[selected.es],y:[selected.mean],type:'scatter',mode:'markers',name:'Selected · observed ES',marker:{color:orange,size:11},hovertemplate:'Historical ES: %{x:.2f}%<br>Mean: %{y:.3f}%<extra>Selected allocation</extra>'});
   if(strategy==='gaussian'){
    traces.push({x:[selected.normal_es,selected.es],y:[selected.mean,selected.mean],type:'scatter',mode:'lines+markers',showlegend:false,marker:{color:orange,symbol:['circle-open','circle'],size:[11,0]},line:{color:orange,dash:'dot',width:1},hovertemplate:'ES: %{x:.2f}%<extra>Normal estimate → historical risk</extra>'});
   }
 }
 if(strategy!=='equal')l.shapes=[{type:'line',x0:b,x1:b,y0:0,y1:1,yref:'paper',line:{color:orange,dash:'dot',width:1}}];
 Plotly.react('frontier',traces,l,config);
 $('frontier-note').textContent=`Each curve is the estimated efficient boundary for its own risk model: higher mean at a given ES budget. The vertical line is the training budget. ${strategy==='gaussian'?'The hollow orange point uses normal ES; the filled point shows historical ES of the SAME weights. Their horizontal gap measures tail-risk underestimation for this allocation.':strategy==='equal'?'The diamond is the observed risk and mean of the 10%-each benchmark.':'The orange point shows the chosen allocation. Blue uses observed losses; dashed gray uses only means and covariances.'}`;
 $('results').replaceChildren();
 if(selected)for(const row of [['Training · 2010–2019',selected.mean,selected.es,selected.normal_es],['Evaluation · 2020–2025',selected.test_mean,selected.test_es,null]]){
  const tr=document.createElement('tr');for(const [i,v] of row.entries()){const td=document.createElement('td');td.textContent=i===0?v:v===null?'—':v.toFixed(i===1?3:2)+'%';tr.append(td);}$('results').append(tr);
 }
}
for(const id of ['alpha','strategy','budget'])$(id).addEventListener('input',render);
render();new ResizeObserver(()=>parent.postMessage({type:'financeLabResize',height:Math.ceil(document.body.getBoundingClientRect().height)+32},location.origin)).observe(document.body);
})();
