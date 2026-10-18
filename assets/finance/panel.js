(function(){
'use strict';
if(typeof module!=='undefined'){module.exports=require('./math.js');return;}
const F=window.FinanceMath,D=window.FINANCE_DATA,$=id=>document.getElementById(id);
if(!F||!D||!window.Plotly){$('error').textContent='Charts could not load. The data download remains available below.';return;}
const {mean,sd,invnorm,normalTail}=F;
for(const c of D.meta.columns){const o=document.createElement('option');o.value=c;o.textContent=D.meta.names[c];$('asset').append(o);}
$('asset').value='HiTec';
const blue='#237c9c',gray='#697c89';
const cfg={responsive:true,displaylogo:false,toImageButtonOptions:{format:'svg',filename:'heavy-tail-diagnostic'}};
const layout=(x,y)=>({margin:{l:60,r:12,t:35,b:60},font:{family:'system-ui',size:12,color:'#314252'},xaxis:{title:x,automargin:true,gridcolor:'#e8edf1'},yaxis:{title:y,automargin:true,gridcolor:'#e8edf1'},legend:{orientation:'h',x:0,y:1.16,font:{size:11}}});
function render(){
 const j=D.meta.columns.indexOf($('asset').value)+1,p=$('period').value,h=+$('threshold').value;
 const rows=D.rows.filter(r=>p==='all'||r[0].startsWith(p)),r=rows.map(r=>r[j]),m=mean(r),s=sd(r),z=r.map(x=>(x-m)/s),n=r.length;
 const hits=z.filter(v=>Math.abs(v)>h).length,expected=n*normalTail(h);
 $('summary').textContent=`${n.toLocaleString()} days · daily SD ${s.toFixed(2)}% · moves beyond ±${h} SD: ${hits} observed vs ${expected.toFixed(1)} under normality · ${(hits/expected).toFixed(1)}× the normal frequency.`;
 const ordered=z.map((v,i)=>({v,date:rows[i][0]})).sort((a,b)=>a.v-b.v);
 Plotly.react('qq',[{x:[-4,4],y:[-4,4],type:'scatter',mode:'lines',name:'Normal',line:{color:gray,dash:'dash'}},{x:ordered.map((_,i)=>invnorm((i+.5)/n)),y:ordered.map(o=>o.v),text:ordered.map(o=>o.date),type:'scatter',mode:'markers',name:'Observed',marker:{color:blue,size:3},hovertemplate:'%{text}<br>Normal: %{x:.2f}<br>Observed: %{y:.2f}<extra></extra>'}],layout('Normal quantile','Observed standardized return'),cfg);
 const tx=Array.from({length:101},(_,i)=>Math.max(...z.map(Math.abs))*i/100),nx=Array.from({length:101},(_,i)=>i/20),tl=layout('Threshold x (SD units)','P(|Z| > x)');
 tl.yaxis.type='log';tl.yaxis.range=[Math.log10(1/(3*n)),0];tl.yaxis.tickvals=[1,.1,.01,.001,.0001];tl.yaxis.ticktext=['100%','10%','1%','0.1%','0.01%'];
 tl.shapes=[{type:'line',x0:h,x1:h,y0:0,y1:1,yref:'paper',line:{color:'#bd7345',dash:'dot',width:1}}];
 Plotly.react('tail',[{x:tx,y:tx.map(t=>z.filter(v=>Math.abs(v)>t).length/n||null),mode:'lines',name:'Observed',line:{color:blue}},{x:nx,y:nx.map(normalTail),mode:'lines',name:'Normal',line:{color:gray,dash:'dash'}}],tl,cfg);
 $('tail-comment').textContent=`At ${h} SD: observed ${(100*hits/n).toFixed(2)}%, normal ${(100*normalTail(h)).toFixed(3)}%. The y-axis is logarithmic. Zero observed exceedances are omitted.`;
 const worst=r.map((v,i)=>({v,z:z[i],date:rows[i][0]})).sort((a,b)=>a.v-b.v).slice(0,5);
 $('extremes').replaceChildren(...worst.map(o=>{const tr=document.createElement('tr');for(const v of [o.date,o.v.toFixed(2),o.z.toFixed(2)]){const td=document.createElement('td');td.textContent=v;tr.append(td);}return tr;}));
}
for(const id of ['asset','period','threshold'])$(id).addEventListener('input',render);
render();new ResizeObserver(()=>parent.postMessage({type:'financeLabResize',height:Math.ceil(document.body.getBoundingClientRect().height)+32},location.origin)).observe(document.body);
})();
