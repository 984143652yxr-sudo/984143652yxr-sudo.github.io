(function(){
const frame=document.getElementById('finance-lab');
if(!frame)return;
window.addEventListener('message',e=>{
if(e.origin!==location.origin||e.source!==frame.contentWindow||e.data?.type!=='financeLabResize')return;
const h=Number(e.data.height);
if(Number.isFinite(h)&&h>300&&h<5000)frame.style.height=h+'px';
});
})();
