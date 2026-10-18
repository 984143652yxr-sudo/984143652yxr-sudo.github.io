(function(){
window.addEventListener('message',e=>{
if(e.origin!==location.origin||e.data?.type!=='financeLabResize')return;
const frame=[...document.querySelectorAll('iframe[data-finance-lab]')].find(f=>e.source===f.contentWindow);
const h=Number(e.data.height);
if(frame&&Number.isFinite(h)&&h>300&&h<6000)frame.style.height=h+'px';
});
})();
