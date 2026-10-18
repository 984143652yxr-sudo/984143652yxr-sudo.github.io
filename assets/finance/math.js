/* Historical diagnostics, using percent returns throughout. No analytics or API calls. */
(function(){
'use strict';
const mean=a=>a.reduce((x,y)=>x+y,0)/a.length;
const sd=a=>{const m=mean(a);return Math.sqrt(a.reduce((s,x)=>s+(x-m)**2,0)/(a.length-1));};
function quantile(a,p){const s=[...a].sort((x,y)=>x-y),h=(s.length-1)*p,i=Math.floor(h);return s[i]+(s[Math.min(i+1,s.length-1)]-s[i])*(h-i);}
function es(a,alpha){const l=a.map(x=>-x).sort((x,y)=>y-x),mass=l.length*(1-alpha),k=Math.floor(mass);return (l.slice(0,k).reduce((s,x)=>s+x,0)+(mass-k)*(l[k]||0))/mass;}
function tailContributions(a,b,w,alpha){
const ordered=a.map((v,i)=>({i,loss:-(w*v+(1-w)*b[i])})).sort((x,y)=>y.loss-x.loss);
const mass=a.length*(1-alpha);let remain=mass,ca=0,cb=0;const indices=[];
for(const o of ordered){if(remain<=1e-12)break;const f=Math.min(1,remain);ca+=f*(-w*a[o.i]);cb+=f*(-(1-w)*b[o.i]);indices.push(o.i);remain-=f;}
return {a:ca/mass,b:cb/mass,indices};
}
function bestAllocation(risks,means,budget){let best=-1;for(let i=0;i<risks.length;i++)if(risks[i]<=budget+1e-12&&(best<0||means[i]>means[best]+1e-12))best=i;return best;}
// Acklam inverse-normal approximation, adequate for plotting positions in this sample.
function invnorm(p){const a=[-39.6968302866538,220.946098424521,-275.928510446969,138.357751867269,-30.6647980661472,2.50662827745924],b=[-54.4760987982241,161.585836858041,-155.698979859887,66.8013118877197,-13.2806815528857],c=[-.00778489400243029,-.322396458041136,-2.40075827716184,-2.54973253934373,4.37466414146497,2.93816398269878],d=[.00778469570904146,.32246712907004,2.445134137143,3.75440866190742];let q,r;
if(p<.02425){q=Math.sqrt(-2*Math.log(p));return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5])/((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);}
if(p>1-.02425)return -invnorm(1-p);q=p-.5;r=q*q;return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q/(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);}
function normalTail(x){const z=x/Math.sqrt(2),t=1/(1+.3275911*z);return (((((1.061405429*t-1.453152027)*t)+1.421413741)*t-.284496736)*t+.254829592)*t*Math.exp(-z*z);}
const gaussianES=(a,alpha)=>-mean(a)+sd(a)*Math.exp(-(invnorm(alpha)**2)/2)/Math.sqrt(2*Math.PI)/(1-alpha);
// Expose pure functions for reproducible cross-language checks.
if(typeof module!=='undefined')module.exports={mean,sd,quantile,es,invnorm,normalTail,gaussianES,tailContributions,bestAllocation};
if(typeof window!=='undefined')window.FinanceMath={mean,sd,quantile,es,invnorm,normalTail,gaussianES,tailContributions,bestAllocation};
})();
