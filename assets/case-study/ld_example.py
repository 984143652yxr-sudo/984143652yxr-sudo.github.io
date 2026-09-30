"""Reproduce a deterministic LD illustration. This does not run CASE or use cohort data.
Dependencies: numpy, matplotlib. Run from any directory; outputs sit beside this file.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path(__file__).resolve().parent
R=np.array([[1.,.8],[.8,1.]])
cases=[('Different causal SNPs',np.array([[1.,0.],[0.,1.]])),
       ('Opposing effects in cell type 2',np.array([[1.,1.],[0.,-1.]]))]
fig,axs=plt.subplots(2,2,figsize=(9,6.6),layout='constrained')
for row,(label,B) in enumerate(cases):
 marginal=R@B
 assert np.allclose(np.linalg.solve(R,marginal),B)
 for col,(title,A) in enumerate([('Joint effects B',B),('Expected marginal effects RB',marginal)]):
  ax=axs[row,col];ax.imshow(A,vmin=-1,vmax=1,cmap='RdBu_r',aspect='auto')
  ax.set_xticks([0,1],['Cell type 1','Cell type 2']);ax.set_yticks([0,1],['SNP 1','SNP 2'])
  ax.set_title(title,fontsize=12,pad=10)
  for i in range(2):
   for j in range(2):ax.text(j,i,f'{A[i,j]:.1f}',ha='center',va='center',fontsize=21,color='white' if abs(A[i,j])>.65 else '#172b4d')
  for spine in ax.spines.values():spine.set_visible(False)
 axs[row,0].set_ylabel(label,fontsize=11,labelpad=12)
fig.suptitle('LD mixes joint effects into marginal associations\nTwo SNPs with correlation r = 0.8 · deterministic teaching example',fontsize=14)
fig.savefig(out/'ld-effects.png',dpi=190)
fig.savefig(out/'ld-effects.svg')
print('R condition number:',np.linalg.cond(R))
for label,B in cases: print(label,'\nB =',B,'\nRB =',R@B)
