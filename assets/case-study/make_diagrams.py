"""Original schematic figures for the CASE study. No cohort inference is performed."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
P=Path(__file__).resolve().parent
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'text.color':'#263748','axes.labelcolor':'#263748','svg.fonttype':'none'})
def save(fig,name):
 for ext in ['svg','png']:fig.savefig(P/f'{name}.{ext}',dpi=170,bbox_inches='tight',pad_inches=.18)
 plt.close(fig)
def box(ax,x,y,w,h,text,fill='#f0f6fb',fs=11):
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.012',edgecolor='#b9cbd9',facecolor=fill,lw=1))
 ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=fs)
def arrow(ax,a,b): ax.annotate('',xy=b,xytext=a,arrowprops=dict(arrowstyle='->',color='#74899a',lw=1.3))
# A single LD example: small, light table diagram.
fig,ax=plt.subplots(figsize=(9,2.6));ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
for x,title,vals in [(0.03,'Joint effects B',[[1,0],[0,1]]),(.60,'Marginal associations RB',[[1,.8],[.8,1]])]:
 ax.text(x+.17,.94,title,ha='center',weight='bold',fontsize=12)
 ax.text(x+.13,.75,'Type 1',ha='center');ax.text(x+.27,.75,'Type 2',ha='center')
 for i in range(2):
  ax.text(x-.012,.53-i*.28,f'SNP {i+1}',ha='right',va='center',fontsize=10)
  for j in range(2):box(ax,x+.07+j*.14,.42-i*.28,.12,.21,f'{vals[i][j]:g}',fill='#eef5fa' if vals[i][j] else '#ffffff',fs=15)
arrow(ax,(.40,.42),(.56,.42));ax.text(.48,.65,'LD\nr = 0.8',ha='center',fontsize=11)
ax.text(.18,-.04,'One active cell type per SNP',ha='center',fontsize=10)
ax.text(.78,-.04,'Each SNP associates with both types',ha='center',fontsize=10)
assert np.allclose(np.array([[1,.8],[.8,1]])@np.eye(2),[[1,.8],[.8,1]])
save(fig,'ld-simple')
# A toy exact posterior distribution with one causal SNP in A, null in B.
fig,axs=plt.subplots(1,2,figsize=(9,3.0),layout='constrained')
for ax,p,title,col in zip(axs,[[.60,.36,.04],[.02,.02,.01]],['Cell type A: credible set {SNP 1, SNP 2}','Cell type B: no credible set'],['#b9d8ee','#e0e8ef']):
 ax.barh([2,1,0],p,color=col,height=.52)
 ax.set_yticks([2,1,0],['SNP 1','SNP 2','SNP 3']);ax.set_xlim(0,1);ax.set_xlabel('Posterior inclusion probability');ax.set_title(title,fontsize=11)
 for y,v in zip([2,1,0],p):ax.text(v+.02,y,f'{v:.2f}',va='center',fontsize=10)
 for sp in ['top','right','left']:ax.spines[sp].set_visible(False)
 ax.grid(axis='x',alpha=.15);ax.set_axisbelow(True)
axs[0].text(.5,-.56,'0.60 + 0.36 = 0.96; |r₁₂| = 0.8\neGene in A',transform=axs[0].transAxes,ha='center',fontsize=10)
axs[1].text(.5,-.56,'Total PIP = 0.05\nNo eGene call in B',transform=axs[1].transAxes,ha='center',fontsize=10)
save(fig,'credible-set-example')
# Downstream analysis as three explicit transformations.
fig,ax=plt.subplots(figsize=(10,4));ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
rows=[('Gene calls / specificity','Match genes to markers\nor pathways','2 × 2 enrichment table\nOdds ratio + Fisher test'),('Fine-mapped SNP sets','Build SNP annotations\n+ LD scores + GWAS','Stratified LD-score\nregression\nHeritability enrichment'),('eQTL locus + GWAS locus','Harmonize alleles and LD\nFine-map both signals','Colocalization posterior\nShared or distinct\ncausal signal')]
for y,(a,b,c) in zip([.72,.39,.06],rows):
 box(ax,.015,y,.25,.21,a,fs=11);box(ax,.365,y,.27,.21,b,fill='#f7f8fa',fs=10.5);box(ax,.735,y,.25,.21,c,fs=10.5)
 arrow(ax,(.28,y+.105),(.35,y+.105));arrow(ax,(.65,y+.105),(.72,y+.105))
ax.text(.50,.985,'Input → statistical construction → interpretable output',ha='center',weight='bold',fontsize=13)
save(fig,'downstream-process')
print('Wrote three diagrams; LD algebra checked.')
