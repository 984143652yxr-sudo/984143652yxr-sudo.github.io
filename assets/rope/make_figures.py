from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

ROOT = Path(__file__).parent
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11, 'svg.fonttype': 'none', 'svg.hashsalt': 'rope-study'})
blue, teal, gray = '#3779ad', '#259b91', '#526170'

def save(fig, name):
    fig.savefig(ROOT / (name+'.svg'), bbox_inches='tight', metadata={'Date': None})
    fig.savefig(ROOT / (name+'.png'), bbox_inches='tight', dpi=170)
    p=ROOT/(name+'.svg')
    p.write_text('\n'.join(line.rstrip() for line in p.read_text().splitlines())+'\n')
    plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(10, 4.3))
for ax, m, n in zip(axes, [1, 4], [3, 6]):
    ax.add_patch(plt.Circle((0,0), 1, fill=False, color='#d7e0e6', lw=1.5))
    for p, col, lab in [(m, blue, 'query'), (n, teal, 'key')]:
        ang=p*np.pi/6
        xy=(np.cos(ang),np.sin(ang))
        ax.add_patch(FancyArrowPatch((0,0),xy,arrowstyle='-|>',mutation_scale=18,color=col,lw=2.5))
        ax.text(xy[0]*1.2,xy[1]*1.2,lab,color=col,ha='center',va='center')
    arc=np.linspace(m*np.pi/6,n*np.pi/6,80)
    ax.plot(.42*np.cos(arc),.42*np.sin(arc),color=gray)
    ax.text(0,-1.25,'Gap = 2 tokens     Angle = 60°\nDot product = cos(60°) = 0.5',ha='center',va='top',linespacing=1.8)
    ax.set_title(f'Positions m = {m}, n = {n}',pad=17)
    ax.set(xlim=(-1.6,1.6),ylim=(-1.8,1.5),aspect='equal')
    ax.axis('off')
fig.suptitle('A shared position shift preserves the relative angle',fontsize=15,y=1.04)
fig.text(.5,-.03,'Fixed unrotated query = key = (1, 0); rotation frequency θ = π/6 per token.',ha='center',color=gray)
fig.tight_layout()
save(fig,'relative-rotation')

fig, axes=plt.subplots(2,1,figsize=(9,5.8),layout='constrained')
t=np.arange(101)
for f,col in [(1.,blue),(.1,teal),(.01,'#b58342')]:
    axes[0].plot(t,np.cos(t*f),label=f'θ = {f:g}',color=col,lw=1.7)
axes[0].set(title='Different coordinate pairs rotate at different speeds',ylabel='cos(position × θ)',xlabel='Token position')
axes[0].legend(loc='upper right',ncol=3,framealpha=.95)
r=np.linspace(0,24,300)
axes[1].plot(r,np.cos(r*np.pi/6),color=blue,lw=2)
axes[1].scatter([0,6,12,18,24],[1,-1,1,-1,1],color=blue,s=25,zorder=3)
axes[1].set(title='One pair: the score can rise again at a larger distance',xlabel='Relative position n − m',ylabel='Rotated dot product',ylim=(-1.25,1.25))
for ax in axes:
    ax.spines[['top','right']].set_visible(False)
    ax.grid(alpha=.15)
fig.get_layout_engine().set(h_pad=.18)
save(fig,'frequencies')
