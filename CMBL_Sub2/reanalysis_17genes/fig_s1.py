import json,sys,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
sys.stdout.reconfigure(encoding='utf8')
plt.rcParams['font.family']='Arial' if any('Arial'==f.name for f in font_manager.fontManager.ttflist) else 'DejaVu Sans'
plt.rcParams['axes.unicode_minus']=True
if plt.rcParams['font.family']==['Arial']:
    plt.rcParams.update({'mathtext.fontset':'custom','mathtext.rm':'Arial','mathtext.it':'Arial:italic','mathtext.bf':'Arial:bold'})
NPG=dict(dark='#3C5488',light='#4DBBD5',salmon='#F39B7F',red='#E64B35',grey='#C9CEE0')
A=np.load('perm_arrays_17.npz'); TC=json.load(open('table_counts.json'))['17']
obs={'Raw':66,'Hes1-normalized':120,'Kcnj11-normalized':84,'Ins1-normalized':80}
labels={'Raw':(None,'Raw','expression'),'Hes1-normalized':('Hes1','-','normalized'),'Kcnj11-normalized':('Kcnj11','-','normalized'),'Ins1-normalized':('Ins1','-','normalized')}
from matplotlib.offsetbox import TextArea, HPacker, VPacker, AnnotationBbox
def cat_axis(ax,entries,fs=7.5):
    """two-line category labels; a gene symbol is set in italic (digits included, which mathtext would keep upright)"""
    ax.set_xticklabels([])
    for x,(gene,l1,l2) in enumerate(entries):
        tp=dict(fontsize=fs)
        first=HPacker(children=[TextArea(gene,textprops=dict(fontsize=fs,fontstyle='italic')),TextArea(l1,textprops=tp)],align='baseline',pad=0,sep=0) if gene else TextArea(l1,textprops=tp)
        box=VPacker(children=[first,TextArea(l2,textprops=tp)],align='center',pad=0,sep=2)
        ax.add_artist(AnnotationBbox(box,(x,0),xycoords=('data','axes fraction'),xybox=(0,-7),boxcoords='offset points',box_alignment=(0.5,1.0),frameon=False,pad=0))
fig=plt.figure(figsize=(6.8,5.9),dpi=300)
gs=fig.add_gridspec(2,2,width_ratios=[2.2,1],height_ratios=[1,1],hspace=0.55,wspace=0.28,left=0.1,right=0.985,top=0.94,bottom=0.1)
def violin(ax,x,arr,observed,label_dx=0.18):
    arr=np.asarray(arr)
    if arr.std()>0:
        parts=ax.violinplot([arr],positions=[x],widths=0.8,showextrema=False)
        for b in parts['bodies']: b.set_facecolor(NPG['grey']); b.set_alpha(0.85); b.set_edgecolor('none')
        q1,q3=np.percentile(arr,[25,75]); lo,hi=np.percentile(arr,[2.5,97.5])
        ax.plot([x,x],[lo,hi],color=NPG['dark'],lw=1.0,zorder=3); ax.plot([x,x],[q1,q3],color=NPG['dark'],lw=3.2,zorder=3)
    else:
        ax.plot([x-0.25,x+0.25],[arr[0],arr[0]],color=NPG['dark'],lw=1.2,zorder=3)
    ax.scatter([x],[observed],marker='D',s=42,facecolor='white',edgecolor='black',linewidth=1.3,zorder=5)
    ax.text(x+label_dx,observed,str(observed),va='center',ha='left',fontsize=7.5,zorder=6)
# (a) left
ax=fig.add_subplot(gs[0,0])
for i,k in enumerate(obs): violin(ax,i,A[k],obs[k],label_dx=(0.32 if A[k].std()==0 else 0.18))
ax.set_xticks(range(4)); cat_axis(ax,[labels[k] for k in obs])
ax.set_ylabel('Pairs with pooled $\\mathit{p}$ < 0.05',fontsize=8,fontweight='bold'); ax.set_ylim(55,128)
ax.set_title('Observed counts vs group-preserving null',fontsize=8.5,fontweight='bold')
ax.grid(axis='y',color='#E6E6E6',lw=0.8); ax.set_axisbelow(True)
for s in ('top','right','left'): ax.spines[s].set_visible(False)
ax.tick_params(axis='y',labelsize=7.5,length=0)
fig.text(0.005,0.965,'(a)',fontsize=10,fontweight='bold')
# (a) right
ax=fig.add_subplot(gs[0,1])
violin(ax,0,A['cons'],21); violin(ax,1,A['cons_sign'],21)
ax.set_xticks([0,1]); ax.set_xticklabels(['All four\nschemes','+ same\nsign'],fontsize=7.5)
ax.set_ylabel('Consensus pairs',fontsize=8,fontweight='bold'); ax.set_ylim(12,34)
ax.set_title('Consensus counts',fontsize=8.5,fontweight='bold')
ax.grid(axis='y',color='#E6E6E6',lw=0.8); ax.set_axisbelow(True)
for s in ('top','right','left'): ax.spines[s].set_visible(False)
ax.tick_params(axis='y',labelsize=7.5,length=0)
# (b)
ax=fig.add_subplot(gs[1,:])
groups=[('Raw',labels['Raw']),('Hes1-normalized',labels['Hes1-normalized']),('Kcnj11-normalized',labels['Kcnj11-normalized']),('Ins1-normalized',labels['Ins1-normalized']),('all',(None,'All four','schemes'))]
series=[('pooled_p','pooled $\\mathit{p}$ < 0.05',NPG['dark']),('pooled_q','pooled BH $\\mathit{q}$ < 0.05',NPG['light']),('adj_p','adjusted $\\mathit{p}$ < 0.05',NPG['salmon']),('adj_q','adjusted BH $\\mathit{q}$ < 0.05',NPG['red'])]
w=0.19
for si,(key,lab,col) in enumerate(series):
    xs=[gi+(si-1.5)*w for gi in range(len(groups))]; vals=[TC[g][key] for g,_ in groups]
    ax.bar(xs,vals,width=w*0.94,color=col,label=lab,zorder=3)
    for x,v in zip(xs,vals): ax.text(x,v+2,str(v),rotation=90,ha='center',va='bottom',fontsize=6.8)
ax.set_xticks(range(len(groups))); cat_axis(ax,[l for _,l in groups])
ax.set_ylim(0,175); ax.set_ylabel('Significant pairs',fontsize=8,fontweight='bold')
ax.set_title('Significant pairs before and after adjustment for treatment ($\\mathit{df}$ = 8); only |$\\mathit{r}$| ≥ 0.90 can pass BH at $\\mathit{q}$ < 0.05',fontsize=8,fontweight='bold',pad=18)
ax.legend(ncol=4,frameon=False,fontsize=7,loc='lower center',bbox_to_anchor=(0.5,1.0),handlelength=1.2,columnspacing=1.4)
ax.grid(axis='y',color='#E6E6E6',lw=0.8); ax.set_axisbelow(True)
for s in ('top','right','left'): ax.spines[s].set_visible(False)
ax.tick_params(axis='y',labelsize=7.5,length=0)
fig.text(0.005,0.485,'(b)',fontsize=10,fontweight='bold')
fig.savefig('fig_s1_17genes.png',dpi=300,facecolor='white'); print('saved fig_s1_17genes.png')
