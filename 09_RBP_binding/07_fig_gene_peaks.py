# Figure 1: how many different RBPs bind each position of DDX21, SMC1A, RPS3 (reproducible eCLIP peaks, hg19)
import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
O='./rbp/out/'
P=pd.read_csv(O+'gene_eclip_peaks_hg19.tsv',sep='\t'); P=P[P.reproducible==1]
J=pd.read_csv(O+'leafcutter_junctions_hg19.tsv',sep='\t')
tx={}
for l in open('./rbp/ref/refGene.txt'):
    f=l.split('\t')
    if f[1] in ('NM_004728','NM_006306','NM_001005') and '_' not in f[2]:
        tx[f[12]]=(int(f[4]),int(f[5]),f[3],int(f[6]),int(f[7]),[int(x) for x in f[9].strip(',').split(',')],[int(x) for x in f[10].strip(',').split(',')])
plt.rcParams.update({'font.size':9,'axes.edgecolor':'#52514e','axes.labelcolor':'#0b0b0b','xtick.color':'#52514e','ytick.color':'#52514e'})
fig,axes=plt.subplots(3,1,figsize=(11,8),facecolor='#fcfcfb')
for ax,g in zip(axes,['DDX21','SMC1A','RPS3']):
    s,e,strand,cs,ce,es,ee=tx[g]; x=np.arange(s,e); cov=np.zeros(e-s)
    for rbp,d in P[P.gene==g].groupby('rbp'):
        m=np.zeros(e-s,bool)
        for a,b in zip(d.start,d.end): m[max(a-s,0):min(b-s,e-s)]=True
        cov+=m
    ax.set_facecolor('#fcfcfb')
    ax.fill_between(x/1e3,cov,step='mid',color='#2a78d6',lw=0)
    top=cov.max(); y0=-0.12*top
    ax.hlines(y0,s/1e3,e/1e3,color='#52514e',lw=1)
    for a,b in zip(es,ee):
        h=0.06*top if (b<=cs or a>=ce) else 0.1*top
        ax.add_patch(plt.Rectangle((a/1e3,y0-h/2),(b-a)/1e3,h,color='#52514e',lw=0))
    for v in set(J[J.gene==g].intron_start_hg19)|set(J[J.gene==g].intron_end_hg19):
        ax.axvline(v/1e3,color='#eb6834',lw=0.8,ls=(0,(2,2)))
    ax.set_ylim(y0-0.12*top,top*1.1); ax.set_xlim(s/1e3,e/1e3)
    for sp in ['top','right']: ax.spines[sp].set_visible(False)
    ax.set_title(f'{g}  ({"chr"+str(P[P.gene==g].chrom.iloc[0])[3:]}, {strand} strand, hg19)   {P[P.gene==g].rbp.nunique()} RBPs with reproducible peaks',loc='left',fontsize=10,color='#0b0b0b')
    ax.set_ylabel('RBPs bound'); ax.ticklabel_format(axis='x',useOffset=False,style='plain')
axes[-1].set_xlabel('genomic position (kb)')
fig.text(0.01,0.005,'Blue = number of different RBPs with a reproducible eCLIP peak at that position (log2FC ≥ 2, -log10 p ≥ 2, same strand). Grey = exons (tall = coding). Orange dashed = splice sites of the LeafCutter junctions.',fontsize=8,color='#52514e')
plt.tight_layout(rect=(0,0.02,1,1)); plt.savefig(O+'Fig1_RBP_binding_on_splicing_hits.png',dpi=200,facecolor='#fcfcfb')
