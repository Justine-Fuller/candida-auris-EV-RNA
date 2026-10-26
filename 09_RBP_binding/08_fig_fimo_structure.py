# Figure 2: (A) motif hits in EV-enriched vs comparison RNAs; (B) structure at sites vs rest of RNA, EV vs comparison
import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
O='./rbp/out/'
F=pd.read_csv(O+'fimo_motif_enrichment.tsv',sep='\t'); S=pd.read_csv(O+'structure_per_motif.tsv',sep='\t')
plt.rcParams.update({'font.size':9,'axes.edgecolor':'#52514e','xtick.color':'#52514e','ytick.color':'#52514e'})
fig,(a,b)=plt.subplots(1,2,figsize=(11,5),facecolor='#fcfcfb')
for ax in (a,b):
    ax.set_facecolor('#fcfcfb'); [ax.spines[s].set_visible(False) for s in ['top','right']]
m=max(F.EV_hits_per_kb.max(),F.comparison_hits_per_kb.max())*1.08
a.plot([0,m],[0,m],color='#52514e',lw=1,ls=(0,(3,3)))
a.scatter(F.comparison_hits_per_kb,F.EV_hits_per_kb,s=40,color='#2a78d6',edgecolor='#fcfcfb',lw=1.5)
for r in F.nlargest(2,'EV_hits_per_kb').itertuples(): a.annotate(r.rbp,(r.comparison_hits_per_kb,r.EV_hits_per_kb),xytext=(5,-3),textcoords='offset points',color='#0b0b0b',fontsize=8)
a.set_xlim(0,m); a.set_ylim(0,m); a.set_xlabel('hits per kb, comparison RNAs'); a.set_ylabel('hits per kb, EV-enriched RNAs')
a.set_title('A. Human RBP motifs (FIMO) in C. auris RNAs\n75 motifs; none differ (all q > 0.6)',loc='left',fontsize=10)
lo=min(S.median_delta_EV.min(),S.median_delta_comparison.min())-0.03; hi=max(S.median_delta_EV.max(),S.median_delta_comparison.max())+0.03
b.plot([lo,hi],[lo,hi],color='#52514e',lw=1,ls=(0,(3,3))); b.axhline(0,color='#d0cfca',lw=0.8); b.axvline(0,color='#d0cfca',lw=0.8)
b.scatter(S.median_delta_comparison,S.median_delta_EV,s=40,color='#2a78d6',edgecolor='#fcfcfb',lw=1.5)
for r in S.loc[[S.median_delta_EV.idxmax(),S.median_delta_EV.idxmin()]].itertuples(): b.annotate(r.motif_id.split('_m')[0],(r.median_delta_comparison,r.median_delta_EV),xytext=(5,-3),textcoords='offset points',fontsize=8)
b.set_xlim(lo,hi); b.set_ylim(lo,hi)
b.set_xlabel('site minus rest, comparison RNAs'); b.set_ylabel('site minus rest, EV-enriched RNAs')
b.set_title('B. Unpaired probability at motif sites vs rest of RNA\nabove 0 = more open, below 0 = more structured',loc='left',fontsize=10)
fig.text(0.01,0.01,'Each dot is one MEME motif. Dots on the dashed diagonal behave the same in EV-enriched and comparison RNAs.',fontsize=8,color='#52514e')
plt.tight_layout(rect=(0,0.04,1,1)); plt.savefig(O+'Fig2_FIMO_and_structure.png',dpi=200,facecolor='#fcfcfb')
