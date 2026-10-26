# Figure 4: per motif, unpaired probability at the motif site minus the rest of the same RNA, across all RNAs found in EVs
import pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
A='./rbp/all/'
P=pd.read_csv(A+'structure_per_motif_ev_all.tsv',sep='\t').sort_values('median_delta',ascending=False)
mk=pd.read_csv('./rbp/out/motifs_kept.tsv',sep='\t'); mk['motif_id']=mk.rbp+'_m'+mk.meme_n.astype(str); cons=dict(zip(mk.motif_id,mk.consensus))
col=['#d0cfca' if v=='no difference' else ('#eb6834' if d<0 else '#2a78d6') for v,d in zip(P.verdict,P.median_delta)]
fig,ax=plt.subplots(figsize=(8,13),facecolor='#fcfcfb'); ax.set_facecolor('#fcfcfb')
y=range(len(P)); ax.barh(list(y),P.median_delta,height=0.75,color=col)
ax.set_yticks(list(y)); ax.set_yticklabels([f'{r} {cons[m]}' for r,m in zip(P.rbp,P.motif_id)],fontsize=7,family='monospace')
ax.axvline(0,color='#52514e',lw=0.8); ax.tick_params(axis='y',length=0)
for s in ['top','right','left']: ax.spines[s].set_visible(False)
ax.xaxis.grid(True,color='#e6e5e0',lw=0.6); ax.set_axisbelow(True)
ax.set_xlabel('unpaired probability: motif site minus rest of RNA (median across EV RNAs)')
ax.set_title('Are motif sites in EV RNAs more folded than the rest of the RNA?\nOrange = site MORE structured (more paired)   Blue = site LESS structured (more open)   Grey = no difference',loc='left',fontsize=9)
plt.tight_layout(); plt.savefig(A+'Fig4_structure_at_sites_all_EV_RNAs.png',dpi=200,facecolor='#fcfcfb')
