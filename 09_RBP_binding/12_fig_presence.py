# Figure 3: share of the 7,152 RNAs found in EVs that carry at least one site for each human RBP motif
import pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
A='./rbp/all/'
P=pd.read_csv(A+'motif_presence_per_motif.tsv',sep='\t').sort_values('pct_of_EV_RNAs')
P=P[P.EV_RNAs_with_site>0]
fig,ax=plt.subplots(figsize=(8,13),facecolor='#fcfcfb'); ax.set_facecolor('#fcfcfb')
y=range(len(P)); ax.barh(list(y),P.pct_of_EV_RNAs,height=0.75,color='#2a78d6')
ax.set_yticks(list(y)); ax.set_yticklabels([f'{r} {c}' for r,c in zip(P.rbp,P.consensus)],fontsize=7,family='monospace')
for i,(v,n) in enumerate(zip(P.pct_of_EV_RNAs,P.EV_RNAs_with_site)): ax.text(v+0.4,i,f'{n:,}',va='center',fontsize=6.5,color='#52514e')
for s in ['top','right','left']: ax.spines[s].set_visible(False)
ax.tick_params(axis='y',length=0); ax.xaxis.grid(True,color='#e6e5e0',lw=0.6); ax.set_axisbelow(True)
ax.set_xlabel('% of RNAs found in EVs with ≥ 1 motif site (FIMO p < 1e-4)')
ax.set_title('Human RBP motifs present in C. auris EV RNAs\n7,152 RNAs (5,146 mRNA, 2,006 lncRNA); numbers = RNAs with a site',loc='left',fontsize=10)
plt.tight_layout(); plt.savefig(A+'Fig3_motif_presence_all_EV_RNAs.png',dpi=200,facecolor='#fcfcfb')
