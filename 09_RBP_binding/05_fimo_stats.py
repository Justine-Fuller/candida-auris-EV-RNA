# Step 4b: are human-RBP motifs more common in EV-enriched RNAs than in comparison RNAs?
# per motif: logistic regression  has_hit ~ EV_enriched + log(length) + GC ; BH across motifs
import pandas as pd, numpy as np, statsmodels.formula.api as smf, warnings
from statsmodels.stats.multitest import multipletests
warnings.filterwarnings('ignore')
O='./rbp/out/'; EV='./rbp/ev/'
def fa(f):
    d={};k=None
    for l in open(f):
        l=l.strip()
        if l.startswith('>'): k=l[1:].split()[0]; d[k]=''
        else: d[k]+=l.upper()
    return d
seqs={**fa(EV+'ev_enriched.fa'),**fa(EV+'comparison.fa')}
sets=pd.read_csv(EV+'ev_sets.tsv',sep='\t').set_index('model_id')
genes=pd.DataFrame({'len':{k:len(v) for k,v in seqs.items()},'gc':{k:(v.count('G')+v.count('C'))/len(v) for k,v in seqs.items()}})
genes['ev']=(sets.loc[genes.index,'set']=='EV_enriched').astype(int); genes['loglen']=np.log(genes['len'])
hits=pd.concat([pd.read_csv(O+f'fimo_{s}.tsv',sep='\t') for s in ['ev_enriched','comparison']])
hits.columns=[c.lstrip('#').strip() for c in hits.columns]
shuf=pd.read_csv(O+'fimo_ev_enriched_shuffled.tsv',sep='\t'); shuf.columns=[c.lstrip('#').strip() for c in shuf.columns]
motifs=pd.read_csv(O+'motifs_kept.tsv',sep='\t'); motifs['motif_id']=motifs.rbp+'_m'+motifs.meme_n.astype(str)
rows=[]
nEV=genes.ev.sum(); nC=len(genes)-nEV; kbEV=genes.len[genes.ev==1].sum()/1e3; kbC=genes.len[genes.ev==0].sum()/1e3
for m in motifs.itertuples():
    h=hits[hits.motif_id==m.motif_id]; g=genes.copy(); g['hit']=g.index.isin(h.sequence_name).astype(int)
    hev=h[h.sequence_name.isin(g.index[g.ev==1])]; hc=h[h.sequence_name.isin(g.index[g.ev==0])]
    try:
        fit=smf.logit('hit ~ ev + loglen + gc',data=g).fit(disp=0); OR=np.exp(fit.params['ev']); p=fit.pvalues['ev']
    except Exception: OR,p=np.nan,np.nan
    rows.append(dict(motif_id=m.motif_id,rbp=m.rbp,consensus=m.consensus,meme_E=m.E,
        EV_RNAs_with_hit=g.hit[g.ev==1].sum(),EV_frac=round(g.hit[g.ev==1].mean(),3),
        comparison_frac=round(g.hit[g.ev==0].mean(),3),
        EV_hits_per_kb=round(len(hev)/kbEV,3),comparison_hits_per_kb=round(len(hc)/kbC,3),
        shuffled_EV_hits_per_kb=round((shuf.motif_id==m.motif_id).sum()/kbEV,3),
        adj_odds_ratio=round(OR,3),p=p))
R=pd.DataFrame(rows); ok=R.p.notna()
R.loc[ok,'q']=multipletests(R.p[ok],method='fdr_bh')[1]
R['obs_over_shuffled']=(R.EV_hits_per_kb/R.shuffled_EV_hits_per_kb.replace(0,np.nan)).round(2)
R=R.sort_values('p'); R.to_csv(O+'fimo_motif_enrichment.tsv',sep='\t',index=False)
print('RNAs: EV',nEV,'comparison',nC); print('motifs tested',len(R),'q<0.05',(R.q<0.05).sum(),' up',((R.q<0.05)&(R.adj_odds_ratio>1)).sum(),' down',((R.q<0.05)&(R.adj_odds_ratio<1)).sum())
print(R.head(15).to_string())
