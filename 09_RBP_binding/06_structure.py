# Step 5: are RBP motif sites in EV-enriched RNAs more structured (more base-paired) than the rest of the RNA?
# RNAplfold-style local folding (ViennaRNA): W=150, L=100, unpaired probability per nucleotide
import RNA, pandas as pd, numpy as np, multiprocessing as mp
from scipy.stats import wilcoxon, mannwhitneyu
O='./rbp/out/'; EV='./rbp/ev/'
def fa(f):
    d={};k=None
    for l in open(f):
        l=l.strip()
        if l.startswith('>'): k=l[1:].split()[0]; d[k]=''
        else: d[k]+=l.upper()
    return d
def unpaired(item):
    k,s=item
    up=RNA.pfl_fold_up(s.replace('T','U'),1,150,100)   # u=1 -> per-nt unpaired prob
    return k,np.array([up[i][1] for i in range(1,len(s)+1)])
if __name__=='__main__':
    seqs={**fa(EV+'ev_enriched.fa'),**fa(EV+'comparison.fa')}
    sets=pd.read_csv(EV+'ev_sets.tsv',sep='\t').set_index('model_id')['set']
    with mp.Pool(4) as p: UP=dict(p.map(unpaired,seqs.items(),chunksize=8))
    np.savez_compressed(O+'unpaired_prob.npz',**UP)
    hits=pd.concat([pd.read_csv(O+f'fimo_{s}.tsv',sep='\t') for s in ['ev_enriched','comparison']])
    hits.columns=[c.lstrip('#').strip() for c in hits.columns]
    rows=[]
    for (mid,seq),h in hits.groupby(['motif_id','sequence_name']):
        u=UP[seq]; mask=np.zeros(len(u),bool)
        for s,e in zip(h.start,h.stop): mask[s-1:e]=True
        rows.append(dict(motif_id=mid,rbp=mid.rsplit('_m',1)[0],model_id=seq,set=sets[seq],n_sites=len(h),
            site_unpaired=u[mask].mean(),rest_unpaired=u[~mask].mean()))
    S=pd.DataFrame(rows); S['delta']=S.site_unpaired-S.rest_unpaired   # negative = site MORE structured than rest
    S.to_csv(O+'structure_per_site.tsv',sep='\t',index=False)
    # per RNA (all motifs pooled): one value per RNA, so each molecule counts once
    out=[]
    for st in ['EV_enriched','comparison']:
        x=S[S.set==st]
        per_rna=x.groupby('model_id').apply(lambda d:pd.Series({'site':np.average(d.site_unpaired,weights=d.n_sites),'rest':d.rest_unpaired.mean()}))
        w=wilcoxon(per_rna.site,per_rna.rest)
        out.append(dict(level='all motifs pooled',set=st,n_RNAs=len(per_rna),median_site_unpaired=per_rna.site.median().round(3),
            median_rest_unpaired=per_rna.rest.median().round(3),median_delta=(per_rna.site-per_rna.rest).median().round(3),p_wilcoxon=w.pvalue))
    ev=S[S.set=='EV_enriched']; cp=S[S.set=='comparison']
    out.append(dict(level='EV vs comparison delta (per site-RNA pair)',set='EV_enriched vs comparison',n_RNAs=f'{ev.model_id.nunique()} / {cp.model_id.nunique()}',
        median_delta=f'{ev.delta.median():.3f} vs {cp.delta.median():.3f}',p_wilcoxon=mannwhitneyu(ev.delta,cp.delta).pvalue))
    pd.DataFrame(out).to_csv(O+'structure_summary.tsv',sep='\t',index=False)
    # per RBP motif in EV-enriched RNAs
    from statsmodels.stats.multitest import multipletests
    pr=[]
    for mid,x in ev.groupby('motif_id'):
        if len(x)<10: continue
        c=cp[cp.motif_id==mid]
        pr.append(dict(motif_id=mid,n_EV_RNAs=len(x),median_delta_EV=x.delta.median(),p_site_vs_rest=wilcoxon(x.delta).pvalue,
            median_delta_comparison=c.delta.median() if len(c) else np.nan,
            p_EV_vs_comparison=mannwhitneyu(x.delta,c.delta).pvalue if len(c)>=10 else np.nan))
    P=pd.DataFrame(pr); P['q_site_vs_rest']=multipletests(P.p_site_vs_rest,method='fdr_bh')[1]
    ok=P.p_EV_vs_comparison.notna(); P.loc[ok,'q_EV_vs_comparison']=multipletests(P.p_EV_vs_comparison[ok],method='fdr_bh')[1]
    P.sort_values('p_site_vs_rest').to_csv(O+'structure_per_motif.tsv',sep='\t',index=False)
    print(pd.DataFrame(out).to_string()); print(P.sort_values('p_site_vs_rest').head(15).to_string())
    print('motifs: site more structured (q<0.05, delta<0):',((P.q_site_vs_rest<0.05)&(P.median_delta_EV<0)).sum(),' less structured:',((P.q_site_vs_rest<0.05)&(P.median_delta_EV>0)).sum(),'of',len(P))
