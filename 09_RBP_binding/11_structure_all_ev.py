# Structure at motif sites vs the rest of the same RNA, for every RNA found in EVs (same folding as 06_structure.py)
import RNA, pandas as pd, numpy as np, multiprocessing as mp
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
A='./rbp/all/'
def fa(f):
    d={};k=None
    for l in open(f):
        l=l.strip()
        if l.startswith('>'): k=l[1:].split()[0]; d[k]=''
        else: d[k]+=l.upper()
    return d
def unpaired(item):
    k,s=item; up=RNA.pfl_fold_up(s.replace('T','U'),1,150,100)
    return k,np.array([up[i][1] for i in range(1,len(s)+1)],dtype=np.float32)
if __name__=='__main__':
    seqs=fa(A+'ev_all.fa')
    with mp.Pool(4) as p: UP=dict(p.imap_unordered(unpaired,sorted(seqs.items(),key=lambda x:-len(x[1])),chunksize=4))
    np.savez_compressed(A+'unpaired_prob_ev_all.npz',**UP)
    h=pd.read_csv(A+'fimo_ev_all.tsv',sep='\t'); h.columns=[c.lstrip('#').strip() for c in h.columns]
    rows=[]
    for (mid,seq),x in h.groupby(['motif_id','sequence_name']):
        u=UP[seq]; m=np.zeros(len(u),bool)
        for s,e in zip(x.start,x.stop): m[s-1:e]=True
        rows.append((mid,seq,len(x),float(u[m].mean()),float(u[~m].mean())))
    S=pd.DataFrame(rows,columns=['motif_id','model_id','n_sites','site_unpaired','rest_unpaired']); S['delta']=S.site_unpaired-S.rest_unpaired
    S.to_csv(A+'structure_per_site_ev_all.tsv.gz',sep='\t',index=False)
    per=S.groupby('model_id').apply(lambda d:pd.Series({'site':np.average(d.site_unpaired,weights=d.n_sites),'rest':d.rest_unpaired.mean()}))
    w=wilcoxon(per.site,per.rest)
    summ=dict(n_RNAs=len(per),median_site_unpaired=per.site.median(),median_rest_unpaired=per.rest.median(),median_delta=(per.site-per.rest).median(),
        pct_RNAs_site_more_structured=100*(per.site<per.rest).mean(),p_wilcoxon=w.pvalue)
    pd.Series(summ).to_csv(A+'structure_summary_ev_all.tsv',sep='\t',header=False)
    pr=[]
    for mid,x in S.groupby('motif_id'):
        if len(x)<10: continue
        pr.append(dict(motif_id=mid,rbp=mid.rsplit('_m',1)[0],n_EV_RNAs=len(x),median_site_unpaired=x.site_unpaired.median(),median_rest_unpaired=x.rest_unpaired.median(),
            median_delta=x.delta.median(),pct_RNAs_site_more_structured=round(100*(x.delta<0).mean(),1),p=wilcoxon(x.delta).pvalue))
    P=pd.DataFrame(pr); P['q']=multipletests(P.p,method='fdr_bh')[1]
    P['verdict']=np.where(P.q>=0.05,'no difference',np.where(P.median_delta<0,'site MORE structured','site LESS structured'))
    P.sort_values('median_delta').to_csv(A+'structure_per_motif_ev_all.tsv',sep='\t',index=False)
    print(summ); print(P.verdict.value_counts()); print(P.sort_values('median_delta').head(8).to_string())
