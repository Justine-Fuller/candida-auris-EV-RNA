# Presence of human RBP motifs (FIMO p < 1e-4, given strand) in every RNA found in EVs
import pandas as pd
A='./rbp/all/'
h=pd.read_csv(A+'fimo_ev_all.tsv',sep='\t'); h.columns=[c.lstrip('#').strip() for c in h.columns]
s=pd.read_csv(A+'sets.tsv',sep='\t'); ev=s[s.set=='EV_detected'].set_index('model_id')
enr=set(pd.read_csv('./rbp/ev/ev_sets.tsv',sep='\t').query("set=='EV_enriched'").model_id)
mk=pd.read_csv('./rbp/out/motifs_kept.tsv',sep='\t'); mk['motif_id']=mk.rbp+'_m'+mk.meme_n.astype(str)
h['rbp']=h.motif_id.str.rsplit('_m',n=1).str[0]; h['rna_class']=h.sequence_name.map(ev.rna_class); h['EV_enriched_6F']=h.sequence_name.isin(enr)
sites=h[['sequence_name','rna_class','EV_enriched_6F','rbp','motif_id','start','stop','strand','matched_sequence','score','p-value']].rename(columns={'sequence_name':'model_id','stop':'end'})
sites=sites.sort_values(['model_id','start']); sites.to_csv(A+'motif_sites_all_EV_RNAs.tsv.gz',sep='\t',index=False)
kb=ev.len_whole.sum()/1e3; n=len(ev)
per_motif=[]
for m in mk.itertuples():
    x=h[h.motif_id==m.motif_id]; st=x[x['p-value']<1e-5]
    per_motif.append(dict(motif_id=m.motif_id,rbp=m.rbp,consensus=m.consensus,meme_E=m.E,
        EV_RNAs_with_site=x.sequence_name.nunique(),pct_of_EV_RNAs=round(100*x.sequence_name.nunique()/n,1),
        mRNAs_with_site=x[x.rna_class=='PCG_CDS'].sequence_name.nunique(),lncRNAs_with_site=x[x.rna_class=='lncRNA'].sequence_name.nunique(),
        total_sites=len(x),sites_per_kb=round(len(x)/kb,3),EV_RNAs_with_strong_site_p1e5=st.sequence_name.nunique()))
P=pd.DataFrame(per_motif).sort_values('EV_RNAs_with_site',ascending=False); P.to_csv(A+'motif_presence_per_motif.tsv',sep='\t',index=False)
per_rna=h.groupby('sequence_name').agg(n_sites=('motif_id','size'),n_motifs=('motif_id','nunique'),n_RBPs=('rbp','nunique'),RBPs=('rbp',lambda r:','.join(sorted(set(r))))).reindex(ev.index)
per_rna['n_sites']=per_rna.n_sites.fillna(0).astype(int); per_rna['n_motifs']=per_rna.n_motifs.fillna(0).astype(int); per_rna['n_RBPs']=per_rna.n_RBPs.fillna(0).astype(int)
per_rna.insert(0,'rna_class',ev.rna_class); per_rna.insert(1,'length',ev.len_whole); per_rna.insert(2,'EV_enriched_6F',per_rna.index.isin(enr))
per_rna.index.name='model_id'; per_rna.to_csv(A+'motif_presence_per_RNA.tsv',sep='\t')
print('EV RNAs',n,'with >=1 site',(per_rna.n_sites>0).sum(),'median RBPs per RNA',per_rna.n_RBPs.median(),'median sites',per_rna.n_sites.median())
print(P.head(12).to_string()); print(P.tail(5).to_string())
