import pandas as pd, gseapy as gp, json, sys
sys.path.insert(0,"."); from ora import load_maps
go,kg=load_maps("emapper.tsv"); gn=json.load(open("go_names.json"))
t=pd.read_csv("stranded_pcg_calls.csv").dropna(subset=["stat"])
rnk=t[["Feature","stat"]].sort_values("stat",ascending=False)
uni=set(rnk.Feature)
for name,m in [("GO",go),("KEGG",kg)]:
    m=m[m.gene.isin(uni)]
    gs={k:list(v) for k,v in m.groupby("term").gene}
    gs={k:v for k,v in gs.items() if 3<=len(v)<=500}
    r=gp.prerank(rnk=rnk,gene_sets=gs,permutation_num=1000,min_size=3,max_size=500,seed=1,threads=8,outdir=None,verbose=False).res2d
    r["FDR q-val"]=r["FDR q-val"].astype(float); r["NES"]=r["NES"].astype(float)
    s=r[r["FDR q-val"]<0.05].copy()
    s["name"]=s.Term.map(lambda x: gn.get(x) or x)
    s.to_csv(f"GSEA_{name}_stranded.csv",index=False)
    print(name,"sig",len(s),"pos",(s.NES>0).sum(),"neg",(s.NES<0).sum())
    print(" top+:", "; ".join(s[s.NES>0].sort_values("NES",ascending=False).name.head(12)))
    print(" top-:", "; ".join(s[s.NES<0].sort_values("NES").name.head(12)))
