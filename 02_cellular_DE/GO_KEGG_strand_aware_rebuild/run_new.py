import pandas as pd, numpy as np, sys, json
sys.path.insert(0,".")
from ora import load_maps, enricher
U="./"
go,kg=load_maps("emapper.tsv")
R=U+"Differential_expression_analysis/03_Cellular_DE_analysis/mRNA_DE_analysis/STRANDED_recount_2026-10-14/results/"
b=pd.read_csv(R+"new_stranded_B8441_results.csv"); m=pd.read_csv(R+"new_stranded_MMC1_results.csv")
t=b[["Feature","log2FoldChange","padj","stat"]].merge(m[["Feature","log2FoldChange","padj"]],on="Feature",suffixes=("_B8441","_MMC1"))
t=t[t.Feature.str.startswith("B9J08_")]
t=t.rename(columns={"log2FoldChange_B8441":"LFC_B8441","log2FoldChange_MMC1":"LFC_MMC1"})
def call(l,p): return np.where((p<0.05)&(l>0),"up",np.where((p<0.05)&(l<0),"down","ns"))
t["B"]=call(t.LFC_B8441,t.padj_B8441.fillna(1)); t["M"]=call(t.LFC_MMC1,t.padj_MMC1.fillna(1))
t.to_csv("stranded_pcg_calls.csv",index=False)
c=pd.crosstab(t.B,t.M); print(c); print("n",len(t))
uGO=sorted(set(t.Feature)&set(go.gene)); uKG=sorted(set(t.Feature)&set(kg.gene)); print("universe",len(uGO),len(uKG))
sets={"B8441_up":t[t.B=="up"],"B8441_down":t[t.B=="down"],"MMC1_up":t[t.M=="up"],"MMC1_down":t[t.M=="down"],
"Shared_up":t[(t.B=="up")&(t.M=="up")],"Shared_down":t[(t.B=="down")&(t.M=="down")],
"Opp_upB_downM":t[(t.B=="up")&(t.M=="down")],"Opp_downB_upM":t[(t.B=="down")&(t.M=="up")]}
gn=json.load(open("go_names.json"))
kn={}
for f in ["PCG_annotation/Cellular_RNA_analysis/03_Enrichment_results/ORA_KEGG_Upregulated.csv","PCG_annotation/Cellular_RNA_analysis/03_Enrichment_results/ORA_KEGG_Downregulated.csv","PCG_annotation/Cellular_RNA_analysis/03_Enrichment_results/GSEA_KEGG.csv"]:
    d=pd.read_csv(U+f); kn.update(dict(zip(d.ID,d.Description)))
old=pd.read_csv(U+"Differential_expression_analysis/03_Cellular_DE_analysis/mRNA_DE_analysis/EXPLORATORY_2026-09-25_contrast_comparison/06_ORA_all_sets.csv")
kn.update(dict(zip(old[old.db=="KEGG"].ID,old[old.db=="KEGG"].Description)))
out=[]
for n,d in sets.items():
    g=d.Feature
    a=enricher(g,go,uGO); a["db"]="GO"; a["Description"]=a.ID.map(lambda x: gn.get(x) or x)
    k=enricher(g,kg,uKG); k["db"]="KEGG"; k["Description"]=k.ID.map(lambda x: kn.get(x,x))
    print(f"\n== {n} genes={len(g)} GO={len(a)} KEGG={len(k)}")
    print(" GO:", "; ".join(f"{r.Description} ({r.Count}, {r._5:.1e})" for r in a.head(14).itertuples()))
    print(" KEGG:", "; ".join(f"{r.Description} ({r.Count}, {r._5:.1e})" for r in k.head(10).itertuples()))
    for x in (a,k):
        if len(x): x["set"]=n; out.append(x)
pd.concat(out).to_csv("ORA_stranded_all_sets.csv",index=False)
