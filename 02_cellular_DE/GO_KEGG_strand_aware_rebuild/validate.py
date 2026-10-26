import pandas as pd, numpy as np, sys
sys.path.insert(0,".")
from ora import load_maps, enricher
go,kg=load_maps("emapper.tsv")
t=pd.read_csv("./old_contrasts.csv.gz")
t=t[t.Feature.str.startswith("B9J08_")]
uGO=sorted(set(t.Feature)&set(go.gene)); uKG=sorted(set(t.Feature)&set(kg.gene))
s=lambda p: p.notna()&(p<0.05)
sets={"B8441_up":t.Feature[s(t.padj_B8441)&(t.LFC_B8441>0)],"B8441_down":t.Feature[s(t.padj_B8441)&(t.LFC_B8441<0)],
"MMC1_up":t.Feature[s(t.padj_MMC1)&(t.LFC_MMC1>0)],"MMC1_down":t.Feature[s(t.padj_MMC1)&(t.LFC_MMC1<0)],
"Conserved_up":t.Feature[s(t.padj_B8441)&s(t.padj_MMC1)&(t.LFC_B8441>0)&(t.LFC_MMC1>0)],
"Conserved_down":t.Feature[s(t.padj_B8441)&s(t.padj_MMC1)&(t.LFC_B8441<0)&(t.LFC_MMC1<0)],
"Opp_up":t.Feature[s(t.padj_B8441)&s(t.padj_MMC1)&(t.LFC_B8441>0)&(t.LFC_MMC1<0)],
"Opp_down":t.Feature[s(t.padj_B8441)&s(t.padj_MMC1)&(t.LFC_B8441<0)&(t.LFC_MMC1>0)]}
print("universe",len(uGO),len(uKG))
for n,g in sets.items():
    a=enricher(g,go,uGO); b=enricher(g,kg,uKG)
    print(n,len(g),len(a),len(b))
