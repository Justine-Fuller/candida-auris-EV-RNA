#!/usr/bin/env python3
"""Two pie charts: RNA identity of reads (>= 18 nt, incl. multimappers; each read counted once) in cellular sRNA and EV
libraries (pooled over 12 libraries each), from the read census class_summary.tsv. Usage: plot_census_two_pies.py <class_summary.tsv>"""
import csv, sys, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
plt.rcParams["font.family"] = "Arial"
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
GROUPS = {"rRNA": ["rDNA_lncRNA", "rRNA_5S"], "tRNA": ["tRNA"], "Protein-coding (CDS + UTR)": ["PCG_CDS", "PCG_UTR"],
          "lncRNA": ["lncRNA"], "Antisense to an annotated gene": ["antisense_only"],
          "Other (intergenic, unannotated, other ncRNA)": ["StringTie_only", "intergenic", "OTHER_NC"]}
COL = {"rRNA": "#a6a6a6", "tRNA": "#2a78d6", "Protein-coding (CDS + UTR)": "#eb6834", "lncRNA": "#1baf7a",
       "Antisense to an annotated gene": "#eda100", "Other (intergenic, unannotated, other ncRNA)": "#d9d9d9"}
fig, ax = plt.subplots(1, 2, figsize=(12, 6.8))
for a, (pool, title) in zip(ax, (("Cell_sRNA", "Cells (sRNA library)"), ("EV", "Extracellular vesicles (EVs)"))):
    d = {r["class"]: int(r["reads_total"]) for r in rows if r["pool"] == pool}
    v = {g: sum(d[c] for c in cs) for g, cs in GROUPS.items()}; tot = sum(v.values())
    a.pie(list(v.values()), colors=[COL[g] for g in v], startangle=90, counterclock=False, wedgeprops=dict(edgecolor="white", linewidth=1.5))
    a.legend(handles=[Patch(color=COL[g], label=f"{g}  {100 * n / tot:.1f}%") for g, n in v.items()], loc="upper center",
             bbox_to_anchor=(0.5, 0.0), frameon=False, fontsize=11.5, handlelength=1.2)
    a.set_title(f"{title}\n{tot:,} reads, 12 libraries", fontsize=14, fontweight="bold")
fig.text(0.5, -0.12, "Reads ≥ 18 nt, including reads that map to several places, each counted once; pooled over 2 strains × 2 treatments × 3 replicates.\n"
         "rRNA = ribosomal DNA array on chromosome 2 + 5S rRNA. Antisense = only on the opposite strand of an annotated gene.", ha="center", fontsize=9.5, color="#444444")
plt.tight_layout(); plt.savefig("figures/read_census_two_pies.png", dpi=300, bbox_inches="tight"); print("saved")
