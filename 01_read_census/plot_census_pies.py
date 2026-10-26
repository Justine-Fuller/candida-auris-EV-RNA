#!/usr/bin/env python3
"""Pie charts of RNA identity of reads (>= 18 nt, incl. multimappers; each read counted once) per compartment, pooled over
12 libraries, from the read census class_summary.tsv (r2: Cell sRNA, Cell poly(A) mRNA, EV). Left: rRNA vs all other RNA.
Right: the non-rRNA reads by class (legend values = % of ALL reads). Usage: plot_census_pies.py <class_summary.tsv>"""
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
OTHER_ALL = "#4d4d4d"
PANELS = (("Cell_sRNA", "Cellular sRNA"), ("Cell_polyA", "Cellular poly(A) mRNA"), ("EV", "Extracellular vesicles (EV)"))
fmt = lambda p: "<0.1%" if 0 < p < 0.05 else f"{p:.1f}%"
fig, ax = plt.subplots(len(PANELS), 2, figsize=(12, 15))
out = []
for i, (pool, title) in enumerate(PANELS):
    d = {r["class"]: int(r["reads_total"]) for r in rows if r["pool"] == pool}
    v = {g: sum(d[c] for c in cs) for g, cs in GROUPS.items()}; tot = sum(v.values()); pct = {g: 100 * n / tot for g, n in v.items()}
    out.append((title, tot, pct)); other = 100 - pct["rRNA"]
    a = ax[i][0]
    a.pie([pct["rRNA"], other], colors=[COL["rRNA"], OTHER_ALL], startangle=90, counterclock=False,
          wedgeprops=dict(edgecolor="white", linewidth=2))
    a.legend(handles=[Patch(color=COL["rRNA"], label=f"rRNA  {fmt(pct['rRNA'])}"), Patch(color=OTHER_ALL, label=f"All other RNA  {fmt(other)}")],
             loc="upper center", bbox_to_anchor=(0.5, 0.02), ncol=2, frameon=False, fontsize=12, handlelength=1.2)
    a.set_title(f"{title}: all reads\n{tot:,} reads, 12 libraries", fontsize=13, fontweight="bold")
    b = ax[i][1]; keys = [g for g in GROUPS if g != "rRNA"]
    b.pie([pct[g] for g in keys], colors=[COL[g] for g in keys], startangle=90, counterclock=False,
          wedgeprops=dict(edgecolor="white", linewidth=2, width=0.42))
    b.text(0, 0, f"non-rRNA\n{fmt(other)}\nof all reads", ha="center", va="center", fontsize=11, color="#333333")
    b.legend(handles=[Patch(color=COL[g], label=f"{g}  {fmt(pct[g])}") for g in keys], loc="center left", bbox_to_anchor=(0.98, 0.5),
             frameon=False, fontsize=11, handlelength=1.2, title="% of ALL reads", title_fontsize=11)
    b.set_title(f"{title}: the non-rRNA reads", fontsize=13, fontweight="bold")
fig.text(0.5, 0.012, "Reads ≥ 18 nt including multimappers, each read counted once; pooled over 3 replicates × 2 strains × 2 treatments.\n"
         "rRNA = Chr2 rDNA array (LNCRNA.1946/.1057 models) + 5S rRNA. Antisense = only on the opposite strand of an annotated gene.\n"
         "Poly(A) mRNA libraries are poly(A)-selected (rRNA and tRNA largely removed by design), so they are not a like-for-like comparison with the small-RNA libraries.",
         ha="center", fontsize=9.5, color="#444444")
plt.tight_layout(rect=(0, 0.05, 1, 1), h_pad=3.5, w_pad=1)
plt.savefig("figures/read_census_pies.png", dpi=300, bbox_inches="tight"); print("saved")
with open("figures/read_census_pies_values.tsv", "w") as f:
    f.write("compartment\treads\t" + "\t".join(GROUPS) + "\n")
    for t, tot, p in out: f.write(f"{t}\t{tot}\t" + "\t".join(f"{p[g]:.3f}" for g in GROUPS) + "\n")
print(open("figures/read_census_pies_values.tsv").read())
