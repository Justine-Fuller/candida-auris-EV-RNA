#!/usr/bin/env python3
"""Four single full-length pies for Results.docx (Justine 2026-10-12): cells = Cell poly(A) library, EVs; TPM with rRNA (registered).
(1) cells full-length vs fragment, (2) cells: RNA types of all detected RNAs, (3) EV full-length vs fragment, (4) EV RNA types of all detected RNAs.
r3 (2026-10-12, Justine): cells = BOTH cell libraries combined: an RNA is detected in cells if detected in the poly(A) OR the
sRNA library, full-length if full-length in either (each library judged on its own; reads are never pooled). Built from unit_table.
Usage: plot_fl_single_pies.py <FL_r3 unit_table.tsv.gz> <outdir>"""
import csv, sys, os, gzip, collections, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
plt.rcParams["font.family"] = "Arial"
S, OUT = sys.argv[1:3]; os.makedirs(OUT, exist_ok=True)
U = collections.defaultdict(dict)
for x in csv.DictReader(gzip.open(S, "rt"), delimiter="\t"):
    U[x["unit"]][x["compartment"]] = (x["class"], bool(x["detected_groups"]), bool(x["full_length_groups_TPM_with_rRNA"]))
def counts(comps):
    n = collections.Counter()
    for u, d in U.items():
        cls = d[comps[0]][0]; det = any(d[c][1] for c in comps); fl = any(d[c][2] for c in comps)
        if det: n[(cls, "full_length" if fl else "fragment_only")] += 1
    return n
CLS = [("PCG_CDS", "mRNA (protein-coding)"), ("tRNA_family", "tRNA (anticodon families)"), ("rRNA_5S", "5S rRNA"), ("lncRNA", "lncRNA")]
COL = {"tRNA_family": "#2a78d6", "PCG_CDS": "#eb6834", "lncRNA": "#1baf7a", "rRNA_5S": "#a6a6a6"}
def pie(vals, cols, labels, title, fn):
    fig, a = plt.subplots(figsize=(5.2, 6.4))
    a.pie(vals, colors=cols, startangle=90, counterclock=False, wedgeprops=dict(edgecolor="white", linewidth=1.5))
    a.legend(handles=[Patch(color=c, label=l) for c, l in zip(cols, labels)], loc="upper center", bbox_to_anchor=(0.5, 0.02),
             frameon=False, fontsize=11.5, handlelength=1.2)
    a.set_title(title, fontsize=13.5, fontweight="bold"); plt.tight_layout()
    plt.savefig(f"{OUT}/{fn}", dpi=300, bbox_inches="tight"); plt.close(fig)
for k, (comps, name, tag) in enumerate(((("Cell_polyA", "Cell_sRNA"), "Cells (poly(A) + sRNA libraries)", "cells_combined"), (("EV",), "EVs", "EV"))):
    n = counts(comps)
    fl = {c: n.get((c, "full_length"), 0) for c, _ in CLS}; fr = {c: n.get((c, "fragment_only"), 0) for c, _ in CLS}
    tf, tr = sum(fl.values()), sum(fr.values())
    pie([tf, tr], ["#5b3a9e", "#d9d9d9"], [f"Full-length  {tf:,} ({100 * tf / (tf + tr):.1f}%)", f"Fragment only  {tr:,} ({100 * tr / (tf + tr):.1f}%)"],
        f"{name}\n{tf + tr:,} detected RNAs", f"pie{2 * k + 1}_{tag}_full_length_vs_fragment.png")
    # r2 (2026-10-12, Justine): pies 2 and 4 show the RNA types of ALL detected RNAs; legend gives how many of each are full-length
    td = tf + tr; keep = [(c, l) for c, l in CLS if fl[c] + fr[c]]
    pc = lambda x: "<0.1%" if 0 < 100 * x / td < 0.1 else f"{100 * x / td:.1f}%"
    pie([fl[c] + fr[c] for c, _ in keep], [COL[c] for c, _ in keep],
        [f"{l}  {fl[c] + fr[c]:,} ({pc(fl[c] + fr[c])});  full-length {fl[c]:,}" for c, l in keep],
        f"{name}\nRNA types of all {td:,} detected RNAs", f"pie{2 * k + 2}_{tag}_detected_RNA_types.png")
    if any(fl[c] == 0 and fr[c] for c, _ in CLS):
        print(name, "no full-length:", {l: fr[c] for c, l in CLS if fl[c] == 0 and fr[c]})
print("saved")
