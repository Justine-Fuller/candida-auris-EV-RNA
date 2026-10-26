#!/usr/bin/env python3
"""Full-length vs fragment pies from FL r3 (main full-length result from 2026-10-12, Justine). Registered pie criteria:
unit = transcript / tRNA anticodon family / pooled 5S rRNA; detected = >= 1 read in >= 2 of 3 replicates of >= 1 group;
full-length = TPM > 100 + coverage rule in >= 2 of 3 replicates of >= 1 group. Cells = Cell sRNA library (registered rule).
One figure per TPM version (the version is not chosen yet).
Usage: plot_fl_pies.py <FL_r3 summary.tsv> <outdir> [polyA]
r2 (2026-10-12): optional "polyA" draws cells from the Cell poly(A) library (Justine), TPM with rRNA only (registered 2026-10-12b)."""
import csv, sys, os, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
plt.rcParams["font.family"] = "Arial"
S, OUT = sys.argv[1:3]; PA = len(sys.argv) > 3 and sys.argv[3] == "polyA"; os.makedirs(OUT, exist_ok=True)
CELL = ("Cell_polyA", "Cells (poly(A) mRNA library)") if PA else ("Cell_sRNA", "Cells (sRNA library)")
rows = list(csv.DictReader(open(S), delimiter="\t"))
CLS = [("tRNA_family", "tRNA (anticodon families)"), ("PCG_CDS", "Protein-coding (CDS)"), ("lncRNA", "lncRNA"), ("rRNA_5S", "5S rRNA")]
COL = {"tRNA_family": "#2a78d6", "PCG_CDS": "#eb6834", "lncRNA": "#1baf7a", "rRNA_5S": "#a6a6a6"}
vals = []
VERS = (("with_rRNA", "TPM computed with rRNA"),) if PA else (("without_rRNA", "TPM computed without rRNA"), ("with_rRNA", "TPM computed with rRNA"))
for ver, vlab in VERS:
    fig, ax = plt.subplots(2, 2, figsize=(12, 11.5))
    for i, (comp, title) in enumerate((CELL, ("EV", "Extracellular vesicles (EVs)"))):
        n = {(r["class"], r["status"]): int(r["units"]) for r in rows if r["compartment"] == comp and r["TPM_version"] == ver}
        fl = {c: n.get((c, "full_length"), 0) for c, _ in CLS}; fr = {c: n.get((c, "fragment_only"), 0) for c, _ in CLS}
        tf, tr = sum(fl.values()), sum(fr.values())
        for c, lab in CLS: vals.append((ver, comp, c, fl[c], fr[c]))
        a = ax[i][0]
        a.pie([tf, tr], colors=["#5b3a9e", "#d9d9d9"], startangle=90, counterclock=False, wedgeprops=dict(edgecolor="white", linewidth=1.5))
        a.legend(handles=[Patch(color="#5b3a9e", label=f"Full-length  {tf:,} ({100 * tf / (tf + tr):.1f}%)"),
                          Patch(color="#d9d9d9", label=f"Fragment only  {tr:,} ({100 * tr / (tf + tr):.1f}%)")],
                 loc="upper center", bbox_to_anchor=(0.5, 0.0), frameon=False, fontsize=11.5, handlelength=1.2)
        a.set_title(f"{title}\n{tf + tr:,} detected RNAs", fontsize=14, fontweight="bold")
        a = ax[i][1]; keep = [(c, lab) for c, lab in CLS if fl[c]]
        a.pie([fl[c] for c, _ in keep], colors=[COL[c] for c, _ in keep], startangle=90, counterclock=False,
              wedgeprops=dict(edgecolor="white", linewidth=1.5))
        a.legend(handles=[Patch(color=COL[c], label=f"{lab}  {fl[c]:,} of {fl[c] + fr[c]:,} detected") for c, lab in CLS],
                 loc="upper center", bbox_to_anchor=(0.5, 0.0), frameon=False, fontsize=11.5, handlelength=1.2)
        a.set_title(f"{title}\nwhat the {tf:,} full-length RNAs are", fontsize=14, fontweight="bold")
    fig.text(0.5, -0.06, "Full-length = TPM > 100 and every 100-nt bin with mean depth ≥ 3× (tRNA: mean depth ≥ 3× over the mature tRNA "
             "plus ≥ 1 CCA read; 5S rRNA: 118 nt, bins 0–100 and 100–118), in ≥ 2 of 3 replicates of at least one strain × treatment group.\n"
             f"Detected = ≥ 1 read in ≥ 2 of 3 replicates of at least one group. 5S rRNA = all 13 copies pooled; tRNAs counted as 43 anticodon families. {vlab}.\n"
             "18S / 5.8S / 25S rRNA are not annotated in the genome and are not included.", ha="center", fontsize=9.5, color="#444444")
    plt.tight_layout(); plt.savefig(f"{OUT}/full_length_pies_FL_r3_{ver}{'_cells_polyA' if PA else ''}.png", dpi=300, bbox_inches="tight"); plt.close(fig)
with open(f"{OUT}/full_length_pies_FL_r3_values{'_cells_polyA' if PA else ''}.tsv", "w") as w:
    w.write("TPM_version\tcompartment\tclass\tfull_length\tfragment_only\n"); [w.write("\t".join(map(str, v)) + "\n") for v in vals]
print("saved")
