#!/usr/bin/env python3
"""Abundance-weighted pies for Results.docx (Justine 2026-10-12: weight by TPM; cells = Cell sRNA library).
Abundance = TPM WITH rRNA (registered), per library, averaged over the 12 libraries of a compartment (each library sums to 1e6).
Units: PCG CDS (5,330 incl. the 26 V3 mRNAs), lncRNA, 43 tRNA families, pooled 5S, and the rDNA array as ONE unsplit rRNA unit
(read-census rDNA reads / 8,130 nt, as in fl_classify r2/r3; 18S/5.8S/25S are not annotated). Reads outside these units
(UTRs, antisense, intergenic, StringTie-only) are not in any TPM.
Pies 1/3: TPM share in full-length RNAs vs fragment-only RNAs (incl. RNAs below detection, ~0.1%) vs unsplit rRNA (not classified).
Pies 2/4: TPM share by RNA type. r2 (2026-10-12, Justine): 5S merged into one rRNA slice with the unsplit rDNA array.
Usage: plot_abundance_pies.py <tpm_with_rRNA_per_library.tsv.gz> <FL_r3 unit_table.tsv.gz> <outdir>"""
import csv, sys, os, gzip, collections, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
plt.rcParams["font.family"] = "Arial"
TP, UT, OUT = sys.argv[1:4]; os.makedirs(OUT, exist_ok=True)
st = {}
for x in csv.DictReader(gzip.open(UT, "rt"), delimiter="\t"):
    st[(x["compartment"], x["unit"])] = "full_length" if x["full_length_groups_TPM_with_rRNA"] else "fragment"
cls = collections.defaultdict(collections.Counter); sta = collections.defaultdict(collections.Counter); libs = collections.defaultdict(set)
for x in csv.DictReader(gzip.open(TP, "rt"), delimiter="\t"):
    c = x["compartment"]
    if c not in ("EV", "Cell_sRNA"): continue
    libs[c].add(x["library"]); v = float(x["TPM_with_rRNA"]); k = x["class"]; cls[c]["rRNA_all" if k in ("rRNA_rDNA_array", "rRNA_5S") else k] += v
    sta[c]["rRNA" if k == "rRNA_rDNA_array" else st[(c, x["unit"])]] += v
# r3 (2026-10-12, Justine): one blue scheme for all pies (dark blue -> lighter shades); plain labels
TYPES = [("tRNA_family", "tRNA", "#08306b"), ("rRNA_all", "rRNA (18S / 5.8S / 25S + 5S)", "#2171b5"),
         ("PCG_CDS", "mRNA (protein-coding)", "#6baed6"), ("lncRNA", "lncRNA", "#c6dbef")]
STAT = [("full_length", "Full-length", "#08306b"), ("fragment", "Fragment", "#6baed6"),
        ("rRNA", "rRNA 18S / 5.8S / 25S (not classified)", "#c6dbef")]
def pie(vals, cols, labels, title, fn):
    fig, a = plt.subplots(figsize=(5.2, 6.4))
    a.pie(vals, colors=cols, startangle=90, counterclock=False, wedgeprops=dict(edgecolor="white", linewidth=1.5))
    a.legend(handles=[Patch(color=c, label=l) for c, l in zip(cols, labels)], loc="upper center", bbox_to_anchor=(0.5, 0.02),
             frameon=False, fontsize=11.5, handlelength=1.2)
    a.set_title(title, fontsize=13.5, fontweight="bold"); plt.tight_layout()
    plt.savefig(f"{OUT}/{fn}", dpi=300, bbox_inches="tight"); plt.close(fig)
rows = []
for k, (c, name, tag) in enumerate((("Cell_sRNA", "Cells (sRNA library)", "cells_sRNA"), ("EV", "EVs", "EV"))):
    n = len(libs[c]); assert n == 12
    s = {key: v / n / 1e4 for key, v in sta[c].items()}; t = {key: v / n / 1e4 for key, v in cls[c].items()}
    assert abs(sum(s.values()) - 100) < 1e-6 and abs(sum(t.values()) - 100) < 1e-6
    pie([s[a] for a, _, _ in STAT], [col for _, _, col in STAT], [f"{l}  {s[a]:.1f}%" for a, l, _ in STAT],
        f"{name}\nabundance: full-length vs fragmented", f"pie{2 * k + 1}_{tag}_abundance_full_length_vs_fragment.png")
    pie([t[a] for a, _, _ in TYPES], [col for _, _, col in TYPES], [f"{l}  {t[a]:.1f}%" for a, l, _ in TYPES],
        f"{name}\nabundance by RNA type", f"pie{2 * k + 2}_{tag}_abundance_RNA_types.png")
    rows += [(c, "status", a, round(s[a], 3)) for a, _, _ in STAT] + [(c, "type", a, round(t[a], 3)) for a, _, _ in TYPES]
with open(f"{OUT}/abundance_pies_values.tsv", "w") as w:
    w.write("compartment\tpie\tslice\tmean_TPM_share_percent\n"); [w.write("\t".join(map(str, r)) + "\n") for r in rows]
print("\n".join("\t".join(map(str, r)) for r in rows))
