#!/usr/bin/env python3
"""plot_ev_vs_cell_AB.py -- 2026-10-13 (Justine): figure drafts for Results "Are EV hotspots also enriched in cells?".
Input (read-only): ../results_E1_main_merged/<group>.ev_vs_cell_bins.tsv.gz and category_counts.tsv (q < 0.05).
A: EV hotspot bins per group split by outcome (share of all EV hotspots; n on top).
B: volcano per group of all tested bins: x = log2 fold enrichment of the observed share (EV / Cell), y = -log10 q (dataset BH);
   dark blue = EV hotspots with a significantly higher EV share (EV-selective), light blue = other bins with a higher EV share,
   orange = higher Cell share, grey = not significant. Bins with a Cell share of 0 are drawn at the right edge (triangles)."""
import csv, gzip, math, os, collections, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
plt.rcParams["font.family"] = "Arial"
H = os.path.dirname(os.path.abspath(__file__)); E1 = f"{H}/../results_E1_main_merged"
G = [("B8441_untreated", "B8441 untreated"), ("B8441_caspofungin", "B8441 caspofungin"), ("MMC1_untreated", "MMC1 untreated"), ("MMC1_caspofungin", "MMC1 caspofungin")]
CAT = [("EV share significantly higher\n(consistent with selective packaging)", ["EV_hotspot__EV_share_higher", "EV_hotspot__also_Cell_hotspot__EV_share_higher"], "#08306b"),
       ("Hotspot in both, no difference\n(consistent with abundance)", ["EV_hotspot__also_Cell_hotspot__no_significant_difference"], "#4292c6"),
       ("No Cell hotspot, no difference\n(inconclusive)", ["EV_hotspot__no_Cell_hotspot__no_significant_difference"], "#c6dbef"),
       ("Cell share significantly higher", ["EV_hotspot__Cell_share_higher"], "#e6550d"),
       ("Not testable (too few reads)", ["EV_hotspot__comparison_not_tested"], "#d9d9d9")]
cc = collections.Counter()
for r in csv.DictReader(open(f"{E1}/category_counts.tsv"), delimiter="\t"):
    if r["level"] == "0.05": cc[(r["group"], r["category"])] += int(r["bins"])
# ---- A
fig, ax = plt.subplots(figsize=(8.6, 4.6))
x = range(len(G)); bottom = [0] * len(G)
for lab, keys, col in CAT:
    n = [sum(cc[(g, k)] for k in keys) for g, _ in G]; tot = [sum(cc[(g, k)] for _, ks, _ in CAT for k in ks) for g, _ in G]
    v = [100 * a / b for a, b in zip(n, tot)]
    ax.bar(x, v, bottom=bottom, color=col, edgecolor="white", width=0.62, label=lab)
    for i, (vi, ni) in enumerate(zip(v, n)):
        if vi >= 6: ax.text(i, bottom[i] + vi / 2, str(ni), ha="center", va="center", fontsize=9, color="white" if col in ("#08306b", "#4292c6", "#e6550d") else "black")
    bottom = [b + vi for b, vi in zip(bottom, v)]
tot = [sum(cc[(g, k)] for _, ks, _ in CAT for k in ks) for g, _ in G]
for i, t in enumerate(tot): ax.text(i, 101.5, f"n = {t:,}", ha="center", va="bottom", fontsize=9.5)
ax.set_xticks(list(x)); ax.set_xticklabels([l for _, l in G], fontsize=10); ax.set_ylim(0, 108); ax.set_ylabel("% of EV hotspot bins", fontsize=10.5)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=9, handlelength=1.2, labelspacing=0.9)
ax.set_title("A. What happens to each EV hotspot in the direct EV-versus-Cell test (q < 0.05)", fontsize=11, fontweight="bold", loc="left")
plt.tight_layout(); plt.savefig(f"{H}/panel_A_categories.png", dpi=300); plt.close(fig)
# ---- B
fig, axs = plt.subplots(2, 2, figsize=(9.5, 8.2))
for ax, (g, lab) in zip(axs.flat, G):
    P = {"ns": [], "ev": [], "sel": [], "cell": []}; inf = {"ns": [], "ev": [], "sel": [], "cell": []}
    for r in csv.DictReader(gzip.open(f"{E1}/{g}.ev_vs_cell_bins.tsv.gz", "rt"), delimiter="\t"):
        if r["status"] != "tested": continue
        q = float(r["q_dataset"]); fe = r["fold_enrichment_EV_over_Cell"]; d = r["direction_observed"]
        y = -math.log10(max(q, 1e-300))
        k = "ns" if q >= 0.05 else ("sel" if d == "higher_in_EV" and r["category_q0.05"].endswith("_EV_share_higher") else "ev" if d == "higher_in_EV" else "cell" if d == "higher_in_Cell" else "ns")
        if fe in ("Inf", "inf"): inf[k].append(y)
        elif fe not in ("NA", "") and float(fe) > 0: P[k].append((math.log2(float(fe)), y))
    allx = [a for v in P.values() for a, _ in v]; xmax = max(abs(min(allx)), max(allx)) + 0.5
    for k, col, sz, z in (("ns", "#bdbdbd", 4, 1), ("cell", "#e6550d", 10, 2), ("ev", "#9ecae1", 10, 2), ("sel", "#08306b", 16, 3)):
        if P[k]: ax.scatter([a for a, _ in P[k]], [b for _, b in P[k]], s=sz, c=col, lw=0, zorder=z)
        if inf[k]: ax.scatter([xmax] * len(inf[k]), inf[k], s=sz + 6, c=col, marker=">", lw=0, zorder=z)
    ax.axhline(-math.log10(0.05), color="#555555", ls="--", lw=0.8); ax.axvline(0, color="#555555", lw=0.6)
    ax.set_xlim(-xmax - 0.3, xmax + 0.6); ax.set_title(lab, fontsize=10.5, fontweight="bold", loc="left")
    ax.set_xlabel("log$_2$ (EV share / Cell share) of the bin", fontsize=9.5); ax.set_ylabel("−log$_{10}$ q", fontsize=9.5)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    n = sum(len(v) for v in P.values()) + sum(len(v) for v in inf.values())
    ax.text(0.02, 0.97, f"{n:,} bins tested", transform=ax.transAxes, fontsize=8.5, va="top")
fig.legend(handles=[Line2D([], [], marker="o", ls="", color=c, markersize=6, label=l) for c, l in
                    (("#08306b", "EV hotspot, EV share higher (EV-selective)"), ("#9ecae1", "not an EV hotspot, EV share higher"),
                     ("#e6550d", "Cell share higher"), ("#bdbdbd", "not significant"))] +
           [Line2D([], [], color="#555555", ls="--", label="q = 0.05")], loc="lower center", ncol=3, frameon=False, fontsize=9)
fig.suptitle("Every tested bin: EV-versus-Cell share and significance (q < 0.05)", fontsize=11, fontweight="bold", x=0.02, ha="left")
plt.tight_layout(rect=(0, 0.07, 1, 0.97)); plt.savefig(f"{H}/panel_B_volcano.png", dpi=300); plt.close(fig)
print("DONE")
