#!/usr/bin/env python3
"""plot_ev_vs_cell_D.py -- 2026-10-13 draft panel D: read depth along the example transcript (from coverage_track.py on the HPC,
frozen read rules), EV vs Cell sRNA. Each library is scaled to the same total depth over the transcript (so only the SHAPE is
compared, as in the test); thin lines = replicates, thick = mean. Shaded = EV-selective bin (E1 main merged, q < 0.05).
Usage: plot_ev_vs_cell_D.py <group> <transcript_id> <label> <cov dir>"""
import sys, csv, gzip, os, glob, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Arial"
G, TX, LAB, CD = sys.argv[1:5]; H = os.path.dirname(os.path.abspath(__file__))
R = [r for r in csv.DictReader(gzip.open(f"{H}/../results_E1_main_merged/{G}.ev_vs_cell_bins.tsv.gz", "rt"), delimiter="\t") if r["transcript_id"] == TX]
fig, ax = plt.subplots(figsize=(8.6, 3.6))
for r in R:
    if r["category_q0.05"].endswith("_EV_share_higher"): ax.axvspan(int(r["bin_tx_start0"]), int(r["bin_tx_end"]), color="#fdd0a2", alpha=0.6, lw=0)
for pat, col, name in (("*_sRNA_*", "#9ecae1", "Cell sRNA"), ("*_EV_*", "#08306b", "EV")):
    tr = []
    for f in sorted(glob.glob(f"{CD}/{pat}.cov.tsv")):
        d = [int(l.split("\t")[1]) for l in open(f) if l[0].isdigit()]; t = sum(d); tr.append([100 * x / t for x in d])
        ax.plot(range(len(d)), tr[-1], color=col, lw=0.5, alpha=0.5)
    m = [sum(v) / len(v) for v in zip(*tr)]; ax.plot(range(len(m)), m, color=col, lw=1.8, label=f"{name} (mean of {len(tr)})")
for r in R: ax.axvline(int(r["bin_tx_start0"]), color="#dddddd", lw=0.5, zorder=0)
ax.set_xlim(0, len(m)); ax.set_xlabel("position along the coding sequence (nt, 5′ → 3′); grey lines = 100-nt bins", fontsize=10)
ax.set_ylabel("depth, % of the\ntranscript's total", fontsize=10)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.legend(frameon=False, fontsize=9, loc="upper left"); ax.set_title(f"D. {LAB}: read coverage, EV vs cell (shaded = EV-selective bin)", fontsize=10.5, fontweight="bold", loc="left")
plt.tight_layout(); plt.savefig(f"{H}/panel_D_example_coverage.png", dpi=300); print("DONE")
