#!/usr/bin/env python3
"""plot_ev_vs_cell_C.py -- 2026-10-13 draft panel C: each 100-nt bin's share of its transcript (bin count / transcript physical
reads) in the 3 EV and 3 Cell sRNA libraries of one group, for one example transcript (E1 main merged); EV-selective bin(s) marked.
Usage: plot_ev_vs_cell_C.py <group> <transcript_id> <label>"""
import sys, csv, gzip, os, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Arial"
G, TX, LAB = sys.argv[1:4]; H = os.path.dirname(os.path.abspath(__file__))
R = [r for r in csv.DictReader(gzip.open(f"{H}/../results_E1_main_merged/{G}.ev_vs_cell_bins.tsv.gz", "rt"), delimiter="\t") if r["transcript_id"] == TX]
R.sort(key=lambda r: int(r["bin_number"]))
libs = [k[len("count_"):] for k in R[0] if k.startswith("count_")]
ev = [l for l in libs if "_EV_" in l]; ce = [l for l in libs if "_EV_" not in l]
fig, ax = plt.subplots(figsize=(8.6, 3.8)); xs = list(range(1, len(R) + 1))
for i, r in enumerate(R):
    if r["category_q0.05"].endswith("_EV_share_higher"): ax.axvspan(i + 0.5, i + 1.5, color="#fdd0a2", alpha=0.6, lw=0, zorder=0)
for grp, col, off, name in ((ce, "#9ecae1", -0.15, "Cell sRNA"), (ev, "#08306b", 0.15, "EV")):
    for j, L in enumerate(grp):
        sh = [100 * float(r[f"count_{L}"]) / float(r[f"transcript_reads_{L}"]) for r in R]
        ax.scatter([x + off for x in xs], sh, s=18, color=col, zorder=3, label=name if j == 0 else None)
    m = [sum(100 * float(r[f"count_{L}"]) / float(r[f"transcript_reads_{L}"]) for L in grp) / len(grp) for r in R]
    ax.plot([x + off for x in xs], m, color=col, lw=1.2, zorder=2)
ax.axhline(100 / len(R), color="#777777", ls=":", lw=0.9, label="even spread (1 / number of bins)")
ax.set_xticks(xs); ax.set_xticklabels([str(x) for x in xs], fontsize=8.5); ax.set_xlabel("100-nt bin along the coding sequence (5′ → 3′)", fontsize=10)
ax.set_ylabel("% of the transcript's reads", fontsize=10)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.legend(frameon=False, fontsize=9, loc="upper left")
r = [x for x in R if x["category_q0.05"].endswith("_EV_share_higher")][0]
ax.set_title(f"C. {LAB}: share of each bin\n(shaded = EV-selective bin, q = {float(r['q_dataset']):.1g}; dots = replicates, line = mean)", fontsize=10.5, fontweight="bold", loc="left")
plt.tight_layout(); plt.savefig(f"{H}/panel_C_example_bin_shares.png", dpi=300); print("DONE", len(R), "bins")
