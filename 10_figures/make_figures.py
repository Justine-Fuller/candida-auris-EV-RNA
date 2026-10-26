"""Assemble main thesis figures (Thesis_v5). Run from the CandidaEV folder on the laptop.
Usage: python3 make_figures.py <CandidaEV root> <project-files thesis_v5 dir with inputs> <out dir>"""
import sys, os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont

root, extra, out = sys.argv[1:4]
os.makedirs(out, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
P = lambda *a: os.path.join(root, *a)

def font(size):
    for f in ["/System/Library/Fonts/Supplemental/Arial Bold.ttf", "/Library/Fonts/Arial Bold.ttf",
              "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"]:
        if os.path.exists(f):
            return ImageFont.truetype(f, size)
    return ImageFont.load_default()

def load(path, crop_top=0, crop_bottom=0):
    im = Image.open(path).convert("RGB")
    w, h = im.size
    return im.crop((0, crop_top, w, h - crop_bottom))

def scale_w(im, w):
    return im.resize((w, int(im.size[1] * w / im.size[0])), Image.LANCZOS)

def stack(panels, letters, width=3000, gap=60, name="fig.png"):
    """panels: list of rows; each row a list of images placed side by side."""
    rows = []
    for row in panels:
        n = len(row); pw = (width - gap * (n - 1)) // n
        ims = [scale_w(im, pw) for im in row]
        h = max(i.size[1] for i in ims)
        canvas = Image.new("RGB", (width, h + 70), "white")
        x = 0
        for im in ims:
            canvas.paste(im, (x, 70)); x += pw + gap
        rows.append((canvas, n, pw))
    H = sum(r[0].size[1] for r in rows) + gap * (len(rows) - 1)
    fig = Image.new("RGB", (width, H), "white")
    d = ImageDraw.Draw(fig); f = font(64); y = 0; k = 0
    for canvas, n, pw in rows:
        fig.paste(canvas, (0, y))
        for j in range(n):
            d.text((j * (pw + gap) + 5, y + 2), letters[k], fill="black", font=f); k += 1
        y += canvas.size[1] + gap
    fig.save(os.path.join(out, name), dpi=(300, 300))

# ---- Figure 1: caspofungin response (already A-C) ----
Image.open(P("Differential_expression_analysis/00_Summary_figures/Cell_mRNA_STRANDED_volcano_PCA_2026-10-14.png")).save(os.path.join(out, "Fig1_caspofungin_response.png"), dpi=(300, 300))

# ---- Figure 2: what is in EVs vs cells (new compact bars) ----
col = {"rRNA": "#9e9e9e", "tRNA": "#2a6fdb", "mRNA (protein-coding)": "#e8692f", "lncRNA": "#1aa37a",
       "Antisense to a gene": "#e8a400", "Other": "#d9d9d9",
       "Full-length": "#0b2e66", "Fragment": "#7fb3e6", "rRNA 18S/5.8S/25S (not classified)": "#cfcfcf"}
panels = [
 ("A  Sequencing reads by RNA type", ["rRNA", "tRNA", "mRNA (protein-coding)", "lncRNA", "Antisense to a gene", "Other"],
  {"Cells": [89.9, 5.1, 1.6, 1.1, 1.8, 0.6], "EVs": [91.9, 2.0, 2.5, 1.2, 1.8, 0.6]}),
 ("B  Estimated molecules (TPM) by RNA type", ["tRNA", "rRNA", "mRNA (protein-coding)", "lncRNA"],
  {"Cells": [67.6, 27.7, 2.6, 2.2], "EVs": [58.4, 32.4, 4.6, 4.7]}),
 ("C  Estimated molecules (TPM): whole or in pieces", ["Full-length", "Fragment", "rRNA 18S/5.8S/25S (not classified)"],
  {"Cells": [76.5, 4.6, 18.8], "EVs": [63.7, 8.7, 27.6]}),
]
fig, axes = plt.subplots(3, 1, figsize=(7.0, 5.6))
for ax, (title, cats, data) in zip(axes, panels):
    for yi, comp in enumerate(["EVs", "Cells"]):
        left = 0
        for c, v in zip(cats, data[comp]):
            ax.barh(yi, v, left=left, color=col[c], edgecolor="white", height=0.7)
            if v >= 4:
                ax.text(left + v / 2, yi, f"{v:.1f}%", ha="center", va="center", fontsize=7.5,
                        color="white" if c in ("tRNA", "Full-length", "rRNA") else "black")
            left += v
    ax.set_yticks([0, 1]); ax.set_yticklabels(["EVs", "Cells"]); ax.set_xlim(0, 100)
    ax.set_title(title, loc="left", fontsize=9.5, fontweight="bold")
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=col[c]) for c in cats], labels=cats, ncol=3,
              fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(0, -0.28))
fig.supxlabel("% of total (mean of 12 libraries per compartment)", fontsize=9, y=0.005)
plt.tight_layout(h_pad=2.6, rect=[0, 0.08, 1, 1])
plt.savefig(os.path.join(out, "Fig2_EV_RNA_content.png"), dpi=300); plt.close()

# ---- Figure 3: tRNA families ----
Image.open(P("v2_tRNA_2026-10-09/figures/tRNA_family_EV_vs_Cell_4groups_r1_with_dagger.png")).save(os.path.join(out, "Fig3_tRNA_families.png"), dpi=(300, 300))

# ---- Figure 4: EV-selective fragments (categories + TUB2 example) ----
catA = load(os.path.join(extra, "panel_A_categories_r2.png"), crop_top=90)
covB = load(P("v2_bin100_option1_DRAFT_2026-10-08/figures_ev_vs_cell/panel_D_example_coverage_r2_means_only.png"))  # 2026-10-23: replicate lines removed
stack([[catA], [covB]], "AB", name="Fig4_EV_selective_fragments.png")

# ---- Figure 5: codon ----
Image.open(os.path.join(extra, "codon_r3_figure.png")).save(os.path.join(out, "Fig5_codon_human_genes.png"), dpi=(300, 300))

# ---- Figure 6: host side ----
ddx = load(P("Differential splicing analysis/04_Figures/DDX21_splice_junction_usage_compact_PowerPoint.png"))
rbp = load(os.path.join(extra, "Fig1_RBP_binding_on_splicing_hits.png"))
mot = pd.read_csv(os.path.join(extra, "motif_presence_per_motif.tsv"), sep="\t")
enr = pd.read_csv(os.path.join(extra, "motif_presence_per_motif_EV_enriched.tsv"), sep="\t")
m = mot.merge(enr[["motif_id", "EV_enriched_RNAs_with_site"]], on="motif_id", how="left").fillna(0)
n_all, n_enr = 7152, 327
m["pct_enr"] = 100 * m.EV_enriched_RNAs_with_site / n_enr
m["pct_other"] = 100 * (m.EV_RNAs_with_site - m.EV_enriched_RNAs_with_site) / (n_all - n_enr)
fig, ax = plt.subplots(figsize=(3.6, 3.4))
ax.plot([0, 60], [0, 60], ls="--", color="grey", lw=1)
ax.scatter(m.pct_other, m.pct_enr, s=16, color="#2a6fdb", edgecolor="white", lw=0.4)
offs = [(8, 2), (10, -8), (12, -20)]
for (_, r), o in zip(m.nlargest(3, "pct_enr").iterrows(), offs):
    ax.annotate(r.rbp, (r.pct_other, r.pct_enr), fontsize=7, xytext=o, textcoords="offset points",
                arrowprops=dict(arrowstyle="-", lw=0.5, color="grey"))
ax.set_xlabel("% of other EV RNAs with a site\n(n = 6,825)"); ax.set_ylabel("% of EV-enriched RNAs with a site\n(n = 327)")
ax.set_xlim(0, 60); ax.set_ylim(0, 60)
ax.set_title("Each dot = one human RBP motif (75)", fontsize=8, loc="left")
plt.tight_layout(); motp = "/tmp/_tmp_motif_scatter.png"; plt.savefig(motp, dpi=300); plt.close()
stack([[ddx, Image.open(motp).convert("RGB")], [rbp]], "ABC", name="Fig6_host_splicing_RBPs.png")
print("done", sorted(os.listdir(out)))
# 2026-10-26: thesis Fig 6 replaced by figures_r2/Fig6_host_splicing_RBPs_r2.png (panel A = compact DDX21 built from Poster_CR media/image10.png; A:B widths 1840:1100; B and C cropped from the previous Fig 6).
