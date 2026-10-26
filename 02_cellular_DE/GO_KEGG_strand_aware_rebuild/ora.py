# Python re-implementation of clusterProfiler::enricher (DOSE enricher_internal) as used in
# EXPLORATORY_2026-09-25_contrast_comparison/02_ORA_by_contrast.R:
# pvalueCutoff 0.05 (on p and BH p.adjust), qvalueCutoff 0.2 (Storey), minGSSize 3, maxGSSize 500.
import numpy as np, pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

def load_maps(emapper):
    # emapper file: lines starting with "##" are comments; header line starts with "#query"
    rows = [l.rstrip("\n").split("\t") for l in open(emapper) if not l.startswith("##")]
    hdr = rows[0]; hdr[0] = "query"
    bg = pd.DataFrame(rows[1:], columns=hdr).replace("-", np.nan)
    bg["Gene_ID"] = bg["query"].str.replace(r"-.*", "", regex=True)
    bg = bg.drop_duplicates("Gene_ID")
    go = (bg.dropna(subset=["GOs"])[["Gene_ID", "GOs"]].assign(term=lambda d: d.GOs.str.split(","))
          .explode("term").query("term != ''")[["term", "Gene_ID"]].drop_duplicates())
    kg = (bg.dropna(subset=["KEGG_Pathway"])[["Gene_ID", "KEGG_Pathway"]]
          .assign(term=lambda d: d.KEGG_Pathway.str.split(",")).explode("term"))
    kg = kg[kg.term.str.startswith("map")][["term", "Gene_ID"]].drop_duplicates()
    return go.rename(columns={"Gene_ID": "gene"}), kg.rename(columns={"Gene_ID": "gene"})

def storey_q(p, lambdas=np.arange(0.05, 0.96, 0.05)):
    # qvalue::pi0est default "smoother" (smooth.spline df=3) approximated by a cubic polynomial fit
    p = np.asarray(p); m = len(p)
    if m < 2: return p.copy()
    pi0 = np.array([np.mean(p >= l) / (1 - l) for l in lambdas])
    try:
        fit = np.polyval(np.polyfit(lambdas, pi0, 3), lambdas); pi0_hat = min(max(fit[-1], 0), 1)
    except Exception:
        pi0_hat = 1.0
    if not np.isfinite(pi0_hat) or pi0_hat <= 0: pi0_hat = 1.0
    o = np.argsort(p); q = np.empty(m); run = 1.0
    for rank in range(m - 1, -1, -1):
        i = o[rank]; run = min(run, pi0_hat * m * p[i] / (rank + 1)); q[i] = run
    return np.minimum(q, 1)

def enricher(genes, t2g, universe, pcut=0.05, qcut=0.2, mn=3, mx=500):
    uni = set(universe)
    t2g = t2g[t2g.gene.isin(uni)]
    annotated = set(t2g.gene)
    N = len(annotated)                       # DOSE: N = number of universe genes with any annotation
    q = set(genes) & annotated; n = len(q)
    sizes = t2g.groupby("term").gene.nunique()
    sizes = sizes[(sizes >= mn) & (sizes <= mx)]
    hit = t2g[t2g.gene.isin(q) & t2g.term.isin(sizes.index)].groupby("term").gene.apply(list)
    if len(hit) == 0: return pd.DataFrame()
    k = hit.apply(len); M = sizes[hit.index]
    p = hypergeom.sf(k - 1, N, M, n)
    padj = multipletests(p, method="fdr_bh")[1]; qv = storey_q(p)
    r = pd.DataFrame({"ID": hit.index, "Count": k.values, "setSize": M.values, "pvalue": p, "p.adjust": padj,
                      "qvalue": qv, "geneID": ["/".join(sorted(g)) for g in hit.values], "n_in": n, "N": N})
    r = r[(r.pvalue < pcut) & (r["p.adjust"] < pcut) & (r.qvalue < qcut)]
    return r.sort_values("pvalue")
