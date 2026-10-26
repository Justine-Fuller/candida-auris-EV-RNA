#!/usr/bin/env python3
"""sim6e.py -- Step 6E gates (2026-10-12): synthetic K1-format libraries for the planted-word, GC-shift and null checks.
Usage: sim6e.py <outdir>   -> <outdir>/K1/<lib>.k1.json.gz (PCG_CDS only; other types empty) + <outdir>/manifest.tsv
Groups (3 EV + 3 Cell each, 30,000 reads per library, read length 18-35, per-read GC drawn from a Beta distribution,
G/C and A/T equally split, small per-library jitter of the GC mean (sd 0.01) so replicates differ):
  SIM_planted : EV and Cell same GC (mean 0.45); 3% of EV reads carry GGACTA at a random position  -> GGACTA must be called EV.
  SIM_gcshift : EV GC mean 0.52, Cell GC mean 0.42, no word planted   -> after GC matching no word may be called.
  SIM_null    : EV and Cell identical                               -> no word may be called."""
import sys, os, gzip, json, collections
import numpy as np
OUT = sys.argv[1]; os.makedirs(f"{OUT}/K1", exist_ok=False)
rng = np.random.default_rng(20261012)
def lib(n, gcmean, plant=None, frac=0.0, conc=40.0):
    c = collections.Counter()
    m = min(max(gcmean + rng.normal(0, 0.01), 0.05), 0.95)
    gcs = rng.beta(m * conc, (1 - m) * conc, n); lens = rng.integers(18, 36, n)
    for gc, L in zip(gcs, lens):
        u = rng.random(L); s = np.where(u < gc / 2, "G", np.where(u < gc, "C", np.where(u < gc + (1 - gc) / 2, "A", "T")))
        s = "".join(s)
        if plant and rng.random() < frac:
            i = int(rng.integers(0, L - 5)); s = s[:i] + plant + s[i + 6:]
        c[s] += 1
    return c
man = ["library\tcompartment\tgroup"]
for G, evgc, cegc, plant in (("SIM_planted", 0.45, 0.45, "GGACTA"), ("SIM_gcshift", 0.52, 0.42, None), ("SIM_null", 0.45, 0.45, None)):
    for comp, gc in (("EV", evgc), ("Cell_sRNA", cegc)):
        for i in (1, 2, 3):
            L = f"{G}_{comp}_{i}"; man.append(f"{L}\t{comp}\t{G}")
            t = lib(30000, gc, plant if comp == "EV" else None, 0.03)
            with gzip.open(f"{OUT}/K1/{L}.k1.json.gz", "wt") as w:
                json.dump(dict(library=L, accounting={}, types={"PCG_CDS": dict(t), "lncRNA": {}, "tRNA": {}, "5S_rRNA": {}}), w)
open(f"{OUT}/manifest.tsv", "w").write("\n".join(man) + "\n")
print("SIM_WRITTEN")
