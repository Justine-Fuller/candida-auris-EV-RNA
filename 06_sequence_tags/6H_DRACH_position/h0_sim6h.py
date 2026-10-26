#!/usr/bin/env python3
"""h0_sim6h.py -- Step 6H synthetic libraries in H1 format (gates). 3 EV + 3 Cell per group, 30,000 reads each, length 18-35,
per-read GC ~ Beta (mean 0.45), position class 0-9 (0 = last CDS bin). In class 0, 40% of reads carry GGACT (a DRACH-rich gene end).
  SIM_posconfound: EV 40% of reads in class 0, Cell 10% (rest uniform over 1-9); nothing else differs
                   -> DRACH must be EV-higher WITHOUT position matching and NOT called WITH position matching.
  SIM_tag:         EV and Cell both 10% in class 0; EV reads in every class carry GGACT in an extra 5% -> DRACH called WITH matching.
  SIM_null:        EV = Cell                                                                       -> no call.
r2 (2026-10-15, Justine option a): optional 2nd argument "nojitter" writes ONLY a null group SIM_null_nojitter (no per-library GC
wobble: every library uses GC mean 0.45 exactly); without it the output is identical to r1 (same seed, same groups).
Usage: h0_sim6h.py <outdir> [nojitter]  -> <outdir>/H1/<lib>.h1.json.gz + <outdir>/manifest.tsv"""
import sys, os, gzip, json, collections
import numpy as np
OUT = sys.argv[1]; NOJ = len(sys.argv) > 2 and sys.argv[2] == "nojitter"; os.makedirs(f"{OUT}/H1", exist_ok=False)
rng = np.random.default_rng(20261015 + (1 if NOJ else 0))
def seq(L, gc):
    u = rng.random(L); return "".join(np.where(u < gc / 2, "G", np.where(u < gc, "C", np.where(u < gc + (1 - gc) / 2, "A", "T"))))
def lib(n, p0, extra):
    c = collections.Counter(); m = 0.45 if NOJ else 0.45 + rng.normal(0, 0.01)
    for _ in range(n):
        cl = 0 if rng.random() < p0 else int(rng.integers(1, 10)); L = int(rng.integers(18, 36)); s = seq(L, rng.beta(m * 40, (1 - m) * 40))
        if (cl == 0 and rng.random() < 0.4) or rng.random() < extra:
            i = int(rng.integers(0, L - 4)); s = s[:i] + "GGACT" + s[i + 5:]
        c[f"{s}\t{cl}"] += 1
    return c
man = ["library\tcompartment\tgroup"]
GROUPS = [("SIM_null_nojitter", 0.1, 0.1, 0.0)] if NOJ else [("SIM_posconfound", 0.4, 0.1, 0.0), ("SIM_tag", 0.1, 0.1, 0.05), ("SIM_null", 0.1, 0.1, 0.0)]
for G, evp, cep, extra in GROUPS:
    for comp, p0, ex in (("EV", evp, extra), ("Cell_sRNA", cep, 0.0)):
        for i in (1, 2, 3):
            L = f"{G}_{comp}_{i}"; man.append(f"{L}\t{comp}\t{G}")
            with gzip.open(f"{OUT}/H1/{L}.h1.json.gz", "wt") as w: json.dump(dict(library=L, accounting={}, reads=dict(lib(30000, p0, ex))), w)
open(f"{OUT}/manifest.tsv", "w").write("\n".join(man) + "\n"); print("H0_SIM_WRITTEN")
