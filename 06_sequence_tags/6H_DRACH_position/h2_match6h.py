#!/usr/bin/env python3
"""h2_match6h.py -- Step 6H (v2, 2026-10-15), stage H2: like-with-like matching of Cell sRNA reads to EV reads (plan Q2 a), then DRACH
counts in the Step 6G "k5drach" format (one pooled DRACH row + the 1,006 other 5-mers; the rows partition all 5-mer windows), so the
unchanged g2_edger6g.R can test them.
The matching algorithm is the Step 6E K2 one (kmer_match6e.py) with the strata generalised from "10% GC bin" to:
  gc    = 10% GC bin of the read (the 6E/6G control, on the reads that have a position class);
  pos   = position class (100-nt CDS bin counted from the 3' end, from H1);
  posgc = position class x 10% GC bin.
EV target = mean of the EV libraries' stratum proportions; strata where EV has reads but some Cell library has none are dropped
from both sides; each Cell library is downsampled without replacement per stratum to the largest size with exactly the EV split
(multivariate hypergeometric over that stratum's distinct sequences); seeds SeedSequence([20261015, draw, group#, mode#, library#]);
draw 1 = primary. Reads with position class "NA" are excluded in every mode (so all modes use the same reads).
Usage: h2_match6h.py <manifest (library, compartment, group; header)> <H1 dir> <outdir> <ndraws> <gc|pos|posgc>
Outputs: <group>.PCG_CDS.counts.tsv.gz (word, gc_class, <lib>__all, <lib>__m, <lib>__m<d>), <group>.strata.tsv (reads per library
per stratum: all, matchable, matched draw 1; EV target), H2_summary.tsv."""
import sys, os, gzip, json, csv, itertools, collections
import numpy as np
MAN, H1, OUT, ND, MODE = sys.argv[1:6]; ND = int(ND); assert MODE in ("gc", "pos", "posgc"); os.makedirs(OUT, exist_ok=False)
MI = {"gc": 0, "pos": 1, "posgc": 2}[MODE]
DRACH = {a + b + "AC" + c for a in "AGT" for b in "AG" for c in "ACT"}; assert len(DRACH) == 18
ALL = ["".join(p) for p in itertools.product("ACGT", repeat=5)]
WORDS = ["DRACH"] + [w for w in ALL if w not in DRACH]; WI = {w: (0 if w in DRACH else WORDS.index(w)) for w in ALL}
GCC = ["NA"] + [str(sum(c in "GC" for c in w)) for w in WORDS[1:]]; NW = len(WORDS)
def gcbin(s):
    a = sum(s.count(c) for c in "ACGT")
    return None if a == 0 else min(int(10 * (s.count("G") + s.count("C")) / a), 9)
def kcount(pairs):
    v = np.zeros(NW, dtype=np.int64)
    for s, c in pairs:
        if c == 0: continue
        for i in range(len(s) - 4):
            j = WI.get(s[i:i + 5])
            if j is not None: v[j] += c
    return v
def key(s, cl):
    g = gcbin(s)
    if g is None: return None
    return g if MODE == "gc" else (cl if MODE == "pos" else f"{cl}|{g}")
man = list(csv.DictReader(open(MAN), delimiter="\t")); groups = list(dict.fromkeys(r["group"] for r in man)); summ = []
for gi, G in enumerate(groups):
    ev = [r["library"] for r in man if r["group"] == G and r["compartment"] == "EV"]
    ce = [r["library"] for r in man if r["group"] == G and r["compartment"] == "Cell_sRNA"]
    assert len(ev) == 3 and len(ce) == 3, G
    D = {}
    for L in ev + ce:
        reads = json.load(gzip.open(f"{H1}/{L}.h1.json.gz", "rt"))["reads"]; d = collections.defaultdict(list)
        for k, c in reads.items():
            s, cl = k.split("\t")
            if cl == "NA": continue
            kk = key(s, cl)
            if kk is not None: d[kk].append((s, c))
        D[L] = d
    keys = sorted({k for L in D for k in D[L]}, key=str); K = len(keys)
    N = {L: np.array([sum(c for _, c in D[L].get(k, [])) for k in keys]) for L in ev + ce}
    evl = [L for L in ev if N[L].sum() > 0]; cel = [L for L in ce if N[L].sum() > 0]
    cols = collections.OrderedDict(); srows = []; maxdev = 0.0; kept = {}
    for L in ev + ce: cols[f"{L}__all"] = kcount(p for k in keys for p in D[L].get(k, []))
    P = np.mean([N[L] / N[L].sum() for L in evl], axis=0)
    keep = np.array([not (P[i] > 0 and any(N[L][i] == 0 for L in cel)) for i in range(K)])
    Pk = np.where(keep, P, 0.0); Pk = Pk / Pk.sum()
    for L in ev: cols[f"{L}__m"] = kcount(p for i, k in enumerate(keys) if keep[i] for p in D[L].get(k, [])); kept[L] = int(N[L][keep].sum())
    for li, L in enumerate(ce):
        n = int(np.floor(min(N[L][i] / Pk[i] for i in range(K) if Pk[i] > 0)))
        tgt = np.floor(n * Pk + 1e-9).astype(np.int64); assert all(tgt <= N[L])
        for d in range(1, ND + 1):
            rng = np.random.default_rng(np.random.SeedSequence([20261015, d, gi, MI, li])); pairs = []; got = np.zeros(K, np.int64)
            for i, k in enumerate(keys):
                if tgt[i] == 0: continue
                sq = [s for s, _ in D[L][k]]; cn = np.array([c for _, c in D[L][k]], dtype=np.int64)
                draw = rng.multivariate_hypergeometric(cn, int(tgt[i]), method="marginals" if cn.size > 1000 else "count")
                pairs += [(s, int(x)) for s, x in zip(sq, draw) if x]; got[i] = draw.sum()
            assert (got == tgt).all()                                     # gate: each stratum drawn exactly to target
            if d == 1:
                kept[L] = int(got.sum())
                maxdev = max(maxdev, float(np.abs(got / got.sum() - Pk).max()))
                for i, k in enumerate(keys): srows.append((L, k, "matched_draw1", int(got[i])))
            cols[f"{L}__m" if d == 1 else f"{L}__m{d}"] = kcount(pairs)
    for L in ev + ce:
        for i, k in enumerate(keys): srows += [(L, k, "all", int(N[L][i])), (L, k, "matchable" if keep[i] else "dropped", int(N[L][i]))]
    srows += [("EV_target", k, "proportion", float(Pk[i])) for i, k in enumerate(keys)]
    with gzip.open(f"{OUT}/{G}.PCG_CDS.counts.tsv.gz", "wt") as w:
        w.write("word\tgc_class\t" + "\t".join(cols) + "\n"); M = np.column_stack(list(cols.values()))
        for i, wd in enumerate(WORDS): w.write(f"{wd}\t{GCC[i]}\t" + "\t".join(map(str, M[i])) + "\n")
    with open(f"{OUT}/{G}.strata.tsv", "w") as w:
        w.write("library\tstratum\tset\tvalue\n"); [w.write("\t".join(map(str, r)) + "\n") for r in srows]
    for L in ev + ce:
        summ.append(dict(group=G, mode=MODE, library=L, reads_with_position=int(N[L].sum()), strata=K, strata_dropped=int((~keep).sum()),
                         reads_kept_EV_or_drawn_Cell=kept[L],
                         windows_all=int(cols[f"{L}__all"].sum()), windows_matched=int(cols[f"{L}__m"].sum()),
                         max_proportion_deviation_cell_vs_EV_target=round(maxdev, 6)))
keys_ = list(summ[0])
with open(f"{OUT}/H2_summary.tsv", "w") as w:
    w.write("\t".join(keys_) + "\n"); [w.write("\t".join(str(r[k]) for k in keys_) + "\n") for r in summ]
print(open(f"{OUT}/H2_summary.tsv").read()); print("H2_DONE")
