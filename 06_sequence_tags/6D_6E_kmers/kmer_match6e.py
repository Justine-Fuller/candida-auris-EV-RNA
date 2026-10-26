#!/usr/bin/env python3
"""kmer_match6e.py -- Step 6E (v2, 2026-10-12), stage K2: GC matching (approved Q4 approach 1) + read-weighted 6-mer counts.
Usage: kmer_match6e.py <manifest.tsv (library, compartment EV|Cell_sRNA, group; header)> <K1 dir> <outdir> [ndraws=3]
Per group (3 EV vs 3 Cell sRNA) and RNA type (PCG_CDS, lncRNA, tRNA, 5S_rRNA):
- read GC = (G + C) / (A + C + G + T) of the read; 10 bins of 10% (bin = min(floor(10 * GC), 9)).
- EV target GC distribution = MEAN over the EV libraries (with >= 1 read of the type) of each library's bin proportions (each EV
  library weighs the same, whatever its depth).
- bins where EV has reads but some Cell library has none cannot be matched: they are dropped from BOTH sides for the matched
  counts (reads recorded in <group>.<type>.gc.tsv).
- each Cell sRNA library is DOWNSAMPLED without replacement, per GC bin, to the largest read total n whose bin split equals
  the EV target: n = floor(min_b N_b / p_b) over bins with p_b > 0; bin targets floor(n * p_b); reads drawn as a multivariate
  hypergeometric over that bin's distinct sequences. Fixed seed: SeedSequence([20261012, draw, group#, type#, library#]).
  Draw 1 = primary; draws 2.. = stability check only.
- 6-mers: every overlapping 6-mer (step 1) of every read; windows with any non-ACGT letter skipped; read-weighted (a read seen
  c times counts c times).
Outputs per group/type: <group>.<type>.counts.tsv.gz with columns word, gc_class (G+C count 0-6), then
  <lib>__all (all reads, no matching; EV and Cell), <lib>__m (EV: reads in matchable bins; Cell: matched draw 1), <lib>__m<d>
  (Cell draws d >= 2); and <group>.<type>.gc.tsv (reads per library per GC bin: all, matchable, matched per draw)."""
import sys, os, gzip, json, csv, itertools, collections
import numpy as np
MAN, K1, OUT = sys.argv[1:4]; ND = int(sys.argv[4]) if len(sys.argv) > 4 else 3
os.makedirs(OUT, exist_ok=False)
TYPES = ["PCG_CDS", "lncRNA", "tRNA", "5S_rRNA"]
WORDS = ["".join(p) for p in itertools.product("ACGT", repeat=6)]; WI = {w: i for i, w in enumerate(WORDS)}
GCC = np.array([sum(c in "GC" for c in w) for w in WORDS])
def gcbin(s):
    a = sum(s.count(c) for c in "ACGT")
    return None if a == 0 else min(int(10 * (s.count("G") + s.count("C")) / a), 9)
def kcount(pairs):
    v = np.zeros(4096, dtype=np.int64)
    for s, c in pairs:
        if c == 0: continue
        for i in range(len(s) - 5):
            j = WI.get(s[i:i + 6])
            if j is not None: v[j] += c
    return v
man = list(csv.DictReader(open(MAN), delimiter="\t"))
groups = list(dict.fromkeys(r["group"] for r in man))
cache = {}
def seqs(lib):
    return cache[lib]
summ = []
for gi, G in enumerate(groups):
    ev = [r["library"] for r in man if r["group"] == G and r["compartment"] == "EV"]
    ce = [r["library"] for r in man if r["group"] == G and r["compartment"] == "Cell_sRNA"]
    assert len(ev) == 3 and len(ce) == 3, G
    cache.clear(); cache.update({L: json.load(gzip.open(f"{K1}/{L}.k1.json.gz", "rt"))["types"] for L in ev + ce})
    for ti, T in enumerate(TYPES):
        D = {}   # lib -> {bin: [(seq, count)]}
        for L in ev + ce:
            d = collections.defaultdict(list); nonacgt = 0
            for s, c in seqs(L)[T].items():
                b = gcbin(s)
                if b is None: nonacgt += c; continue
                d[b].append((s, c))
            D[L] = d
        N = {L: np.array([sum(c for _, c in D[L].get(b, [])) for b in range(10)]) for L in ev + ce}
        evl = [L for L in ev if N[L].sum() > 0]; cel = [L for L in ce if N[L].sum() > 0]
        cols = collections.OrderedDict(); gcrows = []
        for L in ev + ce: cols[f"{L}__all"] = kcount(p for b in range(10) for p in D[L].get(b, []))
        if evl and cel:
            P = np.mean([N[L] / N[L].sum() for L in evl], axis=0)
            keep = np.array([not (P[b] > 0 and any(N[L][b] == 0 for L in cel)) for b in range(10)])
            Pk = np.where(keep, P, 0.0); Pk = Pk / Pk.sum()
            for L in ev: cols[f"{L}__m"] = kcount(p for b in range(10) if keep[b] for p in D[L].get(b, []))
            for li, L in enumerate(ce):
                if L not in cel:
                    for d in range(1, ND + 1): cols[f"{L}__m" if d == 1 else f"{L}__m{d}"] = np.zeros(4096, dtype=np.int64)
                    continue
                n = int(np.floor(min(N[L][b] / Pk[b] for b in range(10) if Pk[b] > 0)))
                tgt = np.floor(n * Pk + 1e-9).astype(np.int64); assert all(tgt <= N[L])
                for d in range(1, ND + 1):
                    rng = np.random.default_rng(np.random.SeedSequence([20261012, d, gi, ti, li]))
                    pairs = []
                    for b in range(10):
                        if tgt[b] == 0: continue
                        sq = [s for s, _ in D[L][b]]; cn = np.array([c for _, c in D[L][b]], dtype=np.int64)
                        draw = rng.multivariate_hypergeometric(cn, int(tgt[b]), method="marginals" if cn.size > 1000 else "count")
                        pairs += [(s, int(k)) for s, k in zip(sq, draw) if k]
                        if d == 1: gcrows.append((L, b, "matched_draw1", int(draw.sum())))
                        elif d == 2: gcrows.append((L, b, "matched_draw2", int(draw.sum())))
                    cols[f"{L}__m" if d == 1 else f"{L}__m{d}"] = kcount(pairs)
            for L in ev + ce:
                for b in range(10):
                    gcrows.append((L, b, "all", int(N[L][b]))); gcrows.append((L, b, "matchable" if keep[b] else "dropped_unmatchable_bin", int(N[L][b])))
            gcrows += [("EV_target", b, "proportion", float(Pk[b])) for b in range(10)]
        with gzip.open(f"{OUT}/{G}.{T}.counts.tsv.gz", "wt") as w:
            w.write("word\tgc_class\t" + "\t".join(cols) + "\n")
            M = np.column_stack(list(cols.values()))
            for i, wd in enumerate(WORDS): w.write(f"{wd}\t{GCC[i]}\t" + "\t".join(map(str, M[i])) + "\n")
        with open(f"{OUT}/{G}.{T}.gc.tsv", "w") as w:
            w.write("library\tgc_bin\tset\tvalue\n"); [w.write("\t".join(map(str, r)) + "\n") for r in gcrows]
        for L in ev + ce:
            row = dict(group=G, type=T, library=L, compartment="EV" if L in ev else "Cell_sRNA", reads=int(N[L].sum()),
                       kmers_all=int(cols[f"{L}__all"].sum()), kmers_matched=int(cols[f"{L}__m"].sum()) if f"{L}__m" in cols else "NA")
            summ.append(row); print(json.dumps(row))
keys = list(summ[0])
with open(f"{OUT}/K2_summary.tsv", "w") as w:
    w.write("\t".join(keys) + "\n"); [w.write("\t".join(str(r[k]) for k in keys) + "\n") for r in summ]
print("K2_DONE")
