#!/usr/bin/env python3
"""g4_stretch6g.py -- Step 6G (v2, 2026-10-14), stretch level (plan Q3 a, Q5 a): are GGAC / DRACH denser in EV-selective 100-nt bins
than in the other tested bins of the SAME transcripts?
Sets (identical to Step 6E Q6, k4_candidates_q6.py): PCG_CDS and lncRNA; EV-selective = bin100 E1 main merged category_q0.05 ending
"_EV_share_higher"; background = other TESTED bins of the same transcripts that are EV-selective in no group; "pooled" (a bin is
EV-selective if it is in any group) + each group.
Measure: rate per 100 nt = count * 100 / bin sequence length, for the 256 4-mers (GGAC one of them; ranking reference) and DRACH
(5-mer windows matching [AGT][AG]AC[ACT], pooled); bin sequence = RNA-strand genomic sequence of the bin's blocks (as Q6).
Statistic: mean rate in EV-selective bins minus mean rate in background bins.
Two versions: noGC = raw rates; GC = each unit's rates centred within 5 GC classes of the bins (quantile classes of bin GC over the
set's bins, class = ceil(rank * 5 / n), ties by bin id; as 6F). r2 (2026-10-14, Justine option a): r1 used a straight-line GC
residual and failed the GC-shift gate (1 call, CTCA); kept as g4_stretch6g_r1_linearGC.py.
p: labels shuffled WITHIN transcripts (each transcript keeps its number of EV-selective bins), two-sided, Besag-Clifford sequential
stopping (H exceedances, at most BMAX shuffles); BH over the 257 units of a set; call = q < 0.05.
Gates (mode sim): planted GGACT written once into 50% of the pooled PCG EV-selective bins -> GGAC and DRACH called (both versions);
GC shift (A/T -> G/C at 10% in EV-selective bins) -> 0 calls in the GC version. Mode real also runs the random-label null: NNULL
draws where, per transcript, the same number of fake EV-selective bins is drawn from its background bins (real EV-selective bins
removed); gate = median calls 0 per set and version.
Usage: g4_stretch6g.py <real|sim> <E1 dir> <bin100m annotation> <genome.fasta> <outdir> <seed> <H> <BMAX> <NNULL>"""
import sys, os, gzip, csv, itertools, collections, random
import numpy as np
MODE, E1, ANN, FA, OUT, SEED, H, BMAX, NNULL = sys.argv[1:10]; SEED, H, BMAX, NNULL = int(SEED), int(H), int(BMAX), int(NNULL)
assert MODE in ("real", "sim"); os.makedirs(OUT, exist_ok=False)
GROUPS = ["B8441_untreated", "B8441_caspofungin", "MMC1_untreated", "MMC1_caspofungin"]
K4 = ["".join(p) for p in itertools.product("ACGT", repeat=4)]; UNITS = K4 + ["DRACH"]; UI = {w: i for i, w in enumerate(K4)}
DRACH = {a + b + "AC" + c for a in "AGT" for b in "AG" for c in "ACT"}; assert len(DRACH) == 18
def counts(s):
    v = np.zeros(len(UNITS))
    for i in range(len(s) - 3):
        j = UI.get(s[i:i + 4])
        if j is not None: v[j] += 1
    for i in range(len(s) - 4):
        if s[i:i + 5] in DRACH: v[-1] += 1
    return v
assert counts("GGACTGGACA")[UI["GGAC"]] == 2 and counts("GGACTGGACA")[-1] == 2 and counts("GGACNGGACA")[-1] == 1
G2 = {}; name = None; buf = []
for l in open(FA):
    if l.startswith(">"):
        if name: G2[name] = "".join(buf)
        name = l[1:].split()[0]; buf = []
    else: buf.append(l.strip().upper())
G2[name] = "".join(buf)
RC = str.maketrans("ACGTN", "TGCAN")
ann = {r["bin_id"]: r for r in csv.DictReader(open(ANN), delimiter="\t")}
def bseq(b):
    a = ann[b]; parts = []
    for blk in a["genomic_blocks_5to3"].split(";"):
        s, e = map(int, blk.split("-")); x = G2[a["chrom"]][s:e]
        parts.append(x.translate(RC)[::-1] if a["strand"] == "-" else x)
    return "".join(parts)
st = {}
for G in GROUPS:
    d = {}
    for r in csv.DictReader(gzip.open(f"{E1}/{G}.ev_vs_cell_bins.tsv.gz", "rt"), delimiter="\t"):
        d[r["bin_id"]] = (r["transcript_id"], r["rna_class"], r["status"] == "tested", r["category_q0.05"].endswith("_EV_share_higher"))
    st[G] = d
def sets(groups, cls):   # identical to k4_candidates_q6.py
    fg = {b for G in groups for b, v in st[G].items() if v[3] and v[1] == cls}
    tx = {st[G][b][0] for G in groups for b in fg if b in st[G]}
    bg = {b for G in groups for b, v in st[G].items() if v[0] in tx and v[2] and v[1] == cls} - {b for G in GROUPS for b, v in st[G].items() if v[3]}
    return sorted(fg), sorted(bg)
def txof(b): return next(st[G][b][0] for G in GROUPS if b in st[G])
def test(R, lab, tx, seed):
    """R: bins x units rates; lab: 1 = EV-selective; tx: transcript per bin. Returns diff, p, q, perms used."""
    lab = np.asarray(lab, float); n1, n0 = lab.sum(), (1 - lab).sum()
    w = lab / n1 - (1 - lab) / n0; d0 = R.T @ w
    const = np.all(R == R[0], axis=0)
    grp = collections.defaultdict(list)
    for i, t in enumerate(tx): grp[t].append(i)
    idx = [np.array(v) for v in grp.values() if len(v) > 1 and 0 < lab[v].sum() < len(v)]
    ex = np.zeros(R.shape[1], int); L = np.zeros(R.shape[1], int); act = np.where(~const)[0]; rng = np.random.default_rng(seed)
    while act.size and L.max() < BMAX:
        wl = w.copy()
        for ix in idx: wl[ix] = w[ix][rng.permutation(ix.size)]
        d = R[:, act].T @ wl
        ex[act] += np.abs(d) >= np.abs(d0[act]) * (1 - 1e-9); L[act] += 1; act = act[ex[act] < H]
    p = np.where(const, 1.0, np.where(ex >= H, H / np.maximum(L, 1), (ex + 1) / (L + 1)))
    p = np.minimum(p, 1.0)
    o = np.argsort(p); m = p.size; q = np.empty(m); prev = 1.0
    for r in range(m, 0, -1): prev = min(prev, p[o[r - 1]] * m / r); q[o[r - 1]] = prev
    return d0, p, q, L
def rates(seqs): return np.array([counts(s) * 100.0 / len(s) for s in seqs])
def resid(R, gc, ids):
    n = len(gc); o = sorted(range(n), key=lambda i: (gc[i], ids[i])); cl = np.empty(n, int)
    for r, i in enumerate(o, 1): cl[i] = -(-r * 5 // n)
    out = R.copy()
    for c in np.unique(cl): out[cl == c] -= R[cl == c].mean(axis=0)
    return out
gcf = lambda s: (s.count("G") + s.count("C")) / len(s)
gates = []; out = []; summ = []
def run_set(label, cls, fgs, bgs, seqfg, seqbg, txs, seed, nnull):
    ids = list(fgs) + list(bgs); seqs = seqfg + seqbg; lab = [1] * len(seqfg) + [0] * len(seqbg); R = rates(seqs); gc = np.array([gcf(s) for s in seqs])
    res = {}
    for ver, RR in (("noGC", R), ("GC", resid(R, gc, ids))):
        d, p, q, L = test(RR, lab, txs, seed); res[ver] = (d, p, q)
        mE = R[np.array(lab) == 1].mean(axis=0); mB = R[np.array(lab) == 0].mean(axis=0)
        rk = {u: i + 1 for i, u in enumerate(sorted(K4, key=lambda u: (p[UNITS.index(u)], u)))}
        for i, u in enumerate(UNITS):
            out.append(dict(set=label, rna_class=cls, version=ver, unit=u, EV_selective_bins=len(seqfg), background_bins=len(seqbg),
                            mean_rate_per100_EV_selective=round(mE[i], 6), mean_rate_per100_background=round(mB[i], 6), diff_stat=d[i],
                            p=p[i], q=q[i], perms=L[i], direction="higher_in_EV_selective" if d[i] > 0 else ("higher_in_background" if d[i] < 0 else "none"),
                            rank_among_256_by_p=rk.get(u, "NA")))
        row = dict(set=label, rna_class=cls, version=ver, EV_selective_bins=len(seqfg), background_bins=len(seqbg),
                   transcripts=len(set(txs)), calls=int((q < 0.05).sum()), GGAC_q=q[UI["GGAC"]], GGAC_rank=rk["GGAC"], DRACH_q=q[-1])
        if nnull:
            # random-label null: per transcript, fake EV-selective bins drawn from its background bins
            bgi = collections.defaultdict(list); nfg = collections.Counter(txs[:len(seqfg)])
            for i in range(len(seqfg), len(seqs)): bgi[txs[i]].append(i)
            rr = random.Random(seed + 7); nl = []
            for b in range(nnull):
                pick = []
                for t, k in nfg.items():
                    if len(bgi[t]) > k: pick += rr.sample(bgi[t], k)
                if not pick: continue
                keep = list(range(len(seqfg), len(seqs))); lab2 = [1 if i in set(pick) else 0 for i in keep]
                if sum(lab2) == 0 or sum(lab2) == len(lab2): continue
                RR2 = RR[keep] if ver == "noGC" else resid(R[keep], gc[keep], [ids[i] for i in keep])
                _, _, q2, _ = test(RR2, lab2, [txs[i] for i in keep], seed + 1000 + b); nl.append(int((q2 < 0.05).sum()))
            row.update(null_draws=len(nl), null_median=float(np.median(nl)) if nl else "NA", null_max=max(nl) if nl else "NA",
                       null_draws_with_call=sum(x > 0 for x in nl))
            if nl: gates.append((f"G4_null_median_0_{cls}_{label}_{ver}", float(np.median(nl)) == 0, f"{sum(x > 0 for x in nl)}/{len(nl)} draws with a call; max {max(nl)}"))
        summ.append(row)
    return res
rng = random.Random(SEED)
for cls in ("PCG_CDS", "lncRNA"):
    for label, groups in [("pooled", GROUPS)] + [(G, [G]) for G in GROUPS]:
        if MODE == "sim" and not (cls == "PCG_CDS" and label == "pooled"): continue
        fg, bg = sets(groups, cls)
        if not fg or not bg: summ.append(dict(set=label, rna_class=cls, version="both", EV_selective_bins=len(fg), background_bins=len(bg), calls="not_testable")); continue
        sf = [bseq(b) for b in fg]; sb = [bseq(b) for b in bg]; txs = [txof(b) for b in fg + bg]
        if MODE == "real":
            run_set(label, cls, fg, bg, sf, sb, txs, SEED, NNULL)
        else:
            pick = set(rng.sample(range(len(sf)), round(0.5 * len(sf))))
            sp = [(s[:p] + "GGACT" + s[p + 5:]) if i in pick else s for i, s in enumerate(sf) for p in [rng.randrange(0, len(s) - 4)]]
            res = run_set("pooled_planted", cls, fg, bg, sp, sb, txs, SEED, 0)
            for ver in ("noGC", "GC"):
                q = res[ver][2]; d = res[ver][0]
                gates.append((f"G4_SIM_planted_GGAC_DRACH_called_{ver}", q[UI["GGAC"]] < 0.05 and q[-1] < 0.05 and d[UI["GGAC"]] > 0 and d[-1] > 0,
                              f"GGAC q={q[UI['GGAC']]:.3g}, DRACH q={q[-1]:.3g}"))
            sg = ["".join(rng.choice("GC") if ch in "AT" and rng.random() < 0.1 else ch for ch in s) for s in sf]
            res = run_set("pooled_gcshift", cls, fg, bg, sg, sb, txs, SEED, 0)
            n = int((res["GC"][2] < 0.05).sum())
            gates.append(("G4_SIM_gcshift_0_calls_GC_version", n == 0, f"GC version calls {n}; noGC version calls {int((res['noGC'][2] < 0.05).sum())}"))
def wr(path, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w") as w:
        w.write("\t".join(keys) + "\n"); [w.write("\t".join(str(r.get(k, "")) for k in keys) + "\n") for r in rows]
wr(f"{OUT}/stretch_units.tsv", out); wr(f"{OUT}/stretch_summary.tsv", summ)
with open(f"{OUT}/gates.tsv", "w") as w:
    w.write("gate\tresult\tdetail\n"); [w.write(f"{n}\t{'PASS' if ok else 'FAIL'}\t{d}\n") for n, ok, d in gates]
for r in summ: print(r)
for g in gates: print(*g)
ok = all(g[1] for g in gates); print("G4_ALL_GATES_PASS" if ok else "G4_GATE_FAILURE"); sys.exit(0 if ok else 1)
