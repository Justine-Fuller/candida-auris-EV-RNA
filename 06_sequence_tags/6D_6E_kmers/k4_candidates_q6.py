#!/usr/bin/env python3
"""k4_candidates_q6.py -- Step 6E (v2, 2026-10-12), stage K4: tag candidates (Q5 a) + confirmation inside the same genes (Q6 a).
Usage: k4_candidates_q6.py <K3 dir> <E1_main_merged dir> <bin100m_master_annotation.tsv> <genome.fasta> <outdir>
Q5 a: a word is a tag candidate for an RNA type when it is significant (primary GC-matched test, q < 0.05) with the SAME observed
      direction in >= 3 of the 4 groups. Reported alongside: the GC-class cross-check, Cell draws 2 / 3 (stability), naive test.
Q6 a: for candidate words of PCG_CDS and lncRNA: EV-selective 100-nt bins (v2 bin100 E1 main merged layout; category at q < 0.05 =
      EV hotspot AND significantly higher EV share, i.e. category_q0.05 ending "_EV_share_higher", the Step 5A L2a definition)
      vs the other TESTED bins (status "tested") of the same transcripts. Word present = occurs >= 1 time in the bin's RNA-strand
      genomic sequence (blocks 5'->3', minus strand reverse-complemented). Two-sided Fisher exact test per word; BH over the
      candidate words of that type. Done per group and pooled (a bin counts as EV-selective if it is in any group; background =
      tested bins of those transcripts that are EV-selective in no group)."""
import sys, os, gzip, csv, collections, math
def fisher_exact(t):
    """two-sided Fisher exact test (sum of table probabilities <= the observed one, relative tolerance 1e-7), as R / scipy"""
    (a, b), (c, d) = t; r1, c1, n = a + b, a + c, a + b + c + d
    pr = lambda x: math.comb(c1, x) * math.comb(n - c1, r1 - x) / math.comb(n, r1)
    lo, hi = max(0, r1 + c1 - n), min(r1, c1); p0 = pr(a)
    p = min(1.0, sum(pr(x) for x in range(lo, hi + 1) if pr(x) <= p0 * (1 + 1e-7)))
    o = (a * d / (b * c)) if b * c > 0 else (float("inf") if a * d > 0 else float("nan"))
    return o, p
def bh(ps):
    m = len(ps); o = sorted(range(m), key=lambda i: ps[i]); q = [0.0] * m; prev = 1.0
    for rank in range(m, 0, -1):
        i = o[rank - 1]; prev = min(prev, ps[i] * m / rank); q[i] = prev
    return q
assert abs(fisher_exact([[3, 1], [1, 3]])[1] - 0.4857142857) < 1e-9 and abs(fisher_exact([[8, 2], [1, 5]])[1] - 0.0349650350) < 1e-9
assert [round(x, 6) for x in bh([0.01, 0.04, 0.03, 0.2])] == [0.04, 0.053333, 0.053333, 0.2]
K3, E1, ANN, FA, OUT = sys.argv[1:6]; os.makedirs(OUT, exist_ok=False)
GROUPS = ["B8441_untreated", "B8441_caspofungin", "MMC1_untreated", "MMC1_caspofungin"]; TYPES = ["PCG_CDS", "lncRNA", "tRNA", "5S_rRNA"]
TESTS = ["primary", "draw2", "draw3", "gcclass", "naive"]
fl = lambda x: float("nan") if x in ("NA", "") else float(x)
sig = lambda r, k, d: r[f"q_{k}"] not in ("NA", "") and float(r[f"q_{k}"]) < 0.05 and r[f"direction_{k}"] == d
cand = []; persig = []
for T in TYPES:
    W = {G: {r["word"]: r for r in csv.DictReader(gzip.open(f"{K3}/{G}.{T}.words.tsv.gz", "rt"), delimiter="\t")} for G in GROUPS}
    for G in GROUPS:
        for k in TESTS:
            persig.append((T, G, k, sum(sig(r, k, "higher_in_EV") for r in W[G].values()), sum(sig(r, k, "higher_in_Cell") for r in W[G].values())))
    for w in W[GROUPS[0]]:
        for d in ("higher_in_EV", "higher_in_Cell"):
            n = {k: sum(sig(W[G][w], k, d) for G in GROUPS) for k in TESTS}
            if n["primary"] >= 3:
                row = dict(type=T, word=w, gc_class=W[GROUPS[0]][w]["gc_class"], direction=d,
                           **{f"groups_sig_same_dir_{k}": n[k] for k in TESTS})
                for G in GROUPS:
                    r = W[G][w]; e, c = fl(r["share_EV_primary"]), fl(r["share_Cell_primary"])
                    row[f"q_primary_{G}"] = r["q_primary"]; row[f"FE_EV_over_Cell_{G}"] = (e / c) if c and c == c and c > 0 else "NA"
                    row[f"q_gcclass_{G}"] = r["q_gcclass"]
                cand.append(row)
with open(f"{OUT}/significant_words_per_group.tsv", "w") as w:
    w.write("type\tgroup\ttest\tq05_higher_in_EV\tq05_higher_in_Cell\n"); [w.write("\t".join(map(str, x)) + "\n") for x in persig]
ck = list(cand[0]) if cand else ["type", "word", "direction"]
with open(f"{OUT}/tag_candidates.tsv", "w") as w:
    w.write("\t".join(ck) + "\n"); [w.write("\t".join(str(r.get(k, "")) for k in ck) + "\n") for r in cand]
print("candidates per type/direction:", collections.Counter((r["type"], r["direction"]) for r in cand))
# ---------------- Q6
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
st = {}   # group -> {bin_id: (transcript, rna_class, tested, ev_selective)}
for G in GROUPS:
    d = {}
    for r in csv.DictReader(gzip.open(f"{E1}/{G}.ev_vs_cell_bins.tsv.gz", "rt"), delimiter="\t"):
        d[r["bin_id"]] = (r["transcript_id"], r["rna_class"], r["status"] == "tested", r["category_q0.05"].endswith("_EV_share_higher"))
    st[G] = d
def sets(groups, cls):
    fg = {b for G in groups for b, v in st[G].items() if v[3] and v[1] == cls}
    tx = {st[G][b][0] for G in groups for b in fg if b in st[G]}
    bg = {b for G in groups for b, v in st[G].items() if v[0] in tx and v[2] and v[1] == cls} - {b for G in GROUPS for b, v in st[G].items() if v[3]}
    return sorted(fg), sorted(bg), len(tx)
q6 = []
for T, cls in (("PCG_CDS", "PCG_CDS"), ("lncRNA", "lncRNA")):
    words = [r for r in cand if r["type"] == T]
    for label, groups in [("pooled", GROUPS)] + [(G, [G]) for G in GROUPS]:
        fg, bg, ntx = sets(groups, cls)
        if not words or not fg or not bg: q6.append(dict(type=T, set=label, word="", note=f"not tested: candidates {len(words)}, EV-selective bins {len(fg)}, background bins {len(bg)}")); continue
        fs = [bseq(b) for b in fg]; bs = [bseq(b) for b in bg]; rows = []
        for r in words:
            a = sum(r["word"] in s for s in fs); c = sum(r["word"] in s for s in bs)
            o, p = fisher_exact([[a, len(fs) - a], [c, len(bs) - c]])
            rows.append(dict(type=T, set=label, word=r["word"], word_direction=r["direction"], transcripts=ntx, EV_selective_bins=len(fs),
                             EV_selective_with_word=a, background_bins=len(bs), background_with_word=c,
                             frac_EV_selective=round(a / len(fs), 4), frac_background=round(c / len(bs), 4), odds_ratio=o, p=p))
        qs = bh([x["p"] for x in rows])
        for x, q in zip(rows, qs): x["q_BH_over_candidates"] = q; x["note"] = ""
        q6 += rows
qk = ["type", "set", "word", "word_direction", "transcripts", "EV_selective_bins", "EV_selective_with_word", "background_bins",
      "background_with_word", "frac_EV_selective", "frac_background", "odds_ratio", "p", "q_BH_over_candidates", "note"]
with open(f"{OUT}/q6_fisher_ev_selective_bins.tsv", "w") as w:
    w.write("\t".join(qk) + "\n"); [w.write("\t".join(str(r.get(k, "")) for k in qk) + "\n") for r in q6]
print("Q6 rows:", len(q6), "q<0.05:", sum(1 for r in q6 if r.get("q_BH_over_candidates", 1) < 0.05))
print("K4_DONE")
