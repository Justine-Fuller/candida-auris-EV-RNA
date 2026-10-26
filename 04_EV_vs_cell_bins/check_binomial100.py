#!/usr/bin/env python3
"""check_binomial100.py -- INDEPENDENT Python re-check of binomial_test100.R outputs on the real data (2026-10-08).
Usage: check_binomial100.py <final_tables_dir> <merged_annotation.tsv> <manifest.tsv> <binomial_outdir> <seed>
X1 k, n, p0 and status re-derived from the final count tables + annotation (D2, D5, D6 rules) == R output, every bin
X2 p for a seeded random sample of 1,000 tested bins per library == exact upper tail computed in log space (lgamma)
X3 dataset-wide and within-transcript BH recomputed in Python for ALL tested bins == R
Writes <binomial_outdir>/independent_check.tsv; exit 1 on any failure."""
import sys, csv, gzip, math, random
TD, ANN, MAN, OD, seed = sys.argv[1:6]
rnd = random.Random(int(seed))
ann = {r["bin_id"]: r for r in csv.DictReader(open(ANN), delimiter="\t")}
def tail(k, n, p):
    if k <= 0: return 1.0
    if p >= 1: return 1.0
    lt = [math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) + i * math.log(p) + (n - i) * math.log1p(-p) for i in range(k, n + 1)]
    m = max(lt); return min(1.0, math.exp(m) * sum(math.exp(x - m) for x in lt))
def bh(ps):
    o = sorted(range(len(ps)), key=lambda i: ps[i]); m = len(ps); q = [0.0] * m; run = 1.0
    for r in range(m - 1, -1, -1): run = min(run, ps[o[r]] * m / (r + 1)); q[o[r]] = run
    return q
res = []; ok_all = True
for L in [r["library"] for r in csv.DictReader(open(MAN), delimiter="\t")]:
    ft = {x["bin_id"]: x for x in csv.DictReader(gzip.open(f"{TD}/{L}.bin100_final_table.tsv.gz", "rt"), delimiter="\t")}
    R = list(csv.DictReader(gzip.open(f"{OD}/{L}.binomial_results.tsv.gz", "rt"), delimiter="\t"))
    x1 = True
    for x in R:
        f, a = ft[x["bin_id"]], ann[x["bin_id"]]; n = int(f["transcript_physical_reads"])
        st = ("excluded_rDNA_competitor_only" if a["role"] == "competitor_only" else "excluded_same_strand_overlap" if a["same_strand_binned_model_overlap"]
              else "not_tested_no_reads" if n == 0 else "not_tested_lt10_reads" if n < 10 else "tested")
        x1 &= (x["status"] == st and int(float(x["k"])) == int(f["bin_count"]) and int(float(x["n"])) == n
               and abs(float(x["p0"]) - int(a["bin_length"]) / int(a["transcript_length"])) < 1e-12)
    te = [x for x in R if x["status"] == "tested"]
    smp = rnd.sample(te, min(1000, len(te)))
    d2 = max(abs(float(x["p"]) - tail(int(float(x["k"])), int(float(x["n"])), float(x["p0"]))) / max(1e-300, tail(int(float(x["k"])), int(float(x["n"])), float(x["p0"]))) for x in smp)
    x2 = d2 < 1e-8
    q = bh([float(x["p"]) for x in te]); x3a = max(abs(q[i] - float(x["q_dataset"])) for i, x in enumerate(te)) < 1e-12
    byt = {}
    for i, x in enumerate(te): byt.setdefault(x["transcript_id"], []).append(i)
    x3b = True
    for t, idx in byt.items():
        qq = bh([float(te[i]["p"]) for i in idx]); x3b &= all(abs(qq[j] - float(te[i]["q_within_transcript"])) < 1e-12 for j, i in enumerate(idx))
    res.append((L, len(R), len(te), x1, x2, f"{d2:.2e}", x3a, x3b)); ok_all &= x1 and x2 and x3a and x3b
with open(f"{OD}/independent_check.tsv", "w") as o:
    o.write("library\tbins\tbins_tested\tX1_inputs_status\tX2_p_exact_tail\tX2_max_rel_diff\tX3_BH_dataset\tX3_BH_within\n")
    for r in res: o.write("\t".join(map(str, r)) + "\n")
print(open(f"{OD}/independent_check.tsv").read()); print("INDEPENDENT_CHECK_PASS" if ok_all else "INDEPENDENT_CHECK_FAILED"); sys.exit(0 if ok_all else 1)
