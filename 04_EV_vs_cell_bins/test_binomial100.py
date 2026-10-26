#!/usr/bin/env python3
"""Known-answer test for binomial_test100.R (2026-10-08). Usage: test_binomial100.py <scriptdir> <workdir> <Rscript>
Toy: 3 libraries (one group). T1: 1,000 nt, 10 bins of 100; T2: 1,000 nt, 10 bins; R1: competitor-only rDNA transcript;
O1: lncRNA with bin 2 flagged as overlapping another same-strand transcript.
Expected (exact, computed here with fractions, independent of R):
  T1 bin1, lib A: n=10, k=5, p0=0.1 -> p = sum_{i=5..10} C(10,i) 0.1^i 0.9^(10-i)
  T1 in lib C has n=9 -> every T1 bin 'not_tested_lt10_reads' in C; T2 has 0 reads in C -> 'not_tested_no_reads'
  R1 bins -> 'excluded_rDNA_competitor_only'; O1 bin2 -> 'excluded_same_strand_overlap' (O1 bin1/bin3 tested)
  every tested p equals the exact tail; dataset-wide q = BH over tested bins only (recomputed here)
  group call at 0.10: T1 bin1 called in A and B -> n_replicates_called 2 -> group_call_2of3 TRUE"""
import sys, os, subprocess, gzip, csv
from fractions import Fraction
from math import comb
S, W, RS = sys.argv[1:4]; os.makedirs(W, exist_ok=True)
def tail(k, n, p):
    p = Fraction(p).limit_denominator(10 ** 9)
    return float(sum(comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k, n + 1)))
T = {"T1": ("PCG_CDS", "ordinary", 10, [100] * 10, ""), "T2": ("PCG_CDS", "ordinary", 10, [100] * 10, ""),
     "R1": ("lncRNA", "competitor_only", 3, [100, 100, 150], ""), "O1": ("lncRNA", "ordinary", 3, [100, 100, 120], "bin002")}
ann_cols = ["bin_id", "rna_class", "model_id", "gene_id", "role", "transcript_length", "n_bins", "bin_number", "bin_length",
            "same_strand_binned_model_overlap", "merged_final_bin"]
with open(f"{W}/ann.tsv", "w") as o:
    o.write("\t".join(ann_cols) + "\n")
    for m, (cls, role, nb, lens, ov) in T.items():
        for i, L in enumerate(lens, 1):
            o.write("\t".join([f"{m}|bin{i:03d}", cls, m, m, role, str(sum(lens)), str(nb), str(i), str(L),
                               "OTHER" if ov == f"bin{i:03d}" else "", "yes" if L > 100 else ""]) + "\n")
# counts per library: (transcript physical reads, {bin: count})
C = {"A": {"T1": (10, {1: 5, 2: 2, 3: 3}), "T2": (40, {4: 30, 5: 10}), "R1": (0, {}), "O1": (12, {1: 2, 2: 4, 3: 6})},
     "B": {"T1": (12, {1: 6, 2: 6}), "T2": (20, {4: 20}), "R1": (0, {}), "O1": (10, {3: 10})},
     "C": {"T1": (9, {1: 9}), "T2": (0, {}), "R1": (0, {}), "O1": (11, {1: 11})}}
cols = ["library", "strain", "treatment", "replicate", "transcript_id", "gene_id", "rna_class", "role", "rdna_flag", "bin_id",
        "bin_number", "n_bins", "bin_tx_start0", "bin_tx_end", "bin_length", "final_bin_shorter", "bin_reads_strict_greatest_overlap",
        "bin_tie_plus1", "bin_count", "transcript_physical_reads", "transcript_strict_assigned_reads", "transcript_tie_reads",
        "transcript_summed_bin_counts", "transcript_unassigned_nonadjacent_tie_reads"]
os.makedirs(f"{W}/tables", exist_ok=True)
with open(f"{W}/man.tsv", "w") as mo:
    mo.write("library\tbam\tkind\tcompartment\tgroup\treplicate\n")
    for r, L in enumerate("ABC", 1):
        mo.write(f"{L}\tx\tEV\tEV\tB8441_untreated\t{r}\n")
        with gzip.open(f"{W}/tables/{L}.bin100_final_table.tsv.gz", "wt") as o:
            o.write("\t".join(cols) + "\n")
            for m, (cls, role, nb, lens, ov) in T.items():
                n, bc = C[L][m]; s = 0
                for i, bl in enumerate(lens, 1):
                    o.write("\t".join(map(str, [L, "B8441", "untreated", r, m, m, cls, role, "", f"{m}|bin{i:03d}", i, nb, s, s + bl, bl, "",
                                                bc.get(i, 0), 0, bc.get(i, 0), n, n, 0, n, 0])) + "\n"); s += bl
p = subprocess.run([RS, f"{S}/binomial_test100.R", f"{W}/tables", f"{W}/ann.tsv", f"{W}/man.tsv", f"{W}/out", "20261008"], capture_output=True, text=True)
ok = p.returncode == 0
if not ok: print("R FAILED", p.stdout[-1500:], p.stderr[-1500:])
else:
    R = {L: {x["bin_id"]: x for x in csv.DictReader(gzip.open(f"{W}/out/{L}.binomial_results.tsv.gz", "rt"), delimiter="\t")} for L in "ABC"}
    exp_p = tail(5, 10, 0.1); got = float(R["A"]["T1|bin001"]["p"])
    c1 = abs(got - exp_p) < 1e-12 * max(1, exp_p); print("T1 bin1 lib A p:", got, "exact:", exp_p, c1)
    allp = True
    for L in "ABC":
        for b, x in R[L].items():
            if x["status"] == "tested":
                allp &= abs(float(x["p"]) - tail(int(x["k"]), int(x["n"]), Fraction(int(x["bin_length"]), int(x["transcript_length"])))) < 1e-12
    print("all tested p = exact tail:", allp)
    st = lambda L, b: R[L][b]["status"]
    c3 = (all(st("C", f"T1|bin{i:03d}") == "not_tested_lt10_reads" for i in range(1, 11)) and st("C", "T2|bin001") == "not_tested_no_reads"
          and all(st(L, f"R1|bin{i:03d}") == "excluded_rDNA_competitor_only" for L in "ABC" for i in (1, 2, 3))
          and all(st(L, "O1|bin002") == "excluded_same_strand_overlap" for L in "ABC") and st("A", "O1|bin001") == "tested")
    print("statuses:", c3)
    def bh(ps):
        o = sorted(range(len(ps)), key=lambda i: ps[i]); m = len(ps); q = [0] * m; run = 1.0
        for r in range(m - 1, -1, -1): run = min(run, ps[o[r]] * m / (r + 1)); q[o[r]] = run
        return q
    c4 = True
    for L in "ABC":
        te = [x for x in R[L].values() if x["status"] == "tested"]
        q = bh([float(x["p"]) for x in te]); c4 &= all(abs(q[i] - float(x["q_dataset"])) < 1e-12 for i, x in enumerate(te))
    print("dataset BH recomputed:", c4)
    g = {x["bin_id"]: x for x in csv.DictReader(gzip.open(f"{W}/out/B8441_untreated.group_calls.tsv.gz", "rt"), delimiter="\t")}
    qa, qb = float(R["A"]["T1|bin001"]["q_dataset"]), float(R["B"]["T1|bin001"]["q_dataset"])
    want = sum(q < 0.10 for q in (qa, qb))
    c5 = int(g["T1|bin001"]["n_replicates_called_q0.1"]) == want and (g["T1|bin001"]["group_call_2of3_q0.1"] == "TRUE") == (want >= 2)
    print("group call T1 bin1 (q_A, q_B):", qa, qb, c5)
    ok = c1 and allp and c3 and c4 and c5
print("ALL_BINOMIAL_TESTS_PASS" if ok else "BINOMIAL_TESTS_FAILED"); sys.exit(0 if ok else 1)
