#!/usr/bin/env python3
"""test_g1_6g.py -- known-answer test of g1_count6g.py (units k4 and k5drach) on a tiny synthetic group.
Six identical libraries (3 EV + 3 Cell_sRNA) so GC matching keeps every read; expected counts computed by hand.
Usage: test_g1_6g.py <scripts dir> <workdir>"""
import sys, os, gzip, json, csv, subprocess
SD, WD = sys.argv[1:3]; os.makedirs(f"{WD}/K1", exist_ok=False)
reads = {"GGACTGGACA": 2,      # 4-mers: GGAC x2 per read; DRACH windows: GGACT, GGACA -> 2 per read
         "GGACNGGACA": 1,      # N windows skipped: GGAC x2, GACA x1; DRACH: GGACA x1 (GGACN has N)
         "TAACCGGACG": 3,      # DRACH: TAACC yes; GGACG no (H != G) -> 1 per read; GGAC x1 per read
         "CGACTAAAAAA": 1}     # CGACT: D != C -> no DRACH; six A's -> AAAA x3 (overlapping); GGAC 0
man = ["library\tcompartment\tgroup"]
for comp in ("EV", "Cell_sRNA"):
    for i in (1, 2, 3):
        L = f"T_{comp}_{i}"; man.append(f"{L}\t{comp}\tTEST")
        with gzip.open(f"{WD}/K1/{L}.k1.json.gz", "wt") as fh:
            json.dump(dict(library=L, accounting={}, types={"PCG_CDS": reads, "lncRNA": {}, "tRNA": {}, "5S_rRNA": {}}), fh)
open(f"{WD}/manifest.tsv", "w").write("\n".join(man) + "\n")
exp = {"k4": {"GGAC": 2 * 2 + 2 * 1 + 1 * 3 + 0, "GACA": 2 + 1, "AAAA": 3, "ACNG": None},
       "k5drach": {"DRACH": 2 * 2 + 1 * 1 + 1 * 3 + 0, "GGACG": 3, "CGACT": 1}}
ok = True
for unit, e in exp.items():
    out = f"{WD}/K2_{unit}"
    subprocess.run([sys.executable, f"{SD}/g1_count6g.py", f"{WD}/manifest.tsv", f"{WD}/K1", out, "1", unit], check=True, capture_output=True)
    X = {r["word"]: r for r in csv.DictReader(gzip.open(f"{out}/TEST.PCG_CDS.counts.tsv.gz", "rt"), delimiter="\t")}
    nrow = len(X); want = 256 if unit == "k4" else 1007
    good = nrow == want; print("PASS" if good else "FAIL", unit, "rows", nrow, "expected", want); ok &= good
    for w, v in e.items():
        if v is None:
            good = w not in X; print("PASS" if good else "FAIL", unit, w, "absent (contains N)"); ok &= good; continue
        for col in ("T_EV_1__all", "T_Cell_sRNA_2__all", "T_EV_3__m", "T_Cell_sRNA_1__m"):
            got = int(X[w][col]); good = got == v; ok &= good
            print("PASS" if good else "FAIL", unit, w, col, got, "expected", v)
    if unit == "k5drach":
        good = all(w not in X for w in ("GGACT", "GGACA", "TAACC")) and X["DRACH"]["gc_class"] == "NA"; ok &= good
        print("PASS" if good else "FAIL", "DRACH variants pooled (not separate rows), gc_class NA")
print("ALL_G1_COUNTER_TESTS_PASS" if ok else "G1_COUNTER_TESTS_FAILED"); sys.exit(0 if ok else 1)
