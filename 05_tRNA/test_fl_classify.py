#!/usr/bin/env python3
"""Known-answer test for fl_classify.py. Usage: test_fl_classify.py <scriptdir> <workdir>
EV: X full-length (coverage ok in 2/3), Y fragment (one bin at 1x), tRNA family Gly/GCC full-length (CCA 1). Cell: X full-length,
Y not detected (0 reads), Gly/GCC fragment-only (coverage/TPM ok but 0 CCA). 'With rRNA': 1e8 rDNA reads per library push every
TPM below 100 -> nothing full-length."""
import sys, os, json, gzip, subprocess, csv
S, W = sys.argv[1], sys.argv[2]
for d in ("t1", "cen", "ev", "cell"): os.makedirs(f"{W}/{d}", exist_ok=True)
with open(f"{W}/ann.tsv", "w") as o:
    o.write("bin_id\tmodel_id\trna_class\trole\tbin_length\ttranscript_length\n")
    for m in ("X", "Y"): o.write(f"{m}|bin001\t{m}\tPCG_CDS\tordinary\t100\t200\n{m}|bin002\t{m}\tPCG_CDS\tordinary\t100\t200\n")
    o.write("LNCRNA.1946|bin001\tLNCRNA.1946\tlncRNA\tcompetitor_only\t100\t100\n")
open(f"{W}/t.bed", "w").write("#chrom\tstart0\tend\tname\tscore\tstrand\tctg_family\tfungidb_gene_id\tintron\ttier\nc1\t0\t72\tT|Gly/GCC\t70\t+\tGly/GCC\tNA\tNA\tA_core\n")
libs = [("EV1", "EV"), ("EV2", "EV"), ("EV3", "EV"), ("C1", "Cell_sRNA"), ("C2", "Cell_sRNA"), ("C3", "Cell_sRNA")]
open(f"{W}/man.tsv", "w").write("library\tbam\tkind\tcompartment\tgroup\treplicate\n" + "".join(f"{L}\tx\tsRNA\t{c}\tG\t1\n" for L, c in libs))
for L, c in libs:
    ev = c == "EV"
    bb = {"X|bin001": 400, "X|bin002": (200 if L == "EV3" else 400)}
    if ev: bb.update({"Y|bin001": 400, "Y|bin002": 100})
    json.dump({"bin_aligned_bases": bb, "trna_families": {"Gly/GCC": {"reads": 40, "bases": 360, "CCA_reads": 1 if L == "EV1" else 0}}},
              open(f"{W}/t1/{L}.trna.json", "w"))
    json.dump({"primary_class": {"rDNA_lncRNA|NHgt1": 100000000}}, open(f"{W}/cen/{L}.census.json", "w"))
    with gzip.open(f"{W}/{'ev' if ev else 'cell'}/{L}.bin100_final_table.tsv.gz", "wt") as o:
        o.write("transcript_id\ttranscript_physical_reads\n")
        o.write(f"X\t50\nX\t50\nY\t{50 if ev else 0}\nY\t{50 if ev else 0}\n")
p = subprocess.run([sys.executable, f"{S}/fl_classify.py", f"{W}/t1", f"{W}/man.tsv", f"{W}/ann.tsv", f"{W}/t.bed", f"{W}/cen", f"{W}/ev", f"{W}/cell", f"{W}/out"],
                   capture_output=True, text=True)
if p.returncode: print(p.stdout, p.stderr); print("FL_CLASSIFY_TEST_FAILED"); sys.exit(1)
U = {(r["compartment"], r["unit"]): r for r in csv.DictReader(gzip.open(f"{W}/out/unit_table.tsv.gz", "rt"), delimiter="\t")}
f0 = lambda c, u: U[(c, u)]["full_length_groups_TPM_without_rRNA"] == "G"; f1 = lambda c, u: U[(c, u)]["full_length_groups_TPM_with_rRNA"] == "G"
dt = lambda c, u: U[(c, u)]["detected_groups"] == "G"
chk = {"EV_X_full_length": f0("EV", "X") and dt("EV", "X"), "EV_Y_fragment": dt("EV", "Y") and not f0("EV", "Y"),
       "EV_tRNA_full_length_with_CCA": f0("EV", "Gly/GCC"), "Cell_X_full_length": f0("Cell_sRNA", "X"),
       "Cell_Y_not_detected": not dt("Cell_sRNA", "Y"), "Cell_tRNA_fragment_no_CCA": dt("Cell_sRNA", "Gly/GCC") and not f0("Cell_sRNA", "Gly/GCC"),
       "with_rRNA_nothing_full_length": not any(f1(c, u) for c in ("EV", "Cell_sRNA") for u in ("X", "Y", "Gly/GCC")),
       "rDNA_lncRNA_not_a_unit": ("EV", "LNCRNA.1946") not in U}
for k, v in chk.items(): print(k, v)
print("ALL_FL_CLASSIFY_TESTS_PASS" if all(chk.values()) else "FL_CLASSIFY_TEST_FAILED"); sys.exit(0 if all(chk.values()) else 1)
