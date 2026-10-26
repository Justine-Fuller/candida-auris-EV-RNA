#!/usr/bin/env python3
"""test_fl_extra_count.py -- known-answer test for fl_extra_count.py (synthetic BAM; 2026-10-09).
Usage: test_fl_extra_count.py <scripts_dir> <workdir>"""
import sys, os, json, subprocess, pysam
S, W = sys.argv[1:3]; os.makedirs(W, exist_ok=True)
H = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": "chrA", "LN": 100000}]}
hdr100 = "bin_id\trna_class\tmodel_id\tgene_id\tchrom\tstrand\trole\ttranscript_length\tn_bins\tbin_number\ttx_start0\ttx_end\tbin_length\tgenomic_blocks_5to3\n"
def ann(path, rows):
    with open(path, "w") as f: f.write(hdr100); [f.write("\t".join(map(str, r)) + "\n") for r in rows]
# existing model E (+, 1000-1300, 3 bins of 100); new model N (-, 5000-5250: bins 5'->3' = 5150-5250, 5050-5150, 5000-5050)
ann(f"{W}/a100.tsv", [["E|bin001", "PCG_CDS", "E", "E", "chrA", "+", "ordinary", 300, 3, 1, 0, 100, 100, "1000-1100"],
                      ["E|bin002", "PCG_CDS", "E", "E", "chrA", "+", "ordinary", 300, 3, 2, 100, 200, 100, "1100-1200"],
                      ["E|bin003", "PCG_CDS", "E", "E", "chrA", "+", "ordinary", 300, 3, 3, 200, 300, 100, "1200-1300"]])
ann(f"{W}/n100.tsv", [["N|bin001", "PCG_CDS", "N", "N", "chrA", "-", "ordinary", 250, 3, 1, 0, 100, 100, "5150-5250"],
                      ["N|bin002", "PCG_CDS", "N", "N", "chrA", "-", "ordinary", 250, 3, 2, 100, 200, 100, "5050-5150"],
                      ["N|bin003", "PCG_CDS", "N", "N", "chrA", "-", "ordinary", 250, 3, 3, 200, 250, 50, "5000-5050"]])
ann(f"{W}/a150.tsv", [["E|bin001", "PCG_CDS", "E", "E", "chrA", "+", "ordinary", 300, 2, 1, 0, 150, 150, "1000-1150"],
                      ["E|bin002", "PCG_CDS", "E", "E", "chrA", "+", "ordinary", 300, 2, 2, 150, 300, 150, "1150-1300"]])
ann(f"{W}/n150.tsv", [["N|bin001", "PCG_CDS", "N", "N", "chrA", "-", "ordinary", 250, 2, 1, 0, 150, 150, "5100-5250"],
                      ["N|bin002", "PCG_CDS", "N", "N", "chrA", "-", "ordinary", 250, 2, 2, 150, 250, 100, "5000-5100"]])
with open(f"{W}/5S.bed", "w") as f:
    f.write("#chrom\tstart0\tend\tname\tlength\tstrand\tv2_id\n")
    f.write("chrA\t20000\t20118\tS1\t118\t+\tx\nchrA\t30000\t30118\tS2\t118\t-\ty\n")
def rec(name, pos, cig, rev=False, nh=1, sec=False, seq_len=None):
    a = pysam.AlignedSegment(); a.query_name = name; a.reference_id = 0; a.reference_start = pos; a.cigarstring = cig
    n = sum(l for o, l in a.cigartuples if o in (0, 1, 4)); a.query_sequence = "A" * n; a.flag = (16 if rev else 0) | (256 if sec else 0)
    a.mapping_quality = 255; a.set_tag("NH", nh); return a
R = [rec("r1", 1010, "30M"),                    # E bin1, physical
     rec("r2", 1085, "30M"),                    # E bin1 15 / bin2 15 -> neighbouring tie -> physical (sRNA)
     rec("r3", 1010, "15M"),                    # < 18 nt -> ignored
     rec("r4", 5160, "30M", rev=True),          # N (- strand) bin1 -> physical for N
     rec("r5", 5160, "30M"),                    # wrong strand for N -> no_model
     rec("m1", 20010, "40M", nh=2), rec("m1", 30050, "40M", rev=True, nh=2, sec=True),   # 5S multimapper: S1 primary, offsets 10..49
     rec("m2", 30000, "30M", rev=True),         # 5S S2 (- strand): genomic 30000-30030 = offsets 88..117 -> 12 bases bin1, 18 bin2
     rec("m3", 20090, "20M", rev=True)]         # wrong strand for S1 -> not 5S
with pysam.AlignmentFile(f"{W}/t.bam", "wb", header=H) as o: [o.write(a) for a in R]
for kind in ("sRNA", "polyA"):
    subprocess.run([sys.executable, f"{S}/fl_extra_count.py", f"{W}/t.bam", f"T_{kind}", kind, f"{W}/a100.tsv", f"{W}/n100.tsv",
                    f"{W}/a150.tsv", f"{W}/n150.tsv", f"{W}/5S.bed", W], check=True)
x = json.load(open(f"{W}/T_sRNA.flextra.json")); p = json.load(open(f"{W}/T_polyA.flextra.json"))
chk = [("sRNA E physical 2", x["physical_reads"].get("E") == 2), ("sRNA N physical 1", x["physical_reads"].get("N") == 1),
       ("sRNA E bin1 bases 30+15", x["bin_aligned_bases"].get("E|bin001") == 45), ("sRNA E bin2 bases 15", x["bin_aligned_bases"].get("E|bin002") == 15),
       ("sRNA N bin1 bases 30", x["bin_aligned_bases"].get("N|bin001") == 30),
       ("sRNA 5S reads 2 (m1, m2)", x["rRNA_5S"]["reads"] == 2), ("sRNA 5S NH1 reads 1", x["rRNA_5S"]["reads_NH1"] == 1),
       ("sRNA 5S bin1 40+12", x["rRNA_5S"]["bin_aligned_bases"]["5S_rRNA|bin001"] == 52),
       ("sRNA 5S bin2 18", x["rRNA_5S"]["bin_aligned_bases"]["5S_rRNA|bin002"] == 18),
       ("polyA flip: r5 -> N", p["physical_reads"].get("N") == 1), ("polyA flip: E none", p["physical_reads"].get("E", 0) == 0),
       ("polyA flip: m3 -> 5S S1 (offsets 90..109)", p["rRNA_5S"]["reads"] == 1 and p["rRNA_5S"]["bin_aligned_bases"] == {"5S_rRNA|bin001": 10, "5S_rRNA|bin002": 10})]
for n, v in chk: print(("PASS " if v else "FAIL ") + n)
print("TEST_PASS" if all(v for _, v in chk) else "TEST_FAIL"); sys.exit(0 if all(v for _, v in chk) else 1)
