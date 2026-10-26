#!/usr/bin/env python3
"""Known-answer test for count_bins100.py (C0) and count_bins100_c1.py (C1), Option-1 100-nt branch (2026-10-08).
Usage: test_counters100.py <scriptdir> <workdir>. r4 (2026-10-09): transcript_totals has the r3 excluded_3prime_overhang_reads column. r3: also checks kind "sRNA" (Cell sRNA) = EV outputs and polyA refused. r2: r1 gave the two counters different library names, so the
library column differed (all counts identical); r2 runs both as "ev" in separate subdirectories. Both counters must give exactly the expected answers AND identical files.
Toy models (bins of 100 nt in spliced transcript coordinates, genomic 0-based half-open blocks listed 5'->3'):
  A '+' exons 100-180, 400-470 (L=150): bin1 = 100-180 + 400-420 (crosses the splice junction), bin2 = 420-470 (50 nt, shorter)
  B '+' exon 1000-1250 (L=250): bin1 1000-1100, bin2 1100-1200, bin3 1200-1250 (50 nt)
  C '-' exon 2000-2200 (L=200): bin1 = 2100-2200 (5' end on minus strand), bin2 = 2000-2100
  D '+' competitor_only 3000-3150; E '+' 1220-1400 (same strand, overlaps B); G '+' 9950-10150 (bin1 crosses 10-kb bucket)
  H '+' 5000-5300 (3 full bins)"""
import sys, os, subprocess, json, gzip, pysam
S, W = sys.argv[1], sys.argv[2]; os.makedirs(W, exist_ok=True)
models = {
 "A": ("+", "ordinary", [[(100, 180), (400, 420)], [(420, 470)]]),
 "B": ("+", "ordinary", [[(1000, 1100)], [(1100, 1200)], [(1200, 1250)]]),
 "C": ("-", "ordinary", [[(2100, 2200)], [(2000, 2100)]]),
 "D": ("+", "competitor_only", [[(3000, 3100)], [(3100, 3150)]]),
 "E": ("+", "ordinary", [[(1220, 1320)], [(1320, 1400)]]),
 "G": ("+", "ordinary", [[(9950, 10050)], [(10050, 10150)]]),
 "H": ("+", "ordinary", [[(5000, 5100)], [(5100, 5200)], [(5200, 5300)]]),
}
ann = f"{W}/bins.tsv"
with open(ann, "w") as o:
    o.write("bin_id\trna_class\tmodel_id\tgene_id\tchrom\tstrand\trole\tgenomic_blocks_5to3\n")
    for m, (st, role, bl) in models.items():
        for i, b in enumerate(bl):
            o.write(f"{m}|bin{i+1:03d}\tPCG_CDS\t{m}\t{m}\tc1\t{st}\t{role}\t{';'.join(f'{s}-{e}' for s, e in b)}\n")
hdr = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "c1", "LN": 20000}]}
def bam(path, reads):
    with pysam.AlignmentFile(path, "wb", header=hdr) as f:
        for q, pos, cig, rev, nh, flag in reads:
            a = pysam.AlignedSegment(); a.query_name = q
            if cig is None:
                a.flag = 4; a.query_sequence = "A" * 30; a.reference_id = -1; a.reference_start = -1
            else:
                a.cigarstring = cig; a.query_sequence = "A" * a.infer_query_length()
                a.flag = (16 if rev else 0) | flag; a.reference_id = 0; a.reference_start = pos; a.set_tag("NH", nh)
            f.write(a)
ev = [
 ("r01_user_example_20v50", 1080, "70M", 0, 1, 0),   # B bin1 20, bin2 50 -> B|bin002
 ("r02_50v20", 1050, "70M", 0, 1, 0),                # B bin1 50, bin2 20 -> B|bin001
 ("r03_tie_25v25", 1075, "50M", 0, 1, 0),            # tie, unresolved
 ("r04_tie_via_deletion", 1090, "10M5D10M", 0, 1, 0),  # 10 / 10 (deleted 1100-1105 contributes 0) -> tie
 ("r05_softclip", 1095, "10S20M", 0, 1, 0),          # aligned 1095-1115: 5 / 15 -> B|bin002 (read length 30)
 ("r06_spliced_A", 170, "10M220N30M", 0, 1, 0),      # 170-180 + 400-430: bin1 10+20=30, bin2 10 -> A|bin001
 ("r07_minus_C", 2090, "30M", 1, 1, 0),              # reverse read on '-' model: 2090-2120: bin2 10, bin1 20 -> C|bin001
 ("r08_wrong_strand", 2150, "20M", 0, 1, 0),         # '+' read on '-' model -> no_binned_model
 ("r09_NH2", 1000, "30M", 0, 2, 0),                  # not eligible
 ("r10_17nt", 1000, "17M", 0, 1, 0),                 # not eligible
 ("r11_two_models", 1210, "30M", 0, 1, 0),           # B bin3 30 + E bin1 20 -> ambiguous_model
 ("r12_competitor", 3010, "20M", 0, 1, 0),           # competitor_only
 ("r13_unmapped", 0, None, 0, 1, 0),
 ("r09_NH2", 1200, "30M", 0, 2, 256),                # secondary record, skipped
 ("r15_bucket", 9999, "151M", 0, 1, 0),              # G bin1 51, bin2 100 -> G|bin002 (r1 bucket bug would double bin1 -> 102)
 ("r16_three_bins", 1090, "120M", 0, 1, 0),          # B 10 / 100 / 10 -> B|bin002
 ("r17_tie_3bins", 5000, "300M", 0, 1, 0),           # H 100 / 100 / 100 -> 3-bin tie, unresolved
 ("r18_single_bin", 1010, "20M", 0, 1, 0),           # entirely inside B bin1 -> B|bin001
]
bam(f"{W}/ev.bam", ev)
exp_acc = {"records": 18, "primary_mapped": 16, "eligible_NH1_ge18": 14, "assigned": 8, "no_binned_model": 1,
           "ambiguous_model": 1, "competitor_only": 1, "tie_unresolved_2bin": 2, "tie_unresolved_3plus": 1,
           "assigned_touching_2plus_bins": 7}
exp_bins = {"B|bin001": 2, "B|bin002": 3, "A|bin001": 1, "C|bin001": 1, "G|bin002": 1}
exp_tx = {"A": (1, 1, 0, 0), "B": (5, 5, 2, 0), "C": (1, 1, 0, 0), "G": (1, 1, 0, 0), "H": (0, 0, 1, 0)}   # r4: + excluded_3prime_overhang_reads column (0: no read hangs off a 3' end here)
exp_ties = {"r03_tie_25v25": ("B", "B|bin001:25;B|bin002:25"), "r04_tie_via_deletion": ("B", "B|bin001:10;B|bin002:10"),
            "r17_tie_3bins": ("H", "H|bin001:100;H|bin002:100;H|bin003:100")}
exp_ex = {"r01_user_example_20v50": ("B|bin002", "B|bin002:50;B|bin001:20"), "r02_50v20": ("B|bin001", "B|bin001:50;B|bin002:20"),
          "r05_softclip": ("B|bin002", "B|bin002:15;B|bin001:5"), "r06_spliced_A": ("A|bin001", "A|bin001:30;A|bin002:10"),
          "r07_minus_C": ("C|bin001", "C|bin001:20;C|bin002:10"), "r15_bucket": ("G|bin002", "G|bin002:100;G|bin001:51"),
          "r16_three_bins": ("B|bin002", "B|bin002:100;B|bin001:10;B|bin003:10")}
rd = lambda p: gzip.open(p, "rt").read()
tab = lambda p: [l.split("\t") for l in rd(p).splitlines()[1:]]
ok = True; files = {}
for tag, script, sfx in (("C0", "count_bins100.py", ""), ("C1", "count_bins100_c1.py", ".c1")):
    D = f"{W}/{tag}"; os.makedirs(D, exist_ok=True)      # r2: same library name for both counters, separate dirs
    p = subprocess.run([sys.executable, f"{S}/{script}", ann, f"{W}/ev.bam", "ev", "EV", D], capture_output=True, text=True)
    if p.returncode: print(tag, "FAILED TO RUN", p.stderr); ok = False; continue
    acc = json.load(open(f"{D}/ev.accounting{sfx}.json"))
    bc = {r[0]: int(r[1]) for r in tab(f"{D}/ev.bin_counts{sfx}.tsv.gz")}
    tx = {r[0]: tuple(map(int, r[1:])) for r in tab(f"{D}/ev.transcript_totals{sfx}.tsv.gz")}
    ti = {r[1]: (r[7], r[11]) for r in tab(f"{D}/ev.ties{sfx}.tsv.gz")}
    ex = {r[1]: (r[12], r[11]) for r in tab(f"{D}/ev.multibin_examples{sfx}.tsv.gz")}
    res = dict(accounting={k: acc.get(k, 0) for k in exp_acc} == exp_acc and acc.get("closes", True),
               bins=bc == exp_bins, transcripts=tx == exp_tx, ties=ti == exp_ties, examples=ex == exp_ex)
    print(tag, res)
    if not all(res.values()): ok = False; print("  got", acc, bc, tx, ti, ex, sep="\n  ")
    files[tag] = [rd(f"{D}/ev.{n}{sfx}.tsv.gz") for n in ("bin_counts", "transcript_totals", "ties", "multibin_examples", "overlap_patterns")]
same = files.get("C0") == files.get("C1") and "C0" in files
print("C0_C1_files_identical", same); ok &= same
# r2: the same reads labelled as a Cell sRNA library ("sRNA") must give byte-identical outputs; polyA must be refused
for tag, script, sfx in (("C0", "count_bins100.py", ""), ("C1", "count_bins100_c1.py", ".c1")):
    D = f"{W}/{tag}_sRNA"; os.makedirs(D, exist_ok=True)
    p = subprocess.run([sys.executable, f"{S}/{script}", ann, f"{W}/ev.bam", "ev", "sRNA", D], capture_output=True, text=True)
    srna = p.returncode == 0 and [rd(f"{D}/ev.{n}{sfx}.tsv.gz") for n in ("bin_counts", "transcript_totals", "ties", "multibin_examples", "overlap_patterns")] == files.get(tag)
    os.makedirs(f"{W}/{tag}_polyA", exist_ok=True)
    q = subprocess.run([sys.executable, f"{S}/{script}", ann, f"{W}/ev.bam", "ev", "polyA", f"{W}/{tag}_polyA"], capture_output=True, text=True)
    refused = q.returncode != 0 and "polyA refused" in q.stderr and not os.listdir(f"{W}/{tag}_polyA")
    print(tag, "Cell_sRNA_identical_to_EV", srna, "polyA_refused", refused); ok &= srna and refused
print("ALL_COUNTER100_TESTS_PASS" if ok else "COUNTER100_TESTS_FAILED"); sys.exit(0 if ok else 1)
