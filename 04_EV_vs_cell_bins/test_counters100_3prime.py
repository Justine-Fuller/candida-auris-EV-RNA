#!/usr/bin/env python3
"""Known-answer test of the r3 3'-overhang rule in count_bins100.py (C0) and count_bins100_c1.py (C1), 2026-10-09.
Rule (Caro; definitions Justine): exclude a read when HALF OR MORE of its aligned bases lie beyond its transcript's 3'
end; 5' overhang and intronic bases are ignored; soft clips never count. Usage: test_counters100_3prime.py <scriptdir> <workdir> (r2: also checks the 'any' sensitivity rule)
Toy models (0-based half-open genomic blocks, listed 5'->3'):
  A '+' exons 100-180, 400-470 (3' end 470): bin1 = 100-180 + 400-420, bin2 = 420-470
  B '+' exon 1000-1250 (3' end 1250): bin1 1000-1100, bin2 1100-1200, bin3 1200-1250
  C '-' exon 2000-2200 (3' end = genomic 2000): bin1 2100-2200, bin2 2000-2100"""
import sys, os, subprocess, json, gzip, pysam
S, W = sys.argv[1], sys.argv[2]; os.makedirs(W, exist_ok=True)
models = {"A": ("+", [[(100, 180), (400, 420)], [(420, 470)]]), "B": ("+", [[(1000, 1100)], [(1100, 1200)], [(1200, 1250)]]),
          "C": ("-", [[(2100, 2200)], [(2000, 2100)]])}
ann = f"{W}/bins.tsv"
with open(ann, "w") as o:
    o.write("bin_id\trna_class\tmodel_id\tgene_id\tchrom\tstrand\trole\tgenomic_blocks_5to3\n")
    for m, (st, bl) in models.items():
        for i, b in enumerate(bl):
            o.write(f"{m}|bin{i+1:03d}\tPCG_CDS\t{m}\t{m}\tc1\t{st}\tordinary\t{';'.join(f'{s}-{e}' for s, e in b)}\n")
hdr = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "c1", "LN": 20000}]}
reads = [  # name, pos, cigar, reverse
 ("x1_half_beyond", 1230, "40M", 0),        # B: 20 inside, 20 beyond -> exactly half -> EXCLUDED
 ("x2_most_beyond", 1240, "30M", 0),        # B: 10 inside, 20 beyond -> EXCLUDED
 ("x3_less_beyond", 1220, "50M", 0),        # B: 30 inside, 20 beyond -> kept -> B|bin003
 ("x4_5prime_overhang", 980, "30M", 0),     # B: 20 before the 5' end (ignored), 10 inside -> kept -> B|bin001
 ("x5_minus_beyond", 1980, "30M", 1),       # C: 20 beyond (< 2000), 10 inside -> EXCLUDED
 ("x6_minus_less", 1990, "30M", 1),         # C: 10 beyond, 20 inside -> kept -> C|bin002
 ("x7_intron", 170, "30M", 0),              # A: 10 exon, 20 intron (ignored) -> kept -> A|bin001
 ("x8_A_beyond", 460, "40M", 0),            # A: 10 inside, 30 beyond -> EXCLUDED
 ("x9_spliced_half", 461, "9M100N9M", 0),   # A: 9 inside (461-470), N, 9 beyond (570-579) -> half -> EXCLUDED
 ("x10_softclip", 1235, "15M20S", 0),       # B: 15 aligned inside, soft clip never counts -> kept -> B|bin003
]
with pysam.AlignmentFile(f"{W}/t.bam", "wb", header=hdr) as f:
    for q, pos, cig, rev in reads:
        a = pysam.AlignedSegment(); a.query_name = q; a.cigarstring = cig; a.query_sequence = "A" * a.infer_query_length()
        a.flag = 16 if rev else 0; a.reference_id = 0; a.reference_start = pos; a.set_tag("NH", 1); f.write(a)
exp_acc = {"eligible_NH1_ge18": 10, "excluded_3prime_overhang": 5, "assigned": 5, "no_binned_model": 0, "ambiguous_model": 0}
exp_bins = {"B|bin003": 2, "B|bin001": 1, "C|bin002": 1, "A|bin001": 1}
exp_tx = {"A": (1, 1, 0, 2), "B": (3, 3, 0, 2), "C": (1, 1, 0, 1)}
exp_x3 = {"x1_half_beyond": ("B", 40, 20), "x2_most_beyond": ("B", 30, 20), "x5_minus_beyond": ("C", 30, 20),
          "x8_A_beyond": ("A", 40, 30), "x9_spliced_half": ("A", 18, 9)}
rd = lambda p: gzip.open(p, "rt").read()
tab = lambda p: [l.split("\t") for l in rd(p).splitlines()[1:]]
ok = True; files = {}
for tag, script, sfx in (("C0", "count_bins100.py", ""), ("C1", "count_bins100_c1.py", ".c1")):
    D = f"{W}/{tag}"; os.makedirs(D, exist_ok=True)
    p = subprocess.run([sys.executable, f"{S}/{script}", ann, f"{W}/t.bam", "t", "EV", D], capture_output=True, text=True)
    if p.returncode: print(tag, "FAILED TO RUN", p.stderr[-1500:]); ok = False; continue
    acc = json.load(open(f"{D}/t.accounting{sfx}.json"))
    bc = {r[0]: int(r[1]) for r in tab(f"{D}/t.bin_counts{sfx}.tsv.gz")}
    tx = {r[0]: tuple(map(int, r[1:])) for r in tab(f"{D}/t.transcript_totals{sfx}.tsv.gz")}
    x3 = {r[1]: (r[7], int(r[8]), int(r[9])) for r in tab(f"{D}/t.excluded_3prime{sfx}.tsv.gz")}
    res = dict(accounting={k: acc.get(k, 0) for k in exp_acc} == exp_acc and acc.get("closes", True), bins=bc == exp_bins,
               transcripts=tx == exp_tx, excluded=x3 == exp_x3)
    print(tag, res); ok &= all(res.values())
    if not all(res.values()): print("  got", acc, bc, tx, x3, sep="\n  ")
    files[tag] = [rd(f"{D}/t.{n}{sfx}.tsv.gz") for n in ("bin_counts", "transcript_totals", "ties", "multibin_examples", "overlap_patterns", "excluded_3prime")]
same = "C0" in files and files.get("C0") == files.get("C1"); print("C0_C1_files_identical", same); ok &= same
# r2 (2026-10-09): SENSITIVITY rule 'any' -- exclude a read with ANY aligned base beyond the 3' end (6th argument)
exp_acc_any = {"eligible_NH1_ge18": 10, "excluded_3prime_overhang": 7, "assigned": 3}
exp_bins_any = {"B|bin003": 1, "B|bin001": 1, "A|bin001": 1}                 # x10 (soft clip, 0 beyond), x4 (5' only), x7 (intron)
exp_tx_any = {"A": (1, 1, 0, 2), "B": (2, 2, 0, 3), "C": (0, 0, 0, 2)}
exp_x3_any = dict(exp_x3, x3_less_beyond=("B", 50, 20), x6_minus_less=("C", 30, 10))
fa = {}
for tag, script, sfx in (("C0", "count_bins100.py", ""), ("C1", "count_bins100_c1.py", ".c1")):
    D = f"{W}/{tag}_any"; os.makedirs(D, exist_ok=True)
    p = subprocess.run([sys.executable, f"{S}/{script}", ann, f"{W}/t.bam", "t", "EV", D, "any"], capture_output=True, text=True)
    if p.returncode: print(tag, "any FAILED TO RUN", p.stderr[-1500:]); ok = False; continue
    acc = json.load(open(f"{D}/t.accounting{sfx}.json"))
    bc = {r[0]: int(r[1]) for r in tab(f"{D}/t.bin_counts{sfx}.tsv.gz")}
    tx = {r[0]: tuple(map(int, r[1:])) for r in tab(f"{D}/t.transcript_totals{sfx}.tsv.gz")}
    x3 = {r[1]: (r[7], int(r[8]), int(r[9])) for r in tab(f"{D}/t.excluded_3prime{sfx}.tsv.gz")}
    res = dict(accounting={k: acc.get(k, 0) for k in exp_acc_any} == exp_acc_any and acc.get("overhang_rule") == "any",
               bins=bc == exp_bins_any, transcripts=tx == exp_tx_any, excluded=x3 == exp_x3_any)
    print(tag, "rule=any", res); ok &= all(res.values())
    if not all(res.values()): print("  got", acc, bc, tx, x3, sep="\n  ")
    fa[tag] = [rd(f"{D}/t.{n}{sfx}.tsv.gz") for n in ("bin_counts", "transcript_totals", "excluded_3prime")]
same = "C0" in fa and fa.get("C0") == fa.get("C1"); print("rule=any C0_C1_files_identical", same); ok &= same
print("ALL_3PRIME_TESTS_PASS" if ok else "3PRIME_TESTS_FAILED"); sys.exit(0 if ok else 1)
