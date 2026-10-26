#!/usr/bin/env python3
"""Known-answer test for trna_cca_recount.py (2026-10-13). Usage: test_cca_recount.py <scripts dir> <trna_count.py> <workdir>
Synthetic genome (one chromosome) with 6 tRNA loci of 72 nt; flanks after the 3' end chosen per case:
  L1 Gly/GCC + flank TTTTT  : (1) soft-clipped CCA (old case)                       -> old CCA, new CCA non-templated
  L2 Thr/CGT + flank CGTTT  : (2) 1 aligned C + soft clip CA (Thr/CGT case)          -> old 0, new CCA non-templated
  L3 Arg/TCT + flank CCATG  : (3) CCA fully aligned over genomic CCA (Chr3.trna3)    -> old 0, new CCA templated-ambiguous
  L4 Ser/CGA + flank TCATT  : (4) CCA aligned over TCA with a mismatch (Ser/CGA)     -> old 0, new CCA non-templated
  L5 Leu/CAA + flank TTTAAG : (5) genomic trailer read (aligned over TTTAAG)         -> not CCA (tail 'other')
  L6 Ala/AGC minus strand, flank (RNA orientation) CGTTT : (2) on the minus strand   -> old 0, new CCA non-templated
  plus: (6) read ending 2 nt short of L1 (no reach), (7) clip CC only at L1 (tail CC), (8) same-family multimapper L2 + L7
  (second Thr/CGT locus, flank TTTTT) with soft-clipped CCA at L7 on a SEQ-less secondary record (old CCA too).
Also: the old-rule fields equal trna_count.py run on the same BAM; INVARIANT_MISMATCH = 0."""
import sys, os, subprocess, json, random, pysam
S, TC, W = sys.argv[1:4]; os.makedirs(W, exist_ok=False)
random.seed(3); RC = str.maketrans("ACGT", "TGCA"); rc = lambda s: s.translate(RC)[::-1]
g = [random.choice("ACGT") for _ in range(20000)]
LOC = [("L1", 1000, "+", "Gly/GCC", "TTTTT"), ("L2", 2000, "+", "Thr/CGT", "CGTTT"), ("L3", 3000, "+", "Arg/TCT", "CCATG"),
       ("L4", 4000, "+", "Ser/CGA", "TCATT"), ("L5", 5000, "+", "Leu/CAA", "TTTAAG"), ("L6", 6000, "-", "Ala/AGC", "CGTTT"),
       ("L7", 7000, "+", "Thr/CGT", "TTTTT")]
for n, s, st, f, fl in LOC:
    if st == "+": g[s + 72:s + 72 + len(fl)] = list(fl)
    else: g[s - len(fl):s] = list(rc(fl))
G = "".join(g); open(f"{W}/g.fa", "w").write(">c1\n" + "\n".join(G[i:i + 60] for i in range(0, len(G), 60)) + "\n")
with open(f"{W}/t.bed", "w") as o:
    o.write("#chrom\tstart0\tend\tname\tscore\tstrand\tctg_family\tfungidb_gene_id\tintron\ttier\n")
    for n, s, st, f, fl in LOC: o.write(f"c1\t{s}\t{s + 72}\t{n}|{f}\t60\t{st}\t{f}\tNA\tNA\tA_core\n")
open(f"{W}/bins.tsv", "w").write("bin_id\tchrom\tstrand\tbin_length\tgenomic_blocks_5to3\nM|bin001\tc1\t+\t100\t15000-15100\n")
hdr = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "c1", "LN": 20000}]}
def gp(s, e): return G[s:e]
R = []   # (qname, [(pos, cigar, reverse, secondary)], SEQ as stored (genome-forward), NH)
R.append(("c1_softCCA", [(1050, "22M3S", 0, 0)], gp(1050, 1072) + "CCA", 1))
R.append(("c2_thr", [(2050, "23M2S", 0, 0)], gp(2050, 2072) + "C" + "CA", 1))
R.append(("c3_templ", [(3050, "25M", 0, 0)], gp(3050, 3072) + "CCA", 1))
R.append(("c4_mism", [(4050, "25M", 0, 0)], gp(4050, 4072) + "CCA", 1))
R.append(("c5_trailer", [(5050, "28M", 0, 0)], gp(5050, 5072) + "TTTAAG", 1))
R.append(("c6_minus", [(5999, "2S23M", 1, 0)], rc("CA") + gp(5999, 6022), 1))   # RNA tail: aligned C (genome G at 5999) + clip CA
R.append(("c7_short", [(1040, "30M", 0, 0)], gp(1040, 1070), 1))
R.append(("c8_CC", [(1050, "22M2S", 0, 0)], gp(1050, 1072) + "CC", 1))
R.append(("c9_mm", [(2050, "22M3S", 0, 0), (7050, "22M3S", 0, 1)], gp(2050, 2072) + "CCA", 2))
with pysam.AlignmentFile(f"{W}/t.bam", "wb", header=hdr) as f:
    for q, alns, seq, nh in R:
        for pos, cig, rev, sec in alns:
            a = pysam.AlignedSegment(); a.query_name = q; a.cigarstring = cig; a.flag = (16 if rev else 0) | (256 if sec else 0)
            a.reference_id = 0; a.reference_start = pos; a.query_sequence = None if sec else seq; a.set_tag("NH", nh); f.write(a)
p = subprocess.run([sys.executable, f"{S}/trna_cca_recount.py", f"{W}/t.bam", "T", f"{W}/t.bed", f"{W}/g.fa", W], capture_output=True, text=True)
print(p.stdout, p.stderr[-2000:])
p2 = subprocess.run([sys.executable, TC, f"{W}/t.bam", "T", f"{W}/bins.tsv", f"{W}/t.bed", W], capture_output=True, text=True); print(p2.stderr[-1000:])
n = json.load(open(f"{W}/T.cca2.json")); o = json.load(open(f"{W}/T.trna.json"))
N, O = n["trna_families"], o["trna_families"]; fails = []
def chk(name, ok, d=""):
    print(("PASS " if ok else "FAIL ") + name, d); ok or fails.append(name)
OLDF = ["reads", "reads_NH1", "reach_3prime", "clipped_at_3prime", "CCA_reads", "reads_bigN_excluded_from_CCA_check"]
chk("old_rule_equals_trna_count", all(N.get(F, {}).get(k, 0) == O[F].get(k, 0) for F in O for k in OLDF) and set(N) == set(O),
    {F: {k: (N.get(F, {}).get(k, 0), O[F].get(k, 0)) for k in OLDF} for F in O})
g_ = lambda F, k: N.get(F, {}).get(k, 0)
chk("1_softclipped_CCA", g_("Gly/GCC", "CCA_reads") == 1 and g_("Gly/GCC", "tail_CCA_non_templated") == 1)
chk("6_short_read_no_reach_7_CC", g_("Gly/GCC", "reads") == 3 and g_("Gly/GCC", "reach_3prime_new") == 2 and g_("Gly/GCC", "tail_CC") == 1)
chk("2_aligned_C_plus_clip_CA", g_("Thr/CGT", "tail_CCA_non_templated") == 2 and g_("Thr/CGT", "CCA_reads") == 1, N.get("Thr/CGT"))
chk("8_multimapper_secondary_CCA", g_("Thr/CGT", "reads") == 2)
chk("3_templated_ambiguous", g_("Arg/TCT", "tail_CCA_templated_ambiguous") == 1 and g_("Arg/TCT", "CCA_reads") == 0)
chk("4_mismatch_TCA", g_("Ser/CGA", "tail_CCA_non_templated") == 1 and g_("Ser/CGA", "CCA_reads") == 0)
chk("5_trailer_not_CCA", g_("Leu/CAA", "tail_CCA") == 0 and g_("Leu/CAA", "tail_other") == 1)
chk("2_minus_strand", g_("Ala/AGC", "tail_CCA_non_templated") == 1 and g_("Ala/AGC", "CCA_reads") == 0, N.get("Ala/AGC"))
chk("invariant", n["accounting"].get("INVARIANT_MISMATCH", 0) == 0)
print("ALL_CCA_RECOUNT_TESTS_PASS" if not fails else f"CCA_RECOUNT_TESTS_FAILED {fails}")
