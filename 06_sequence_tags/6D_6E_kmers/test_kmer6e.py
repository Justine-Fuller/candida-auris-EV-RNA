#!/usr/bin/env python3
"""test_kmer6e.py -- Step 6E known-answer tests (2026-10-12). Usage: test_kmer6e.py <scripts dir> <empty workdir>
T1 kmer_reads6e.py on a synthetic BAM: plus / minus strand PCG, antisense, NH = 2, 17 nt, two-model ambiguity, competitor-only,
   lncRNA with soft clip (kept), FL r3 new CDS model, tRNA multimapper in one family (kept), cross-family (excluded), minus-strand
   tRNA, antisense tRNA, 5S multimapper with a SEQ-less secondary record, read with N, unmapped read.
T2 kmer_match6e.py: "__all" 6-mer counts equal a brute-force count (N windows skipped, reads < 6 nt give none); matched Cell
   GC-bin read numbers are floor(n * EV target) with n maximal; unmatchable bins dropped from both sides; reproducible with seed."""
import sys, os, json, gzip, subprocess, random, re, itertools, collections
import pysam
S, W = sys.argv[1:3]; os.makedirs(W, exist_ok=False)
PY = sys.executable; fails = []
def check(name, ok, detail=""):
    print(("PASS " if ok else "FAIL ") + name, detail); (None if ok else fails.append(name))
RC = str.maketrans("ACGTN", "TGCAN"); rc = lambda s: s.translate(RC)[::-1]
random.seed(1)
rs = lambda n: "".join(random.choice("ACGT") for _ in range(n))
# ---------- T1 synthetic annotation
H = "bin_id\trna_class\tmodel_id\tgene_id\tchrom\tstrand\trole\tgenomic_blocks_5to3\n"
ann = [("PCGP|bin001", "PCG_CDS", "PCGP", "+", "ordinary", "1000-1300"), ("PCGM|bin001", "PCG_CDS", "PCGM", "-", "ordinary", "2000-2300"),
       ("LNC|bin001", "lncRNA", "LNC", "+", "ordinary", "3000-3300"), ("COMP|bin001", "lncRNA", "COMP", "+", "competitor_only", "4000-4300"),
       ("PCGP2|bin001", "PCG_CDS", "PCGP2", "+", "ordinary", "1250-1500")]
with open(f"{W}/ann.tsv", "w") as w:
    w.write(H); [w.write(f"{b}\t{c}\t{m}\t{m}\tchrT\t{st}\t{r}\t{blk}\n") for b, c, m, st, r, blk in ann]
with open(f"{W}/new.tsv", "w") as w:
    w.write(H); w.write("NEWP|bin001\tPCG_CDS\tNEWP\tNEWP\tchrT\t+\tordinary\t5000-5300\n")
with open(f"{W}/trna.bed", "w") as w:
    w.write("#chrom\tstart0\tend\tname\tscore\tstrand\tctg_family\tfungidb_gene_id\tintron\ttier\n")
    for i, (s, e, st, f) in enumerate([(6000, 6072, "+", "Ala/AGC"), (7000, 7072, "+", "Ala/AGC"), (8000, 8072, "+", "Gly/GCC"), (9000, 9082, "-", "Leu/CAA")]):
        w.write(f"chrT\t{s}\t{e}\tt{i}|{f}\t50\t{st}\t{f}\tNA\tNA\tA_core\n")
with open(f"{W}/r5.bed", "w") as w:
    w.write("#chrom\tstart0\tend\tname\tlength\tstrand\tv2_id\nchrT\t10000\t10118\tS1\t118\t+\tx\nchrT\t11000\t11118\tS2\t118\t-\tx\n")
hdr = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}, "SQ": [{"SN": "chrT", "LN": 20000}]})
recs = []
def rec(q, orig, pos, rev, cigar, nh, sec=False, noseq=False, unmapped=False):
    a = pysam.AlignedSegment(hdr); a.query_name = q
    if unmapped: a.flag = 4; a.query_sequence = orig; a.set_tag("NH", 0); recs.append(a); return
    a.flag = (16 if rev else 0) | (256 if sec else 0); a.reference_id = 0; a.reference_start = pos; a.cigarstring = cigar
    a.mapping_quality = 255
    if not noseq: a.query_sequence = rc(orig) if rev else orig
    a.set_tag("NH", nh); recs.append(a)
o = {k: rs(25) for k in range(1, 17)}
o[5] = rs(17); o[15] = o[15][:10] + "N" + o[15][11:]
rec("r1", o[1], 1100, False, "25M", 1)
rec("r2", o[2], 2100, True, "25M", 1)
rec("r3", o[3], 2100, False, "25M", 1)
rec("r4", o[4], 1100, False, "25M", 2); rec("r4", o[4], 5100, False, "25M", 2, sec=True)
rec("r5", o[5], 1100, False, "17M", 1)
rec("r6", o[6], 1240, False, "25M", 1)
rec("r7", o[7], 4100, False, "25M", 1)
rec("r8", o[8], 3100, False, "3S22M", 1)
rec("r9", o[9], 5100, False, "25M", 1)
rec("r10", o[10], 6010, False, "25M", 2); rec("r10", o[10], 7010, False, "25M", 2, sec=True)
rec("r11", o[11], 6010, False, "25M", 2); rec("r11", o[11], 8010, False, "25M", 2, sec=True)
rec("r12", o[12], 9010, True, "25M", 1)
rec("r13", o[13], 10010, False, "25M", 2); rec("r13", o[13], 11010, True, "25M", 2, sec=True, noseq=True)
rec("r14", o[14], 6010, True, "25M", 1)
rec("r15", o[15], 1100, False, "25M", 1)
rec("r16", o[16], 0, False, "", 0, unmapped=True)
with pysam.AlignmentFile(f"{W}/t.bam", "wb", header=hdr) as f:
    for a in recs: f.write(a)
os.makedirs(f"{W}/k1")
p = subprocess.run([PY, f"{S}/kmer_reads6e.py", f"{W}/t.bam", "TLIB", f"{W}/ann.tsv", f"{W}/new.tsv", f"{W}/trna.bed", f"{W}/r5.bed", f"{W}/k1"],
                   capture_output=True, text=True); print(p.stdout, p.stderr)
d = json.load(gzip.open(f"{W}/k1/TLIB.k1.json.gz", "rt")); T = d["types"]; acc = d["accounting"]
check("T1_PCG_CDS", T["PCG_CDS"] == {o[1]: 1, o[2]: 1, o[9]: 1, o[15]: 1}, T["PCG_CDS"])
check("T1_lncRNA_softclip_kept", T["lncRNA"] == {o[8]: 1})
check("T1_tRNA", T["tRNA"] == {o[10]: 1, o[12]: 1}, T["tRNA"])
check("T1_5S", T["5S_rRNA"] == {o[13]: 1})
check("T1_accounting", acc.get("model_ambiguous") == 1 and acc.get("model_competitor_only") == 1 and acc.get("tRNA_cross_family_excluded") == 1
      and acc.get("primary_mapped") == 15 and acc.get("eligible_NH1_ge18") == 10, acc)
# ---------- T2 matching + counting
WORDS = ["".join(x) for x in itertools.product("ACGT", repeat=6)]
def brute(c):
    v = collections.Counter()
    for s, n in c.items():
        for i in range(len(s) - 5):
            k = s[i:i + 6]
            if re.fullmatch("[ACGT]{6}", k): v[k] += n
    return v
os.makedirs(f"{W}/k1b"); man = ["library\tcompartment\tgroup"]; data = {}
def gcread(gc, L):
    return "".join(random.choice("GC") if random.random() < gc else random.choice("AT") for _ in range(L))
for comp, gcs in (("EV", (0.3, 0.5, 0.6)), ("Cell_sRNA", (0.35, 0.45, 0.55))):
    for i in (1, 2, 3):
        L = f"G_{comp}_{i}"; man.append(f"{L}\t{comp}\tG"); c = collections.Counter()
        for _ in range(3000): c[gcread(random.choice(gcs), random.randint(18, 30))] += random.choice((1, 1, 2))
        c["ACGTNACGTACG"] += 3; c["ACGTA"] += 5; c["NNNNNNNNNNNNNNNNNN"] += 1
        if comp == "EV" and i == 1: c["GGGGGGGGGGGGGGGGGGGG"] += 4      # EV-only GC bin 9 -> must be dropped if a Cell lib lacks it
        data[L] = c
        json.dump(dict(library=L, accounting={}, types={"PCG_CDS": dict(c), "lncRNA": {}, "tRNA": {}, "5S_rRNA": {}}), gzip.open(f"{W}/k1b/{L}.k1.json.gz", "wt"))
open(f"{W}/man.tsv", "w").write("\n".join(man) + "\n")
for run in ("k2a", "k2b"):
    p = subprocess.run([PY, f"{S}/kmer_match6e.py", f"{W}/man.tsv", f"{W}/k1b", f"{W}/{run}", "2"], capture_output=True, text=True)
    print(p.stdout[-500:], p.stderr[-2000:])
import csv
rows = list(csv.DictReader(gzip.open(f"{W}/k2a/G.PCG_CDS.counts.tsv.gz", "rt"), delimiter="\t"))
ok = len(rows) == 4096
for L, c in data.items():
    b = brute(c); ok = ok and all(int(r[f"{L}__all"]) == b.get(r["word"], 0) for r in rows)
check("T2_all_counts_equal_bruteforce", ok)
check("T2_reproducible_seed", gzip.open(f"{W}/k2a/G.PCG_CDS.counts.tsv.gz").read() == gzip.open(f"{W}/k2b/G.PCG_CDS.counts.tsv.gz").read())
gc = list(csv.DictReader(open(f"{W}/k2a/G.PCG_CDS.gc.tsv"), delimiter="\t"))
P = {int(r["gc_bin"]): float(r["value"]) for r in gc if r["library"] == "EV_target"}
def gcb(s):
    a = sum(s.count(x) for x in "ACGT"); return None if a == 0 else min(int(10 * (s.count("G") + s.count("C")) / a), 9)
ok = True; det = []
for i in (1, 2, 3):
    L = f"G_Cell_sRNA_{i}"; Nb = collections.Counter()
    for s, n in data[L].items():
        b = gcb(s)
        if b is not None: Nb[b] += n
    m = {int(r["gc_bin"]): int(r["value"]) for r in gc if r["library"] == L and r["set"] == "matched_draw1"}
    n = min(Nb[b] / P[b] for b in P if P[b] > 0); n = int(n)
    exp = {b: int(n * P[b] + 1e-9) for b in P if int(n * P[b] + 1e-9) > 0}
    ok = ok and m == exp; det.append((L, n, sum(m.values())))
check("T2_matched_bins_floor_n_times_target", ok, det)
keep = {int(r["gc_bin"]) for r in gc if r["library"] == "G_EV_1" and r["set"] == "matchable"}
n9 = sum(n for s_, n in data["G_EV_1"].items() if gcb(s_) == 9)
check("T2_EV_only_bin9_dropped", 9 not in keep and P.get(9, 0) == 0 and any(r["set"] == "dropped_unmatchable_bin" and r["gc_bin"] == "9"
      and r["library"] == "G_EV_1" and int(r["value"]) == n9 for r in gc), (sorted(keep), n9))
bm = brute({s_: n for s_, n in data["G_EV_1"].items() if gcb(s_) in keep})
check("T2_EV_matchable_counts", all(int(r["G_EV_1__m"]) == bm.get(r["word"], 0) for r in rows))
print("ALL_6E_COUNTER_TESTS_PASS" if not fails else f"6E_COUNTER_TESTS_FAILED {fails}")
