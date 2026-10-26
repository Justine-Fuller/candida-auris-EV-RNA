#!/usr/bin/env python3
"""Known-answer test for census.py r2. Usage: test_census.py <scriptdir> <workdir>
r2: q3 records are written NON-adjacent (as in a coordinate-sorted BAM); a polyA run (strand flipped) is also checked."""
import sys, os, subprocess, json, pysam
S, W = sys.argv[1], sys.argv[2]; os.makedirs(W, exist_ok=True)
with open(f"{W}/iv.tsv", "w") as o:
    o.write("chrom\tstart\tend\tclass\tgene\ttranscript\tstrand\ttx_offset\n")
    for r in [("PCG_CDS", "P", 1001, 1200, "+"), ("PCG_UTR", "P", 1201, 1300, "+"), ("tRNA", "T1", 3001, 3072, "+"),
              ("rRNA", "R", 5001, 5118, "+"), ("LNC", "LNCRNA.1946", 6001, 7000, "-"), ("LNC", "L", 8001, 9000, "+"),
              ("OTHER_NC", "O", 10001, 10500, "+"), ("STRINGTIE", "S", 11001, 11500, "+")]:
        o.write(f"c1\t{r[2]}\t{r[3]}\t{r[0]}\t{r[1]}\t{r[1]}-T\t{r[4]}\t0\n")
with open(f"{W}/t168.bed", "w") as o:
    o.write("#chrom\tstart0\tend\tname\tscore\tstrand\tctg_family\tfungidb_gene_id\tintron\ttier\n")
    o.write("c1\t3000\t3072\tT1|Gly/GCC\t70\t+\tGly/GCC\tT1\tNA\tA_core\nc1\t4000\t4072\tT2|Gly/GCC\t70\t+\tGly/GCC\tNA\tNA\tA_core\n")
hdr = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "c1", "LN": 20000}]}
reads = [("q1", [(1050, 0, 0)], 30, 1), ("q2", [(1210, 0, 0)], 30, 1), ("q3", [(3010, 0, 0), (4010, 0, 256)], 30, 2),
         ("q4", [(3010, 0, 0), (1050, 0, 256)], 30, 2), ("q5", [(1050, 16, 0)], 30, 1), ("q6", [(15000, 0, 0)], 30, 1),
         ("q7", [(1050, 0, 0)], 15, 1), ("q8", None, 30, 0), ("q9", [(5010, 0, 0)], 30, 1), ("q10", [(6010, 16, 0)], 30, 1),
         ("q11", [(8010, 0, 0)], 30, 1), ("q12", [(10010, 0, 0)], 30, 1), ("q13", [(11010, 0, 0)], 30, 1), ("q14", [(4010, 0, 0)], 30, 1)]
recs = []
for q, alns, L, nh in reads:
    if alns is None: recs.append((q, None, L, nh)); continue
    for al in alns: recs.append((q, al, L, nh))
recs = recs[:3] + recs[4:9] + [recs[3]] + recs[9:]                  # move q3's second record away from its first
with pysam.AlignmentFile(f"{W}/t.bam", "wb", header=hdr) as f:
    for q, al, L, nh in recs:
        alns = None if al is None else [al]
        if alns is None:
            a = pysam.AlignedSegment(); a.query_name = q; a.flag = 4; a.query_sequence = "A" * L; a.reference_id = -1; a.reference_start = -1; f.write(a); continue
        for pos, rev, sec in alns:
            a = pysam.AlignedSegment(); a.query_name = q; a.cigarstring = f"{L}M"; a.query_sequence = "A" * L
            a.flag = rev | sec; a.reference_id = 0; a.reference_start = pos; a.set_tag("NH", nh); f.write(a)
p = subprocess.run([sys.executable, f"{S}/census.py", f"{W}/t.bam", "T", f"{W}/iv.tsv", f"{W}/t168.bed", W, "sRNA"], capture_output=True, text=True)
if p.returncode: print(p.stderr); print("CENSUS_TEST_FAILED"); sys.exit(1)
o = json.load(open(f"{W}/T.census.json")); a = o["accounting"]
exp_acc = {"records": 16, "unmapped_records": 1, "mapped_reads_all_lengths": 13, "mapped_reads_lt18": 1, "reads_ge18_NH1": 10,
           "reads_ge18_NHgt1": 2, "primary_records": 13}
exp_prim = {"PCG_CDS|NH1": 1, "PCG_UTR|NH1": 1, "tRNA|NHgt1": 2, "tRNA|NH1": 1, "antisense_only|NH1": 1, "intergenic|NH1": 1,
            "rRNA_5S|NH1": 1, "rDNA_lncRNA|NH1": 1, "lncRNA|NH1": 1, "OTHER_NC|NH1": 1, "StringTie_only|NH1": 1}
checks = {"accounting": all(a.get(k, 0) == v for k, v in exp_acc.items()) and a["closes"] and not a.get("reads_records_ne_NH"),
          "primary_class": o["primary_class"] == exp_prim,
          "multi_class": o["multi_class_combinations"] == {"tRNA+PCG_CDS|NHgt1": 1},
          "trna_membership": o["trna_membership"] == {"NHgt1|GTF=True|168=True": 2, "NH1|GTF=False|168=True": 1},
          "any_overlap_PCG_CDS": o["any_overlap_class"].get("PCG_CDS|NH1") == 1 and o["any_overlap_class"].get("PCG_CDS|NHgt1") == 1}
# polyA: same BAM, strand flipped -> q5 (reverse read on the + gene) becomes PCG_CDS; q1 (forward) becomes antisense_only, etc.
os.makedirs(f"{W}/pa", exist_ok=True)
p2 = subprocess.run([sys.executable, f"{S}/census.py", f"{W}/t.bam", "T", f"{W}/iv.tsv", f"{W}/t168.bed", f"{W}/pa", "polyA"], capture_output=True, text=True)
pa = json.load(open(f"{W}/pa/T.census.json"))["primary_class"] if p2.returncode == 0 else {}
checks["polyA_strand_flip"] = (pa.get("PCG_CDS|NH1") == 1 and pa.get("rDNA_lncRNA|NH1") is None and pa.get("lncRNA|NH1") is None
                               and pa.get("antisense_only|NH1", 0) == 8 and pa.get("antisense_only|NHgt1", 0) == 2 and pa.get("intergenic|NH1") == 1)
for k, v in checks.items(): print(k, v)
if not all(checks.values()): print(json.dumps(o, indent=1))
print("ALL_CENSUS_TESTS_PASS" if all(checks.values()) else "CENSUS_TEST_FAILED"); sys.exit(0 if all(checks.values()) else 1)
