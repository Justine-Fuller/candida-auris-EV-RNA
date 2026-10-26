#!/usr/bin/env python3
"""coverage_track.py -- 2026-10-13, READ-ONLY figure input: per-base read depth along one transcript (RNA 5'->3', CDS blocks of the
bin100 merged annotation) in one library. Frozen v2 read rules: primary, NH = 1, read length >= 18, read strand = transcript strand
(EV / Cell sRNA are sense-stranded). Depth = aligned bases (pysam blocks) per position; also the library's number of eligible reads
(primary, NH = 1, >= 18 nt) for per-million normalisation.
Usage: coverage_track.py <bam> <library> <bin100m_master_annotation.tsv> <transcript_id> <out.tsv>"""
import sys, csv, pysam
BAM, LIB, ANN, TX, OUT = sys.argv[1:6]
rows = [r for r in csv.DictReader(open(ANN), delimiter="\t") if r["model_id"] == TX]
rows.sort(key=lambda r: int(r["bin_number"])); ch, st = rows[0]["chrom"], rows[0]["strand"]
blocks = []
for r in rows:
    for b in r["genomic_blocks_5to3"].split(";"):
        s, e = map(int, b.split("-")); blocks.append((s, e))
pos = []   # transcript coordinate -> genomic position, 5'->3'
for s, e in blocks: pos += list(range(s, e)) if st == "+" else list(range(e - 1, s - 1, -1))
g2t = {g: i for i, g in enumerate(pos)}; lo, hi = min(pos), max(pos) + 1
depth = [0] * len(pos); elig = 0
with pysam.AlignmentFile(BAM, "rb", check_sq=False) as f:
    for a in f.fetch(until_eof=True):
        if a.is_unmapped or a.is_secondary or a.is_supplementary or a.get_tag("NH") != 1: continue
        if (a.infer_read_length() or a.query_length) < 18: continue
        elig += 1
        if a.reference_name != ch or ("-" if a.is_reverse else "+") != st or a.reference_end <= lo or a.reference_start >= hi: continue
        for s, e in a.get_blocks():
            for g in range(max(s, lo), min(e, hi)):
                i = g2t.get(g)
                if i is not None: depth[i] += 1
with open(OUT, "w") as w:
    w.write(f"#library\t{LIB}\ttranscript\t{TX}\teligible_reads\t{elig}\n" + "tx_pos\tdepth\n")
    for i, d in enumerate(depth): w.write(f"{i}\t{d}\n")
print(LIB, elig, sum(depth))
