#!/usr/bin/env python3
"""build_extra_annotation.py -- FL r3 (2026-10-09). Annotation only, read-only inputs.
Builds the extra full-length units requested by Justine (2026-10-09: "run the analysis on all RNAs with meaningful biological
labels"; ncRNA rule = PCG rule; no 18S/5.8S/25S; units left out on purpose stay out):
  (1) the 26 GTF transcripts with no CDS lines that NCBI V3 (GCA_002759435.3) calls mRNA: CDS-only frame from V3 CDS
      (P2 output), unmerged 100-nt and 150-nt bins from the CDS 5' end (spliced coordinates), final bin at its true length;
  (2) the 13 V3 5S rRNA genes (official V3 coordinates and length), pooled later as ONE unit '5S_rRNA'.
Also reports any same-strand overlap of the new CDS models with existing binned models (ordinary or competitor-only),
because such overlap would change the ambiguity status of existing reads.
Usage: build_extra_annotation.py <v3_cds_intervals.tsv> <V3 gff.gz> <V3 assembly_report> <bin100 master (unmerged)>
                                 <bin150 master> <outdir>
Coordinates 0-based half-open."""
import sys, os, csv, gzip, re, collections
CDS, GFF, AR, A100, A150, OUT = sys.argv[1:7]; os.makedirs(OUT, exist_ok=True)
acc2chr = {}
for l in open(AR):
    if l.startswith("#"): continue
    f = l.rstrip("\n").split("\t"); acc2chr[f[4]] = "Chr%s_C_auris_B8441" % f[2]
cds = collections.defaultdict(list); txid = {}
for r in csv.DictReader(open(CDS), delimiter="\t"):
    cds[r["gene"]].append((r["chrom"], int(r["start0"]), int(r["end"]), r["strand"])); txid[r["gene"]] = r["transcript"]
assert len(cds) == 26, len(cds)
HDR = ["bin_id", "rna_class", "model_id", "gene_id", "chrom", "strand", "role", "transcript_length", "n_bins", "bin_number",
       "tx_start0", "tx_end", "bin_length", "genomic_blocks_5to3"]
def bins_for(blocks, strand, size):
    """blocks: genomic [s,e) of the CDS; returns list of (tx_start0, tx_end, [genomic blocks 5'->3'])"""
    bl = sorted(blocks, reverse=(strand == "-")); L = sum(e - s for s, e in bl); out = []
    for t0 in range(0, L, size):
        t1 = min(t0 + size, L); gb = []; pos = 0
        for s, e in bl:
            n = e - s; lo, hi = max(t0, pos), min(t1, pos + n)
            if lo < hi:
                if strand == "+": gb.append((s + lo - pos, s + hi - pos))
                else: gb.append((e - (hi - pos), e - (lo - pos)))
            pos += n
        out.append((t0, t1, gb))
    return L, out
for size, name in ((100, "new_bins100.tsv"), (150, "new_bins150.tsv")):
    with open(f"{OUT}/{name}", "w") as w:
        w.write("\t".join(HDR) + "\n")
        for g in sorted(cds):
            ch = {c for c, _, _, _ in cds[g]}; st = {s for _, _, _, s in cds[g]}; assert len(ch) == len(st) == 1
            ch, st = ch.pop(), st.pop(); L, bb = bins_for([(s, e) for _, s, e, _ in cds[g]], st, size)
            assert L % 3 == 0 and sum(t1 - t0 for t0, t1, _ in bb) == L
            for i, (t0, t1, gb) in enumerate(bb, 1):
                w.write("\t".join(map(str, [f"{txid[g]}|bin{i:03d}", "PCG_CDS", txid[g], g, ch, st, "ordinary", L, len(bb), i, t0, t1,
                                            t1 - t0, ";".join(f"{s}-{e}" for s, e in gb)])) + "\n")
# 5S loci (V3 official)
n5 = 0
with gzip.open(GFF, "rt") as h, open(f"{OUT}/rRNA_5S_loci_V3.bed", "w") as w:
    w.write("#chrom\tstart0\tend\tname\tlength\tstrand\tv2_id\n")
    for l in h:
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "rRNA" or "product=5S ribosomal RNA" not in f[8]: continue
        lt = re.search(r"locus_tag=([^;]+)", f[8]).group(1); v2 = re.search(r"B8441v2:([^;]+)", f[8]).group(1)
        s, e = int(f[3]) - 1, int(f[4]); n5 += 1
        w.write(f"{acc2chr[f[0]]}\t{s}\t{e}\t{lt}\t{e - s}\t{f[6]}\t{v2}\n")
assert n5 == 13, n5
# same-strand overlap of new CDS with existing binned models
ex = collections.defaultdict(list)
for A in (A100, A150):
    for r in csv.DictReader(open(A), delimiter="\t"):
        for blk in r["genomic_blocks_5to3"].split(";"):
            s, e = map(int, blk.split("-")); ex[(r["chrom"], r["strand"])].append((s, e, r["model_id"], r["role"]))
with open(f"{OUT}/new_model_overlaps_with_existing.tsv", "w") as w:
    w.write("new_gene\texisting_model\trole\toverlap_bp_max_block\n"); n = 0
    for g in sorted(cds):
        hits = {}
        for ch, s, e, st in cds[g]:
            for s2, e2, m, role in ex[(ch, st)]:
                o = min(e, e2) - max(s, s2)
                if o > 0: hits[(m, role)] = max(hits.get((m, role), 0), o)
        for (m, role), o in sorted(hits.items()): w.write(f"{g}\t{m}\t{role}\t{o}\n"); n += 1
print(f"new PCG models: 26; 5S loci: {n5}; same-strand overlaps with existing binned models: {n}")
