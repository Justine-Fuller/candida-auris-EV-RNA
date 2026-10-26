#!/usr/bin/env python3
"""s6f_gene_counts.py -- Step 6F (v2, 2026-10-13), stage F2. READ-ONLY per-library read count per gene model (EV / Cell sRNA).
Read rule = the Step 6E PCG/lncRNA rule (kmer_reads6e.py), unchanged: primary, mapped, NH == 1, read length incl. soft clips
>= 18 nt; touched models = models with >= 1 aligned base (pysam blocks) on their 100-nt-bin exonic positions on the read's
strand (sense; EV / Cell sRNA are sense-stranded); exactly 1 touched model and role != competitor_only -> +1 to that model;
0 -> none; > 1 -> ambiguous (not used); competitor-only -> not used. One read = one count (physical read, NH = 1).
Gene models = bin100 master annotation (PCG CDS frame + lncRNA) + the 26 FL r3 CDS models.
Usage: s6f_gene_counts.py <bam> <library> <bin100 master annotation> <FL r3 new_bins100.tsv> <outdir>
Output: <lib>.gene_counts.tsv (model_id, rna_class, reads; every model, zeros included) + <lib>.accounting.json."""
import sys, json, csv, collections, pysam
BAM, LIB, A100, N100, OUT = sys.argv[1:6]
B = 10000
midx = collections.defaultdict(list); bmod = []; mclass = {}; mrole = {}
for p in (A100, N100):
    for r in csv.DictReader(open(p), delimiter="\t"):
        k = len(bmod); bmod.append(r["model_id"]); mclass[r["model_id"]] = r["rna_class"]; mrole[r["model_id"]] = r["role"]
        for blk in r["genomic_blocks_5to3"].split(";"):
            s, e = map(int, blk.split("-"))
            for b in range(s // B, (e - 1) // B + 1): midx[(r["chrom"], r["strand"], b)].append((s, e, k))
assert set(mclass.values()) <= {"PCG_CDS", "lncRNA"}
def models(a):
    hit = set(); st = "-" if a.is_reverse else "+"
    for s, e in a.get_blocks():
        for b in range(s // B, (e - 1) // B + 1):
            for bs, be, k in midx.get((a.reference_name, st, b), ()):
                if min(e, be) > max(s, bs): hit.add(bmod[k])   # set: a hit seen via two buckets counts once
    return hit
acc = collections.Counter(); C = collections.Counter()
with pysam.AlignmentFile(BAM, "rb", check_sq=False) as f:
    for r in f.fetch(until_eof=True):
        acc["records"] += 1
        if r.is_unmapped or r.is_supplementary or r.is_secondary: continue
        acc["primary_mapped"] += 1
        if r.get_tag("NH") != 1 or (r.infer_read_length() or r.query_length) < 18: continue
        acc["eligible_NH1_ge18"] += 1
        ms = models(r)
        if not ms: acc["model_none"] += 1
        elif len(ms) > 1: acc["model_ambiguous"] += 1
        else:
            m = ms.pop()
            if mrole[m] == "competitor_only": acc["model_competitor_only"] += 1
            else: C[m] += 1; acc[mclass[m] + "_reads"] += 1
with open(f"{OUT}/{LIB}.gene_counts.tsv", "w") as w:
    w.write("model_id\trna_class\treads\n")
    for m in sorted(mclass): w.write(f"{m}\t{mclass[m]}\t{C[m]}\n")
json.dump(dict(library=LIB, bam=BAM, accounting=dict(acc)), open(f"{OUT}/{LIB}.accounting.json", "w"), sort_keys=True, indent=1)
print(json.dumps(dict(acc), sort_keys=True))
