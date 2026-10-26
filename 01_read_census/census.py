#!/usr/bin/env python3
"""census.py -- READ-ONLY read census by RNA class (2026-10-09). Where do ALL mapped reads of an EV / Cell sRNA library
fall, including multimappers (NH > 1)? Descriptive only: no rule, no filter decision, nothing frozen is changed.
Usage: census.py <bam> <library> <intervals.tsv (1-based)> <trna168.bed (0-based)> <outdir> <sRNA|polyA>
r2 (2026-10-09): records are grouped per QNAME until NH records are collected (works for unsorted AND coordinate-
sorted BAMs; leftovers counted as reads_records_ne_NH); polyA (TruSeq stranded mRNA, reverse-stranded) flips the read strand.
Unit = one physical read (QNAME). All alignment records of a read (STAR outSAMmultNmax -1, Unsorted: records of a read are
adjacent; checked) are pooled. Read length = infer_read_length (incl. soft clips); main tables use reads >= 18 nt (frozen
v2 length rule); the all-length total is kept for the accounting gate. Strand: EV and Cell sRNA are sense-stranded, so a
feature is "same strand" when its strand equals the read's alignment strand. Overlap = >= 1 aligned base (pysam blocks,
D/N excluded) inside the feature span.
Classes (same strand). A read gets the SET of classes hit by any of its alignments. For reporting, a primary class is the
first class in this fixed REPORTING order (descriptive convenience, not an analysis rule):
  tRNA > rRNA_5S > rDNA_lncRNA (LNCRNA.1946/.1057) > OTHER_NC > PCG_CDS > PCG_UTR > lncRNA > StringTie_only
  > antisense_only (a feature only on the opposite strand) > intergenic (no feature on either strand)
tRNA = union of the GTF tRNA intervals (old 120) and the tRNAscan 168 set; both memberships are also reported separately.
Outputs: <lib>.census.json"""
import sys, json, collections, pysam
BAM, LIB, IV, TB, OUT, KIND = sys.argv[1:7]
assert KIND in ("sRNA", "polyA")
B = 10000
ORDER = ["tRNA", "rRNA_5S", "rDNA_lncRNA", "OTHER_NC", "PCG_CDS", "PCG_UTR", "lncRNA", "StringTie_only"]
idx = collections.defaultdict(list)
def add(ch, s, e, st, cls):
    for b in range(s // B, (e - 1) // B + 1): idx[(ch, b)].append((s, e, st, cls))
with open(IV) as f:
    next(f)
    for line in f:
        ch, s, e, cls, gene, tx, st = line.rstrip("\n").split("\t")[:7]
        s0, e0 = int(s) - 1, int(e)
        c = {"tRNA": "tRNA_GTF", "rRNA": "rRNA_5S", "OTHER_NC": "OTHER_NC", "PCG_CDS": "PCG_CDS", "PCG_UTR": "PCG_UTR",
             "STRINGTIE": "StringTie_only", "LNC": "rDNA_lncRNA" if gene in ("LNCRNA.1946", "LNCRNA.1057") else "lncRNA"}[cls]
        add(ch, s0, e0, st, c)
with open(TB) as f:
    next(f)
    for line in f:
        r = line.rstrip("\n").split("\t"); add(r[0], int(r[1]), int(r[2]), r[5], "tRNA_168")
def classes(aln):
    rs = "-" if aln.is_reverse else "+"
    if KIND == "polyA": rs = "+" if rs == "-" else "-"
    same, anti = set(), set(); seen = set()
    for s, e in aln.get_blocks():
        for b in range(s // B, (e - 1) // B + 1):
            for it in idx.get((aln.reference_name, b), ()):
                if it in seen: continue
                fs, fe, fst, c = it
                if fs < e and s < fe:
                    seen.add(it)
                    (same if fst in (rs, ".") else anti).add(c)
    return same, anti
acc = collections.Counter(); prim = collections.Counter(); anyc = collections.Counter(); lens = collections.Counter()
combo = collections.Counter(); trna = collections.Counter(); nhdist = collections.Counter()
pending = {}
def flush(q, rs):
    acc["mapped_reads_all_lengths"] += 1
    nh = rs[0].get_tag("NH"); L = rs[0].infer_read_length() or rs[0].query_length
    if len(rs) != nh: acc["reads_records_ne_NH"] += 1
    if L < 18: acc["mapped_reads_lt18"] += 1; return
    k = "NH1" if nh == 1 else "NHgt1"; acc[f"reads_ge18_{k}"] += 1
    same, anti = set(), set()
    for a in rs:
        s_, a_ = classes(a); same |= s_; anti |= a_
    t = {"tRNA" if c in ("tRNA_GTF", "tRNA_168") else c for c in same}
    if "tRNA_GTF" in same or "tRNA_168" in same:
        trna[(k, "tRNA_GTF" in same, "tRNA_168" in same)] += 1
    p = next((c for c in ORDER if c in t), "antisense_only" if anti else "intergenic")
    prim[(p, k)] += 1; lens[(p, min(L, 61))] += 1
    for c in t: anyc[(c, k)] += 1
    if len(t) > 1: combo[("+".join(c for c in ORDER if c in t), k)] += 1
    nhdist[(p, min(nh, 21))] += 1
with pysam.AlignmentFile(BAM, "rb", check_sq=False) as f:
    for r in f.fetch(until_eof=True):
        acc["records"] += 1
        if r.is_unmapped: acc["unmapped_records"] += 1; continue
        if r.is_supplementary: acc["supplementary_records"] += 1; continue
        if not r.is_secondary: acc["primary_records"] += 1
        q = r.query_name; pending.setdefault(q, []).append(r)
        if len(pending[q]) == r.get_tag("NH"): flush(q, pending.pop(q))
    for q, rs in pending.items(): flush(q, rs)                      # incomplete groups (records < NH), counted above
out = dict(library=LIB, bam=BAM, accounting=dict(acc),
           primary_class={f"{c}|{k}": n for (c, k), n in sorted(prim.items())},
           any_overlap_class={f"{c}|{k}": n for (c, k), n in sorted(anyc.items())},
           multi_class_combinations={f"{c}|{k}": n for (c, k), n in sorted(combo.items())},
           trna_membership={f"{k}|GTF={g}|168={s}": n for (k, g, s), n in sorted(trna.items())},
           length_by_primary_class={f"{c}|{l}": n for (c, l), n in sorted(lens.items())},
           nh_by_primary_class={f"{c}|{n_}": n for (c, n_), n in sorted(nhdist.items())})
out["accounting"]["closes"] = (acc["reads_ge18_NH1"] + acc["reads_ge18_NHgt1"] == sum(prim.values())
                               and acc["mapped_reads_all_lengths"] == acc["reads_ge18_NH1"] + acc["reads_ge18_NHgt1"] + acc["mapped_reads_lt18"]
                               and acc["mapped_reads_all_lengths"] == acc["primary_records"])
json.dump(out, open(f"{OUT}/{LIB}.census.json", "w"), indent=1, sort_keys=True)
print(json.dumps(out["accounting"], sort_keys=True))
