#!/usr/bin/env python3
"""kmer_reads6e.py -- Step 6E (v2, 2026-10-12), stage K1. READ-ONLY per-library collector of RNA-oriented read sequences per
RNA type (approved plan 00_PLAN_DRAFT_6D_6E_v2_2026-10-12.md: Q3 a = PCG CDS, lncRNA, tRNA family, 5S; no rDNA unit).
Usage: kmer_reads6e.py <bam> <library> <bin100 master annotation (unmerged)> <FL r3 new_bins100.tsv> <trna168.bed>
                       <rRNA_5S_loci_V3.bed> <outdir>
Read rules (frozen v2; same model / locus logic as fl_extra_count.py and trna_count.py):
- every read: read length incl. soft clips (infer_read_length) >= 18 nt.
- PCG_CDS / lncRNA: primary, mapped, NH == 1. Touched models = models with >= 1 aligned base (pysam blocks; D/N excluded) on
  their 100-nt bin exonic positions on the read's strand (sense only; EV / Cell sRNA are sense-stranded). Exactly 1 touched
  model, role != competitor_only -> the read belongs to that model's class (the 26 FL r3 CDS models are PCG_CDS). 0 -> none;
  > 1 -> ambiguous (not used); competitor-only (LNCRNA.1946 / .1057) -> not used. No bin assignment is needed here.
- tRNA: all records of a read are grouped (until NH records are seen); a read is a tRNA read when every alignment touching a
  tRNA locus (same strand, >= 1 aligned base) touches loci of ONE anticodon family (tRNAscan 168 / 43; T1 rule, multimappers
  allowed within a family). Cross-family reads are not used.
- 5S: grouped records; >= 1 alignment touching a V3 5S locus on the same strand (FL r3 pooled-unit rule, multimappers allowed).
- The four types are collected independently (a read may in principle count in two types; counted in accounting).
Sequence: the whole read as sequenced (SEQ incl. soft clips), turned back to the read's own orientation (reverse complement
when the primary record is on the reverse strand). Because only SENSE hits are used, this is the RNA's own 5'->3' sequence.
Output: <lib>.k1.json.gz = {"library", "accounting", "types": {type: {sequence: reads}}} (collapsed identical sequences)."""
import sys, json, gzip, csv, collections, pysam
BAM, LIB, A100, N100, TB, R5, OUT = sys.argv[1:8]
B = 10000
RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")
rc = lambda s: s.translate(RC)[::-1]
midx = collections.defaultdict(list); bmod = []; mclass = {}; mrole = {}
for p in (A100, N100):
    for r in csv.DictReader(open(p), delimiter="\t"):
        k = len(bmod); bmod.append(r["model_id"]); mclass[r["model_id"]] = r["rna_class"]; mrole[r["model_id"]] = r["role"]
        for blk in r["genomic_blocks_5to3"].split(";"):
            s, e = map(int, blk.split("-"))
            for b in range(s // B, (e - 1) // B + 1): midx[(r["chrom"], r["strand"], b)].append((s, e, k))
assert set(mclass.values()) <= {"PCG_CDS", "lncRNA"}
def loadloci(path, famcol):
    loci = []; idx = collections.defaultdict(list)
    for l in open(path):
        if l.startswith("#"): continue
        x = l.rstrip("\n").split("\t"); s, e = int(x[1]), int(x[2])
        k = len(loci); loci.append(dict(ch=x[0], s=s, e=e, st=x[5], fam=x[famcol] if famcol is not None else "5S_rRNA"))
        for b in range(s // B, (e - 1) // B + 1): idx[(x[0], x[5], b)].append(k)
    return loci, idx
tloci, tidx = loadloci(TB, 6)
sloci, sidx = loadloci(R5, None)
strand = lambda a: "-" if a.is_reverse else "+"
def models(a):
    hit = set(); st = strand(a)
    for s, e in a.get_blocks():
        for b in range(s // B, (e - 1) // B + 1):
            for bs, be, k in midx.get((a.reference_name, st, b), ()):
                if min(e, be) > max(s, bs): hit.add(bmod[k])
    return hit
def lochits(a, loci, idx):
    hit = set(); st = strand(a)
    for s, e in a.get_blocks():
        for b in range(s // B, (e - 1) // B + 1):
            for k in idx.get((a.reference_name, st, b), ()):
                if loci[k]["s"] < e and s < loci[k]["e"]: hit.add(k)
    return hit
def readseq(rs):
    p = next((x for x in rs if not x.is_secondary and x.query_sequence), None) or next((x for x in rs if x.query_sequence), None)
    if p is None: return None
    s = p.query_sequence.upper()
    return rc(s) if p.is_reverse else s
acc = collections.Counter(); T = {t: collections.Counter() for t in ("PCG_CDS", "lncRNA", "tRNA", "5S_rRNA")}
used = collections.defaultdict(set)   # qname -> types (only for reads already used by the NH=1 model path, to count overlaps)
def grouped(rs):
    if (rs[0].infer_read_length() or rs[0].query_length) < 18: return
    th = [lochits(a, tloci, tidx) for a in rs]
    fams = {tloci[k]["fam"] for h in th for k in h}
    q = rs[0].query_name
    if fams:
        acc["tRNA_touching_reads"] += 1
        if len(fams) > 1: acc["tRNA_cross_family_excluded"] += 1
        else:
            s = readseq(rs)
            if s is None: acc["tRNA_no_SEQ"] += 1
            else:
                T["tRNA"][s] += 1; acc["tRNA_reads"] += 1
                if q in used: acc["reads_in_tRNA_and_" + "_".join(sorted(used[q]))] += 1
    if any(lochits(a, sloci, sidx) for a in rs):
        s = readseq(rs)
        if s is None: acc["5S_no_SEQ"] += 1
        else:
            T["5S_rRNA"][s] += 1; acc["5S_reads"] += 1
            if q in used: acc["reads_in_5S_and_" + "_".join(sorted(used[q]))] += 1
    used.pop(q, None)
pending = {}
with pysam.AlignmentFile(BAM, "rb", check_sq=False) as f:
    for r in f.fetch(until_eof=True):
        acc["records"] += 1
        if r.is_unmapped or r.is_supplementary: continue
        nh = r.get_tag("NH")
        if not r.is_secondary:
            acc["primary_mapped"] += 1
            if nh == 1 and (r.infer_read_length() or r.query_length) >= 18:
                acc["eligible_NH1_ge18"] += 1
                ms = models(r)
                if not ms: acc["model_none"] += 1
                elif len(ms) > 1: acc["model_ambiguous"] += 1
                else:
                    m = ms.pop()
                    if mrole[m] == "competitor_only": acc["model_competitor_only"] += 1
                    elif not r.query_sequence: acc["model_no_SEQ"] += 1
                    else:
                        s = r.query_sequence.upper(); s = rc(s) if r.is_reverse else s
                        T[mclass[m]][s] += 1; acc[mclass[m] + "_reads"] += 1; used[r.query_name].add(mclass[m])
                        acc[mclass[m] + "_softclip_bases"] += sum(n for op, n in r.cigartuples if op == 4)
                        acc[mclass[m] + "_read_bases"] += len(s)
        q = r.query_name; pending.setdefault(q, []).append(r)
        if len(pending[q]) == nh: grouped(pending.pop(q))
    for rs in pending.values(): acc["reads_records_lt_NH"] += 1; grouped(rs)
for t, c in T.items(): acc[t + "_unique_sequences"] = len(c)
with gzip.open(f"{OUT}/{LIB}.k1.json.gz", "wt") as w:
    json.dump(dict(library=LIB, bam=BAM, accounting=dict(acc), types={t: dict(c) for t, c in T.items()}), w, sort_keys=True)
print(json.dumps(dict(acc), sort_keys=True))
