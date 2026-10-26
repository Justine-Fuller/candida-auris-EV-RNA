#!/usr/bin/env python3
"""trna_count.py -- v2 tRNA branch, step T1 (2026-10-09). READ-ONLY. Per-library tRNA family counts with read positions on the mature
tRNA, plus PCG/lncRNA 100-nt bin coverage (input for the later full-length pies).
Usage: trna_count.py <bam> <library> <bin100_master_annotation.tsv (unmerged)> <trna168.bed> <outdir> [sRNA|polyA]
r2 (2026-10-09): optional kind; polyA (TruSeq stranded mRNA, reverse-stranded) flips the read strand everywhere (default sRNA).
(1) PCG/lncRNA 100-nt bin coverage: for every primary, uniquely mapped (NH = 1) read >= 18 nt (frozen v2 read rules), the aligned
    bases (pysam blocks; D/N excluded) falling in each same-strand bin are summed (10-kb buckets de-duplicated per block). A bin's
    mean depth = summed aligned bases / bin length. (Reads are not assigned to one bin here: this is per-base depth, as in the
    published / v1 coverage criterion.)
(2) tRNA families (tRNAscan 168 set, 43 anticodon families): all alignment records of a read are grouped (until NH records are
    seen). A read >= 18 nt is assigned to a family when every alignment that overlaps a tRNA locus (same strand, >= 1 aligned base)
    overlaps loci of ONE family only (cross-family reads are excluded and counted). For an assigned read:
    - family_reads += 1;
    - family_bases += aligned bases inside the mature tRNA (locus span minus intron) of ONE alignment (the primary alignment if it
      hits the family, otherwise the first that does);
    - 3' end / CCA (v1 Step 6C follow-up (1) method): over the read's alignments that hit the family, if the read has no 'N' gap
      > 1000 nt and an alignment ends exactly at the locus 3' end (+ strand: reference_end == BED end; - strand: reference_start ==
      BED start), it is counted as reaching the 3' end; if a soft clip sits there, it is 'clipped'; if the clipped bases, read in RNA
      orientation, start with CCA, it is a CCA read. Each read is counted at most once per category.
(3) Positions: for each assigned read, on the chosen alignment's locus, the mature-tRNA coordinates of its first and last aligned
    base (0-based, 5'->3' along the mature tRNA with the intron removed; end exclusive) are recorded with the locus mature length,
    the CCA flag and NH = 1 flag: positions["fam|mature_len|start|end|cca|nh1"] = reads. Fragment types are defined later from these.
Output: <lib>.trna.json"""
import sys, json, collections, pysam
BAM, LIB, ANN, TB, OUT = sys.argv[1:6]; KIND = sys.argv[6] if len(sys.argv) > 6 else "sRNA"
assert KIND in ("sRNA", "polyA")
B = 10000
RC = str.maketrans("ACGTN", "TGCAN")
idx = collections.defaultdict(list); binid = []; binlen = []
with open(ANN) as f:
    h = f.readline().rstrip("\n").split("\t"); c = {k: i for i, k in enumerate(h)}
    for line in f:
        r = line.rstrip("\n").split("\t"); k = len(binid); binid.append(r[c["bin_id"]]); binlen.append(int(r[c["bin_length"]]))
        for blk in r[c["genomic_blocks_5to3"]].split(";"):
            s, e = map(int, blk.split("-"))
            for b in range(s // B, (e - 1) // B + 1): idx[(r[c["chrom"]], r[c["strand"]], b)].append((s, e, k))
tidx = collections.defaultdict(list); loci = []
with open(TB) as f:
    next(f)
    for line in f:
        ch, s, e, name, sc, st, fam, fid, intron, tier = line.rstrip("\n").split("\t")
        s, e = int(s), int(e); intr = None if intron == "NA" else tuple(map(int, intron.split("-")))
        k = len(loci); loci.append(dict(ch=ch, s=s, e=e, st=st, fam=fam, intron=intr, name=name))
        for b in range(s // B, (e - 1) // B + 1): tidx[(ch, st, b)].append(k)
def rstrand(a):
    s = "-" if a.is_reverse else "+"
    return ({"+": "-", "-": "+"}[s]) if KIND == "polyA" else s
def mature_bases(a, L):
    n = 0
    for s, e in a.get_blocks():
        lo, hi = max(s, L["s"]), min(e, L["e"])
        if lo < hi:
            n += hi - lo
            if L["intron"]:
                n -= max(0, min(hi, L["intron"][1]) - max(lo, L["intron"][0]))
    return n
def mature_span(a, L):
    """first and last aligned base of a on locus L, as mature coordinates (5'->3', intron removed); None if no overlap"""
    il = (L["intron"][1] - L["intron"][0]) if L["intron"] else 0; offs = []
    for s, e in a.get_blocks():
        for lo, hi in ((max(s, L["s"]), min(e, L["e"])),):
            if lo >= hi: continue
            pieces = [(lo, hi)]
            if L["intron"]:
                i0, i1 = L["intron"]; pieces = [(x, y) for x, y in ((lo, min(hi, i0)), (max(lo, i1), hi)) if x < y]
            for x, y in pieces:
                for g in (x, y - 1):
                    if L["st"] == "+": off = g - L["s"] - (il if L["intron"] and g >= L["intron"][1] else 0)
                    else: off = (L["e"] - 1 - g) - (il if L["intron"] and g < L["intron"][0] else 0)
                    offs.append(off)
    return (min(offs), max(offs) + 1) if offs else None
def trna_hits(a):
    hits = set()
    for s, e in a.get_blocks():
        for b in range(s // B, (e - 1) // B + 1):
            for k in tidx.get((a.reference_name, rstrand(a), b), ()):
                L = loci[k]
                if L["s"] < e and s < L["e"]: hits.add(k)
    return hits
binbases = collections.Counter(); acc = collections.Counter(); posc = collections.Counter()
fam = collections.defaultdict(collections.Counter)
pending = {}
def flush(rs):
    L_ = rs[0].infer_read_length() or rs[0].query_length
    if L_ < 18: return
    hits = {a_i: trna_hits(a) for a_i, a in enumerate(rs)}
    fams = {loci[k]["fam"] for hs in hits.values() for k in hs}
    if not fams: return
    acc["trna_reads"] += 1
    if len(fams) > 1: acc["trna_reads_cross_family_excluded"] += 1; return
    F = fams.pop(); fam[F]["reads"] += 1
    if rs[0].get_tag("NH") == 1: fam[F]["reads_NH1"] += 1
    pi = next((i for i, a in enumerate(rs) if not a.is_secondary and hits[i]), next(i for i in hits if hits[i]))
    k0 = sorted(hits[pi])[0]; fam[F]["bases"] += mature_bases(rs[pi], loci[k0])
    bigN = any(op == 3 and ln > 1000 for a in rs for op, ln in a.cigartuples)
    reach = clip = cca = False
    if not bigN:
        for i, a in enumerate(rs):
            for k in hits[i]:
                L = loci[k]; ct = a.cigartuples; seq = a.query_sequence or ""
                if not seq:                                      # secondary record without SEQ: use the primary's, re-oriented
                    p0 = next((x for x in rs if not x.is_secondary and x.query_sequence), None)
                    if p0 is not None:
                        seq = p0.query_sequence if p0.is_reverse == a.is_reverse else p0.query_sequence.upper().translate(RC)[::-1]
                if L["st"] == "+" and a.reference_end == L["e"]:
                    reach = True
                    if ct[-1][0] == 4:
                        clip = True; cl = seq[len(seq) - ct[-1][1]:]
                        if cl.upper().startswith("CCA"): cca = True
                elif L["st"] == "-" and a.reference_start == L["s"]:
                    reach = True
                    if ct[0][0] == 4:
                        clip = True; cl = seq[:ct[0][1]].upper().translate(RC)[::-1]
                        if cl.startswith("CCA"): cca = True
    else:
        fam[F]["reads_bigN_excluded_from_CCA_check"] += 1
    fam[F]["reach_3prime"] += reach; fam[F]["clipped_at_3prime"] += clip; fam[F]["CCA_reads"] += cca
    L0 = loci[k0]; ml = (L0["e"] - L0["s"]) - ((L0["intron"][1] - L0["intron"][0]) if L0["intron"] else 0); sp = mature_span(rs[pi], L0)
    if sp: posc[f"{F}|{ml}|{sp[0]}|{sp[1]}|{int(cca)}|{int(rs[0].get_tag('NH') == 1)}"] += 1
with pysam.AlignmentFile(BAM, "rb", check_sq=False) as f:
    for r in f.fetch(until_eof=True):
        acc["records"] += 1
        if r.is_unmapped or r.is_supplementary: continue
        nh = r.get_tag("NH")
        if not r.is_secondary:
            acc["primary_mapped"] += 1
            if nh == 1 and (r.infer_read_length() or r.query_length) >= 18:
                acc["eligible_NH1_ge18"] += 1; rs_ = rstrand(r)
                for s, e in r.get_blocks():
                    seen = set()
                    for b in range(s // B, (e - 1) // B + 1):
                        for it in idx.get((r.reference_name, rs_, b), ()):
                            if it in seen: continue
                            seen.add(it); bs, be, k = it; o = min(e, be) - max(s, bs)
                            if o > 0: binbases[k] += o
        q = r.query_name; pending.setdefault(q, []).append(r)
        if len(pending[q]) == nh: flush(pending.pop(q))
    for rs in pending.values(): acc["reads_records_lt_NH"] += 1; flush(rs)
json.dump(dict(library=LIB, bam=BAM, accounting=dict(acc), bin_aligned_bases={binid[k]: n for k, n in sorted(binbases.items())},
               trna_families={F: dict(v) for F, v in sorted(fam.items())}, positions=dict(sorted(posc.items()))),
          open(f"{OUT}/{LIB}.trna.json", "w"), sort_keys=True)
print(json.dumps(acc, sort_keys=True))
