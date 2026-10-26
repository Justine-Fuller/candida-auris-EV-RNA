#!/usr/bin/env python3
"""fl_extra_count.py -- FL r3 (2026-10-09). READ-ONLY per-library counter for the extra full-length units.
Usage: fl_extra_count.py <bam> <library> <sRNA|polyA> <bin100 existing> <bin100 new> <bin150 existing> <bin150 new>
                         <rRNA_5S_loci_V3.bed> <outdir>
Frozen v2 read rules: primary, mapped, NH == 1, read length (incl. soft clips) >= 18 for transcript counts and bin coverage;
RNA strand = read strand (EV / Cell sRNA) or flipped (Cell polyA). Models = existing binned models (ordinary + competitor-only)
PLUS the 26 new CDS models, so ambiguity is evaluated against the full model set.
(1) Transcript physical reads, reproducing the counters the FL r2 inputs came from:
    - every read: touched models = models with >= 1 aligned base (pysam blocks) on their bin exonic positions, same strand;
      0 -> no_model; > 1 -> ambiguous_model; 1 competitor-only -> competitor_only; otherwise the read belongs to that model;
    - sRNA (bin100 Option-1 B4/C4 tables, unmerged, no 3' rule): greatest aligned overlap over the model's 100-nt bins;
      unique maximum or an exact tie of 2 NEIGHBOURING bins -> physical read; non-neighbouring 2-bin or >= 3-bin tie -> unassigned;
    - polyA (v2_bin150 R2 counter): greatest overlap over the model's 150-nt bins; unique or 2-bin tie -> physical read;
      >= 3-bin tie -> unassigned.
(2) 100-nt bin aligned bases (per-base depth input, no assignment), exactly as trna_count.py: for every eligible read, aligned
    bases of each block in every same-strand 100-nt bin (existing + new), 10-kb buckets de-duplicated per block.
(3) 5S rRNA pooled unit (13 V3 loci, one unit, like a tRNA family): all alignment records of a read are grouped (until NH records
    are seen); a read >= 18 nt with >= 1 alignment overlapping a 5S locus (same RNA strand, >= 1 aligned base) is a 5S read.
    Bases: aligned bases of ONE alignment (the primary if it hits 5S, otherwise the first that does) on its first-hit locus,
    mapped to 5S coordinates 5'->3' and summed into the unmerged 100-nt bins of the 118-nt unit: bin1 = 0-100, bin2 = 100-118.
Output: <lib>.flextra.json"""
import sys, json, csv, collections, pysam
BAM, LIB, KIND, A100, N100, A150, N150, R5, OUT = sys.argv[1:10]
assert KIND in ("sRNA", "polyA")
B = 10000; U5 = 118
def load(paths):
    idx = collections.defaultdict(list); bid, bmod, bnum = [], [], []; role = {}
    for p in paths:
        for r in csv.DictReader(open(p), delimiter="\t"):
            k = len(bid); bid.append(r["bin_id"]); bmod.append(r["model_id"]); bnum.append(int(r["bin_number"])); role[r["model_id"]] = r["role"]
            for blk in r["genomic_blocks_5to3"].split(";"):
                s, e = map(int, blk.split("-"))
                for b in range(s // B, (e - 1) // B + 1): idx[(r["chrom"], r["strand"], b)].append((s, e, k))
    return idx, bid, bmod, bnum, role
I1, BID1, BMOD1, BNUM1, ROLE = load([A100, N100])
I5, BID5, BMOD5, BNUM5, ROLE5 = load([A150, N150])
assert set(BMOD1) == set(BMOD5), "100- and 150-nt annotations must hold the same models"
newmods = {r["model_id"] for r in csv.DictReader(open(N100), delimiter="\t")}
loci = []; lidx = collections.defaultdict(list)
for l in open(R5):
    if l.startswith("#"): continue
    ch, s, e, name, ln, st, v2 = l.rstrip("\n").split("\t"); s, e = int(s), int(e)
    k = len(loci); loci.append(dict(ch=ch, s=s, e=e, st=st, name=name))
    for b in range(s // B, (e - 1) // B + 1): lidx[(ch, st, b)].append(k)
def rstrand(a):
    s = "-" if a.is_reverse else "+"
    return ({"+": "-", "-": "+"}[s]) if KIND == "polyA" else s
def overlaps(a, idx):
    ov = collections.Counter(); rs = rstrand(a)
    for s, e in a.get_blocks():
        seen = set()
        for b in range(s // B, (e - 1) // B + 1):
            for it in idx.get((a.reference_name, rs, b), ()):
                if it in seen: continue
                seen.add(it); bs, be, k = it; o = min(e, be) - max(s, bs)
                if o > 0: ov[k] += o
    return ov
acc = collections.Counter(); phys = collections.Counter(); unass = collections.Counter(); binbases = collections.Counter()
r5 = collections.Counter(); b5 = [0, 0]; perlocus = collections.Counter(); pending = {}
def five_s(rs):
    L_ = rs[0].infer_read_length() or rs[0].query_length
    if L_ < 18: return
    hits = {}
    for i, a in enumerate(rs):
        h = set()
        for s, e in a.get_blocks():
            for b in range(s // B, (e - 1) // B + 1):
                for k in lidx.get((a.reference_name, rstrand(a), b), ()):
                    if loci[k]["s"] < e and s < loci[k]["e"]: h.add(k)
        hits[i] = h
    if not any(hits.values()): return
    r5["reads"] += 1
    if rs[0].get_tag("NH") == 1: r5["reads_NH1"] += 1
    pi = next((i for i, a in enumerate(rs) if not a.is_secondary and hits[i]), next(i for i in hits if hits[i]))
    k0 = sorted(hits[pi])[0]; Lc = loci[k0]; perlocus[Lc["name"]] += 1
    for s, e in rs[pi].get_blocks():
        for g in range(max(s, Lc["s"]), min(e, Lc["e"])):
            off = g - Lc["s"] if Lc["st"] == "+" else Lc["e"] - 1 - g
            if off < U5: b5[0 if off < 100 else 1] += 1
            else: r5["bases_beyond_118"] += 1
with pysam.AlignmentFile(BAM, "rb", check_sq=False) as f:
    for r in f.fetch(until_eof=True):
        acc["records"] += 1
        if r.is_unmapped or r.is_supplementary: continue
        nh = r.get_tag("NH")
        if not r.is_secondary:
            acc["primary_mapped"] += 1
            if nh == 1 and (r.infer_read_length() or r.query_length) >= 18:
                acc["eligible_NH1_ge18"] += 1
                ov1 = overlaps(r, I1)
                for k, o in ov1.items(): binbases[k] += o
                mods = {BMOD1[k] for k in ov1}
                if not mods: acc["no_model"] += 1
                elif len(mods) > 1:
                    acc["ambiguous_model"] += 1
                    for m in mods:
                        if m in newmods: unass[(m, "ambiguous_model")] += 1
                else:
                    m = mods.pop()
                    if ROLE[m] == "competitor_only": acc["competitor_only"] += 1
                    else:
                        if KIND == "sRNA":
                            mx = max(ov1.values()); top = sorted(BNUM1[k] for k, o in ov1.items() if o == mx)
                            ok = len(top) == 1 or (len(top) == 2 and top[1] - top[0] == 1)
                        else:
                            ov5 = overlaps(r, I5); mx = max(ov5.values()); top = [k for k, o in ov5.items() if o == mx]
                            ok = len(top) <= 2
                        if ok: phys[m] += 1; acc["physical"] += 1
                        else: unass[(m, "tie_unassigned")] += 1; acc["tie_unassigned"] += 1
        q = r.query_name; pending.setdefault(q, []).append(r)
        if len(pending[q]) == nh: five_s(pending.pop(q))
    for rs in pending.values(): acc["reads_records_lt_NH"] += 1; five_s(rs)
json.dump(dict(library=LIB, kind=KIND, bam=BAM, accounting=dict(acc), physical_reads=dict(sorted(phys.items())),
               new_model_unassigned={f"{m}|{why}": n for (m, why), n in sorted(unass.items())},
               bin_aligned_bases={BID1[k]: n for k, n in sorted(binbases.items())},
               rRNA_5S=dict(reads=r5["reads"], reads_NH1=r5["reads_NH1"], bases_beyond_118=r5["bases_beyond_118"],
                            bin_aligned_bases={"5S_rRNA|bin001": b5[0], "5S_rRNA|bin002": b5[1]}, primary_locus_reads=dict(perlocus))),
          open(f"{OUT}/{LIB}.flextra.json", "w"), sort_keys=True)
print(json.dumps(dict(acc), sort_keys=True), json.dumps(dict(r5)))
