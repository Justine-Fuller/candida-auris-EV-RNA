#!/usr/bin/env python3
"""h1_positions6h.py -- Step 6H (v2, 2026-10-15), stage H1: per library, the Step 6E PCG_CDS reads WITH their position class.
Plan 00_PLAN_DRAFT_DRACH_position_2026-10-15.md (Justine: all recommended). Q1 a: position = the read's 100-nt CDS bin counted from
the 3' end (0 = last bin), merged bin100 layout M1 (as Step 5A / E1); the bin containing the read's 5' end (frozen v2 P1 = 5'-end
rule). Read selection = the Step 6E PCG_CDS rule unchanged (primary, NH = 1, >= 18 nt incl. soft clips, sense, exactly one touched
non-competitor model on the unmerged bin100 master + FL r3 new_bins100). Reads whose model is not in M1 (the 26 FL r3 models) or whose
5' end lies outside the model's M1 bins get class "NA" (kept, counted, excluded from position matching).
Sequence = whole read as sequenced, RNA orientation (as 6E).
Usage: h1_positions6h.py <bam> <library> <bin100 master (unmerged)> <FL r3 new_bins100.tsv> <bin100m merged annotation> <outdir>
Output: <lib>.h1.json.gz = {"library", "accounting", "reads": {"<seq>\\t<class>": count}}"""
import sys, json, gzip, csv, collections, pysam
BAM, LIB, A100, N100, AM, OUT = sys.argv[1:7]
B = 10000; RC = str.maketrans("ACGTNacgtn", "TGCANtgcan"); rc = lambda s: s.translate(RC)[::-1]
midx = collections.defaultdict(list); bmod = []; mclass = {}; mrole = {}
for p in (A100, N100):
    for r in csv.DictReader(open(p), delimiter="\t"):
        k = len(bmod); bmod.append(r["model_id"]); mclass[r["model_id"]] = r["rna_class"]; mrole[r["model_id"]] = r["role"]
        for blk in r["genomic_blocks_5to3"].split(";"):
            s, e = map(int, blk.split("-"))
            for b in range(s // B, (e - 1) // B + 1): midx[(r["chrom"], r["strand"], b)].append((s, e, k))
mbins = collections.defaultdict(list)   # model -> [(start, end, bins_from_3prime_end)]
for r in csv.DictReader(open(AM), delimiter="\t"):
    if r["rna_class"] != "PCG_CDS": continue
    b3 = int(r["n_bins"]) - int(r["bin_number"])
    for blk in r["genomic_blocks_5to3"].split(";"):
        s, e = map(int, blk.split("-")); mbins[r["model_id"]].append((s, e, b3))
def models(a):
    hit = set(); st = "-" if a.is_reverse else "+"
    for s, e in a.get_blocks():
        for b in range(s // B, (e - 1) // B + 1):
            for bs, be, k in midx.get((a.reference_name, st, b), ()):
                if min(e, be) > max(s, bs): hit.add(bmod[k])
    return hit
def posclass(a, m):
    if m not in mbins: return "NA"
    x = (a.reference_end - 1) if a.is_reverse else a.reference_start          # genomic coordinate of the read's 5' end
    for s, e, b3 in mbins[m]:
        if s <= x < e: return str(b3)
    return "NA"
acc = collections.Counter(); R = collections.Counter()
with pysam.AlignmentFile(BAM, "rb", check_sq=False) as f:
    for r in f.fetch(until_eof=True):
        if r.is_unmapped or r.is_supplementary or r.is_secondary: continue
        if r.get_tag("NH") != 1 or (r.infer_read_length() or r.query_length) < 18: continue
        ms = models(r)
        if len(ms) != 1: continue
        m = ms.pop()
        if mrole[m] == "competitor_only" or mclass[m] != "PCG_CDS" or not r.query_sequence: continue
        s = r.query_sequence.upper(); s = rc(s) if r.is_reverse else s
        c = posclass(r, m); R[f"{s}\t{c}"] += 1; acc["PCG_CDS_reads"] += 1; acc["class_NA" if c == "NA" else "class_known"] += 1
        acc["model_not_in_M1" if m not in mbins else "model_in_M1"] += 1
with gzip.open(f"{OUT}/{LIB}.h1.json.gz", "wt") as w:
    json.dump(dict(library=LIB, accounting=dict(acc), reads=dict(R)), w, sort_keys=True)
print(json.dumps(dict(acc), sort_keys=True))
