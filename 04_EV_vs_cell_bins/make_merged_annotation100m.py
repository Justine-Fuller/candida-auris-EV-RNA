#!/usr/bin/env python3
"""make_merged_annotation100m.py -- Option-1 100-nt branch, decision D1 (Justine, 2026-10-08): the final shorter bin of
each transcript is MERGED into the preceding bin (final bin 100-199 nt). Derived from the frozen 100-nt annotation
(FROZEN_BIN100_ANNOT.md5), which is not modified. Coordinates 0-based half-open.
Usage: make_merged_annotation100m.py <bin100_master_annotation.tsv> <outdir>
Merged bin: keeps the bin_id / bin_number of the preceding bin; tx_end = transcript length; genomic blocks = both bins'
blocks in 5'->3' order (contiguous blocks joined); flag columns = union of both bins' values; new column
merged_final_bin = "yes". n_bins of every bin of a merged transcript decreases by 1. Every other row is copied unchanged.
Gates (exit 1 on failure): M1 same transcripts; M2 per transcript the same genomic positions, each exactly once;
M3 all non-final bins exactly 100 nt and consecutive from the 5' end, final bin 100-199 nt when merged, lengths sum to
the transcript length; M4 rows of bins not involved in a merge are byte-identical to the frozen annotation (apart from
n_bins on merged transcripts and the new column); M5 no single-bin transcript created; M6 bin ids unique."""
import sys, os, csv, collections
SRC, OUT = sys.argv[1:3]
os.makedirs(OUT, exist_ok=True)
rows = list(csv.DictReader(open(SRC), delimiter="\t"))
cols = list(rows[0].keys())
FLAGS = ("rdna_flag", "GANNOT", "same_strand_binned_model_overlap", "HC_tRNA_same_strand", "HC_tRNA_opposite_strand",
         "StringTie_overlap_either_strand")
by = collections.OrderedDict()
for r in rows: by.setdefault(r["model_id"], []).append(r)
def blocks(r): return [tuple(map(int, b.split("-"))) for b in r["genomic_blocks_5to3"].split(";")]
def join(bl, strand):
    out = [list(bl[0])]
    for s, e in bl[1:]:
        if strand == "+" and out[-1][1] == s: out[-1][1] = e
        elif strand == "-" and out[-1][0] == e: out[-1][0] = s
        else: out.append([s, e])
    return ";".join(f"{s}-{e}" for s, e in out)
def union(a, b):
    v = sorted({x for x in (a.split(",") + b.split(",")) if x})
    return ",".join(v)
new, gates, merged_models = [], [], 0
for m, rs in by.items():
    rs = sorted(rs, key=lambda r: int(r["bin_number"]))
    if rs[-1]["final_bin_shorter"] == "yes" and len(rs) >= 2:
        merged_models += 1
        a, b = rs[-2], rs[-1]
        mb = dict(a)
        mb["tx_end"] = b["tx_end"]; mb["bin_length"] = str(int(a["bin_length"]) + int(b["bin_length"]))
        mb["final_bin_shorter"] = ""
        mb["genomic_blocks_5to3"] = join(blocks(a) + blocks(b), a["strand"])
        for f in FLAGS:
            if f in mb: mb[f] = union(a[f], b[f]) if f != "rdna_flag" and f != "GANNOT" else (a[f] or b[f])
        mb["merged_final_bin"] = "yes"
        keep = [dict(r, merged_final_bin="") for r in rs[:-2]] + [mb]
        for r in keep: r["n_bins"] = str(len(rs) - 1)
        new += keep
    else:
        new += [dict(r, merged_final_bin="") for r in rs]
nb = collections.OrderedDict()
for r in new: nb.setdefault(r["model_id"], []).append(r)
def pos(rs):
    out = []
    for r in rs:
        for s, e in blocks(r): out.extend(range(s, e))
    return out
gates.append(("M1_same_transcripts", list(nb) == list(by), len(nb)))
bad2 = [m for m in by if sorted(pos(by[m])) != sorted(pos(nb[m])) or len(set(pos(nb[m]))) != len(pos(nb[m]))]
gates.append(("M2_same_positions_each_once", not bad2, bad2[:3]))
bad3 = []
for m, rs in nb.items():
    L = int(rs[0]["transcript_length"]); n = len(rs)
    ok = sum(int(r["bin_length"]) for r in rs) == L and all(int(r["n_bins"]) == n for r in rs)
    for i, r in enumerate(rs, 1):
        ok &= int(r["bin_number"]) == i and int(r["tx_start0"]) == (i - 1) * 100
        ok &= int(r["bin_length"]) == 100 if i < n else (100 <= int(r["bin_length"]) <= 199 if r["merged_final_bin"] else int(r["bin_length"]) == 100)
        ok &= sum(e - s for s, e in blocks(r)) == int(r["bin_length"]) == int(r["tx_end"]) - int(r["tx_start0"])
    if not ok: bad3.append(m)
gates.append(("M3_bin_layout", not bad3, bad3[:3]))
old = {r["bin_id"]: r for r in rows}
bad4 = [r["bin_id"] for r in new if not r["merged_final_bin"] and any(r[c] != old[r["bin_id"]][c] for c in cols if c != "n_bins")]
gates.append(("M4_unmerged_rows_unchanged", not bad4, bad4[:3]))
gates.append(("M5_no_single_bin_transcripts", all(len(v) >= 2 for v in nb.values()), sum(1 for v in nb.values() if len(v) < 2)))
gates.append(("M6_bin_ids_unique", len({r["bin_id"] for r in new}) == len(new), len(new)))
with open(f"{OUT}/bin100m_master_annotation.tsv", "w") as o:
    w = csv.DictWriter(o, fieldnames=cols + ["merged_final_bin"], delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(new)
summ = [("transcripts", len(nb)), ("transcripts_merged", merged_models), ("transcripts_unchanged", len(nb) - merged_models),
        ("bins_before", len(rows)), ("bins_after", len(new))]
for cls in ("PCG_CDS", "lncRNA"):
    fin = [int(r["bin_length"]) for r in new if r["rna_class"] == cls and r["merged_final_bin"]]
    summ += [(f"{cls}_bins", sum(1 for r in new if r["rna_class"] == cls)), (f"{cls}_merged_final_bin_length_min/median/max",
             f"{min(fin)}/{sorted(fin)[len(fin) // 2]}/{max(fin)}")]
summ.append(("lncRNA_bins_same_strand_overlap_after_merge", sum(1 for r in new if r["same_strand_binned_model_overlap"])))
summ.append(("competitor_only_bins_after_merge", sum(1 for r in new if r["role"] == "competitor_only")))
with open(f"{OUT}/merge_checks.tsv", "w") as o:
    for k, v in summ: o.write(f"{k}\t{v}\n")
    for g, ok, d in gates: o.write(f"{g}\t{'PASS' if ok else 'FAIL'}\t{d}\n")
print(open(f"{OUT}/merge_checks.tsv").read())
sys.exit(0 if all(g[1] for g in gates) else 1)
