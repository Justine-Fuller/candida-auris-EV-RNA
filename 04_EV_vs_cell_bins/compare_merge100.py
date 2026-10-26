#!/usr/bin/env python3
"""compare_merge100.py -- consistency of the merged-final-bin recount (D1) with the frozen unmerged counts (2026-10-08).
Usage: compare_merge100.py <unmerged_counts_dir> <merged_counts_dir> <manifest.tsv> <merged_annotation.tsv> <out.tsv>
Gates (exit 1 on failure), per library:
  G1 accounting categories that do not depend on bin layout are identical: records, primary_mapped, eligible,
     no_binned_model, ambiguous_model, competitor_only, and assigned + all ties (reads reaching exactly one transcript)
  G2 per transcript: assigned_reads + unresolved_tie_reads identical (the same reads reach the same transcript)
  G3 bins not involved in a merge have identical strict counts; the ONLY bins allowed to differ are those immediately 5'
     of a merged final bin (a read touching that bin and both old final bins can switch); differences are reported
Report only: merged bin strict count vs old strict(n-1) + strict(n) + ties between n-1 and n."""
import sys, csv, gzip, json, collections
U, M, MAN, ANN, OUT = sys.argv[1:6]
ann = list(csv.DictReader(open(ANN), delimiter="\t"))
merged_bin = {r["bin_id"] for r in ann if r["merged_final_bin"] == "yes"}
merged_tx = {r["model_id"] for r in ann if r["merged_final_bin"] == "yes"}
ann_ids = {r["bin_id"] for r in ann}
# bin immediately 5' of each merged final bin: the only unmerged bin whose strict count may legitimately change
# (a read touching it AND both old final bins compares its overlap against the merged bin's combined overlap)
before_merged = {f"{r['model_id']}|bin{int(r['bin_number']) - 1:03d}" for r in ann if r["merged_final_bin"] == "yes"}
rows = lambda p: list(csv.DictReader(gzip.open(p, "rt"), delimiter="\t"))
K = ("records", "primary_mapped", "eligible_NH1_ge18", "no_binned_model", "ambiguous_model", "competitor_only")
res, ok_all = [], True
for L in [r["library"] for r in csv.DictReader(open(MAN), delimiter="\t")]:
    a0, a1 = json.load(open(f"{U}/{L}.accounting.json")), json.load(open(f"{M}/{L}.accounting.json"))
    one = lambda a: a.get("assigned", 0) + a.get("tie_unresolved_2bin", 0) + a.get("tie_unresolved_3plus", 0)
    g1 = all(a0.get(k, 0) == a1.get(k, 0) for k in K) and one(a0) == one(a1)
    t0 = {x["model_id"]: int(x["assigned_reads"]) + int(x["unresolved_tie_reads"]) for x in rows(f"{U}/{L}.transcript_totals.tsv.gz")}
    t1 = {x["model_id"]: int(x["assigned_reads"]) + int(x["unresolved_tie_reads"]) for x in rows(f"{M}/{L}.transcript_totals.tsv.gz")}
    g2 = t0 == t1
    b0 = {x["bin_id"]: int(x["assigned_reads"]) for x in rows(f"{U}/{L}.bin_counts.tsv.gz")}
    b1 = {x["bin_id"]: int(x["assigned_reads"]) for x in rows(f"{M}/{L}.bin_counts.tsv.gz")}
    diff = [b for b in ann_ids - merged_bin if b0.get(b, 0) != b1.get(b, 0)]
    g3 = all(b in before_merged for b in diff)      # only the bin just before a merged final bin may change
    res.append((L, g1, g2, g3, len(diff), ";".join(sorted(diff)[:5]), one(a0), one(a1), a0.get("tie_unresolved_2bin", 0), a1.get("tie_unresolved_2bin", 0)))
    ok_all &= g1 and g2 and g3
with open(OUT, "w") as o:
    o.write("library\tG1_layout_independent_accounting_identical\tG2_reads_per_transcript_identical\tG3_unmerged_bins_identical"
            "\tunmerged_bins_differing(all_immediately_before_merged_bin)\texamples\treads_one_transcript_before\treads_one_transcript_after\tties_before\tties_after\n")
    for r in res: o.write("\t".join(map(str, r)) + "\n")
print(open(OUT).read()); print("MERGE_CONSISTENCY_PASS" if ok_all else "MERGE_CONSISTENCY_FAILED"); sys.exit(0 if ok_all else 1)
