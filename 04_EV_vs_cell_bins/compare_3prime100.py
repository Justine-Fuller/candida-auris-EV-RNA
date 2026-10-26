#!/usr/bin/env python3
"""compare_3prime100.py -- consistency of the r3 recount (3'-overhang exclusion) with the previous merged-bin counts (M2),
2026-10-09. Usage: compare_3prime100.py <M2_counts_dir> <N2_counts_dir> <manifest.tsv> <merged_annotation.tsv> <out.tsv>
Gates (exit 1 on failure), per library:
  G1 records, primary_mapped, eligible, no_binned_model, ambiguous_model, competitor_only identical (rule acts later)
  G2 library: (assigned + all ties) before - after == excluded_3prime_overhang after
  G3 per transcript: (assigned + ties) before - after == its excluded_3prime_overhang_reads
  G4 every bin's strict count after <= before (exclusion only removes reads)
Report: excluded reads (number, % of reads reaching one transcript), bins whose strict count fell (final vs other bins)."""
import sys, csv, gzip, json
A, N, MAN, ANN, OUT = sys.argv[1:6]
final = {r["bin_id"] for r in csv.DictReader(open(ANN), delimiter="\t") if int(r["bin_number"]) == int(r["n_bins"])}
rows = lambda p: list(csv.DictReader(gzip.open(p, "rt"), delimiter="\t"))
K = ("records", "primary_mapped", "eligible_NH1_ge18", "no_binned_model", "ambiguous_model", "competitor_only")
one = lambda a: a.get("assigned", 0) + a.get("tie_unresolved_2bin", 0) + a.get("tie_unresolved_3plus", 0)
res, ok = [], True
for L in [r["library"] for r in csv.DictReader(open(MAN), delimiter="\t")]:
    a0, a1 = json.load(open(f"{A}/{L}.accounting.json")), json.load(open(f"{N}/{L}.accounting.json"))
    x = a1.get("excluded_3prime_overhang", 0)
    g1 = all(a0.get(k, 0) == a1.get(k, 0) for k in K); g2 = one(a0) - one(a1) == x
    t0 = {r["model_id"]: int(r["assigned_reads"]) + int(r["unresolved_tie_reads"]) for r in rows(f"{A}/{L}.transcript_totals.tsv.gz")}
    t1 = {r["model_id"]: (int(r["assigned_reads"]) + int(r["unresolved_tie_reads"]), int(r["excluded_3prime_overhang_reads"]))
          for r in rows(f"{N}/{L}.transcript_totals.tsv.gz")}
    g3 = all(t0.get(m, 0) - t1.get(m, (0, 0))[0] == t1.get(m, (0, 0))[1] for m in set(t0) | set(t1))
    b0 = {r["bin_id"]: int(r["assigned_reads"]) for r in rows(f"{A}/{L}.bin_counts.tsv.gz")}
    b1 = {r["bin_id"]: int(r["assigned_reads"]) for r in rows(f"{N}/{L}.bin_counts.tsv.gz")}
    g4 = all(b1.get(b, 0) <= b0.get(b, 0) for b in set(b0) | set(b1))
    fell = [b for b in b0 if b1.get(b, 0) < b0[b]]
    res.append((L, g1, g2, g3, g4, x, f"{100 * x / one(a0):.2f}" if one(a0) else "NA", sum(1 for b in fell if b in final),
                sum(1 for b in fell if b not in final), sum(b0[b] - b1.get(b, 0) for b in fell if b in final)))
    ok &= g1 and g2 and g3 and g4
with open(OUT, "w") as o:
    o.write("library\tG1_upstream_accounting_identical\tG2_library_difference_eq_excluded\tG3_transcript_difference_eq_excluded"
            "\tG4_no_bin_increased\texcluded_3prime_reads\tpct_of_reads_reaching_one_transcript\tfinal_bins_decreased"
            "\tother_bins_decreased\tstrict_reads_removed_from_final_bins\n")
    for r in res: o.write("\t".join(map(str, r)) + "\n")
print(open(OUT).read()); print("3PRIME_CONSISTENCY_PASS" if ok else "3PRIME_CONSISTENCY_FAILED"); sys.exit(0 if ok else 1)
