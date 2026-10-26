#!/usr/bin/env python3
"""gate_extra_counts.py -- FL r3 (2026-10-09). Independent check that fl_extra_count.py reproduces the frozen FL r2 inputs
for every EXISTING model and bin (the 26 new models have no same-strand overlap with existing binned models, so nothing
may change):
  X1 physical reads per existing ordinary model == frozen table (EV: B4 / Cell sRNA: C4 transcript_physical_reads;
     Cell polyA: v2_bin150 R2 physical_reads) -- models absent from a table count as 0;
  X2 100-nt bin aligned bases per existing bin == frozen trna_count.py T1 JSON (bin_aligned_bases);
  X3 accounting: primary_mapped and eligible_NH1_ge18 == T1 JSON accounting;
  X4 5S: bases_beyond_118 == 0; reads >= reads_NH1; bin bases <= 118 nt x reads.
Usage: gate_extra_counts.py <flextra_dir> <manifest_all36 (header)> <T1_all dir> <B4 dir> <C4 dir> <bin150 R2 dir> <new_bins100.tsv> <out.tsv> [library ...]"""
import sys, json, csv, gzip
FX, MAN, T1, B4, C4, R2, NEW, OUT = sys.argv[1:9]; only = set(sys.argv[9:])
new = {r["model_id"] for r in csv.DictReader(open(NEW), delimiter="\t")}
rows = []; ok = True
for m in csv.DictReader(open(MAN), delimiter="\t"):
    L = m["library"]
    if only and L not in only: continue
    x = json.load(open(f"{FX}/{L}.flextra.json")); t = json.load(open(f"{T1}/{L}.trna.json"))
    ref = {}
    if m["compartment"] == "Cell_polyA":
        for r in csv.DictReader(gzip.open(f"{R2}/{L}.transcript_totals.tsv.gz", "rt"), delimiter="\t"): ref[r["model_id"]] = int(r["physical_reads"])
    else:
        d = B4 if m["compartment"] == "EV" else C4
        for r in csv.DictReader(gzip.open(f"{d}/{L}.bin100_final_table.tsv.gz", "rt"), delimiter="\t"): ref.setdefault(r["transcript_id"], int(r["transcript_physical_reads"]))
    got = {k: v for k, v in x["physical_reads"].items() if k not in new}
    bad1 = [k for k in set(ref) | set(got) if ref.get(k, 0) != got.get(k, 0)]
    bb = {k: v for k, v in x["bin_aligned_bases"].items() if k.split("|")[0] not in new}
    bad2 = [k for k in set(bb) | set(t["bin_aligned_bases"]) if bb.get(k, 0) != t["bin_aligned_bases"].get(k, 0)]
    a, ta = x["accounting"], t["accounting"]
    g3 = a["primary_mapped"] == ta["primary_mapped"] and a["eligible_NH1_ge18"] == ta["eligible_NH1_ge18"]
    s = x["rRNA_5S"]; g4 = s["bases_beyond_118"] == 0 and s["reads"] >= s["reads_NH1"] and sum(s["bin_aligned_bases"].values()) <= 118 * s["reads"]
    for g, v, d in (("X1_physical_reads_existing_models", not bad1, f"{len(bad1)} mismatches {sorted(bad1)[:5]}"),
                    ("X2_bin_aligned_bases_existing_bins", not bad2, f"{len(bad2)} mismatches {sorted(bad2)[:5]}"),
                    ("X3_accounting", g3, f"{a['primary_mapped']}/{a['eligible_NH1_ge18']}"), ("X4_5S_sanity", g4, json.dumps(s["bin_aligned_bases"]))):
        rows.append((L, g, "PASS" if v else "FAIL", d)); ok &= v
with open(OUT, "w") as w:
    w.write("library\tgate\tresult\tdetail\n"); [w.write("\t".join(map(str, r)) + "\n") for r in rows]
print(open(OUT).read()); print("X_ALL_PASS" if ok else "X_GATE_FAILURE"); sys.exit(0 if ok else 1)
