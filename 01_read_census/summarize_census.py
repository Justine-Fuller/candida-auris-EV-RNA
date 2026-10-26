#!/usr/bin/env python3
"""summarize_census.py -- gates + summary tables for the read census (2026-10-09).
Usage: summarize_census.py <census_dir> <manifest (header: library bam kind compartment group replicate)> <outdir> [r1_census_dir]
r2 (2026-10-09): adds Cell polyA (independent reference = v2_bin150/run/R2_counts); the QNAME-adjacency gate is replaced by
NH-based grouping (records == NH for every read); if r1_census_dir is given, every library present in r1 must give an
identical census (all class tables and accounting, except the r1-only adjacency counter).
Gates per library: accounting closes; records per read == NH for every read; QNAME records adjacent;
mapped reads (all lengths) == primary_mapped and NH=1 >= 18 nt reads == eligible_NH1_ge18 from the independent v2 bin
counters (EV: v2_bin100_option1/run/B2_counts; Cell sRNA: run/C2_counts_cellsRNA)."""
import sys, os, json, csv, collections, statistics
CD, MAN, OUT = sys.argv[1:4]; R1 = sys.argv[4] if len(sys.argv) > 4 else None; os.makedirs(OUT, exist_ok=True)
BR = "/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/run"
man = list(csv.DictReader(open(MAN), delimiter="\t"))
CLS = ["tRNA", "rRNA_5S", "rDNA_lncRNA", "OTHER_NC", "PCG_CDS", "PCG_UTR", "lncRNA", "StringTie_only", "antisense_only", "intergenic"]
gates, libs = [], []; ok = True
pool = collections.defaultdict(collections.Counter); lenpool = collections.defaultdict(collections.Counter)
combo = collections.defaultdict(collections.Counter); trm = collections.defaultdict(collections.Counter)
for r in man:
    L = r["library"]; comp = r["compartment"]; g = r["group"]
    o = json.load(open(f"{CD}/{L}.census.json")); a = o["accounting"]
    refp = {"EV": f"{BR}/B2_counts", "Cell_sRNA": f"{BR}/C2_counts_cellsRNA",
            "Cell_polyA": "/scratch/flljus002/Phase3D_forensic_audit/v2_bin150/run/R2_counts"}[comp]
    ref = json.load(open(f"{refp}/{L}.accounting.json"))
    checks = {"closes": a["closes"], "records_eq_NH": a.get("reads_records_ne_NH", 0) == 0,
              "mapped_eq_independent": a["mapped_reads_all_lengths"] == ref["primary_mapped"],
              "NH1_ge18_eq_independent": a["reads_ge18_NH1"] == ref["eligible_NH1_ge18"]}
    if R1 and os.path.exists(f"{R1}/{L}.census.json"):
        o1 = json.load(open(f"{R1}/{L}.census.json")); a1 = {k: v for k, v in o1["accounting"].items() if k != "qname_not_adjacent"}
        checks["identical_to_r1"] = all(o1[k] == o[k] for k in o if k not in ("accounting", "bam")) and a1 == a
    for k, v in checks.items(): gates.append((L, k, "PASS" if v else "FAIL")); ok &= v
    tot = a["reads_ge18_NH1"] + a["reads_ge18_NHgt1"]
    row = dict(library=L, compartment=comp, group=g, mapped_all_lengths=a["mapped_reads_all_lengths"], reads_ge18=tot,
               pct_NH1=round(100 * a["reads_ge18_NH1"] / tot, 2))
    for c in CLS:
        n = sum(v for k, v in o["primary_class"].items() if k.split("|")[0] == c); row[f"{c}_pct"] = round(100 * n / tot, 2)
    libs.append(row)
    for key in (comp, f"{comp}|{g}"):
        for k, v in o["primary_class"].items(): pool[key][k] += v
        for k, v in o["length_by_primary_class"].items(): lenpool[key][k] += v
        for k, v in o["multi_class_combinations"].items(): combo[key][k] += v
        for k, v in o["trna_membership"].items(): trm[key][k] += v
def W(name, hdr, rows):
    with open(f"{OUT}/{name}", "w") as f: f.write("\t".join(hdr) + "\n"); [f.write("\t".join(map(str, x)) + "\n") for x in rows]
W("gates.tsv", ["library", "gate", "result"], gates)
W("library_summary.tsv", list(libs[0]), [list(x.values()) for x in libs])
rows = []
for key in sorted(pool):
    tot = sum(pool[key].values())
    for c in CLS:
        n1, n2 = pool[key].get(f"{c}|NH1", 0), pool[key].get(f"{c}|NHgt1", 0)
        hist = {int(k.split("|")[1]): v for k, v in lenpool[key].items() if k.split("|")[0] == c}
        med = "NA"
        if hist:                                              # weighted lower median
            half, cum = (sum(hist.values()) + 1) // 2, 0
            for l in sorted(hist):
                cum += hist[l]
                if cum >= half: med = l; break
        frac_30_35 = round(100 * sum(v for l, v in hist.items() if 30 <= l <= 35) / max(1, sum(hist.values())), 1)
        rows.append([key, c, n1, n2, n1 + n2, round(100 * (n1 + n2) / tot, 2), round(100 * n1 / max(1, n1 + n2), 1), med, frac_30_35])
W("class_summary.tsv", ["pool", "class", "reads_NH1", "reads_NHgt1", "reads_total", "pct_of_reads_ge18", "pct_of_class_reads_NH1",
                        "median_read_length", "pct_reads_30_35nt"], rows)
W("multi_class_combinations.tsv", ["pool", "combination|NH", "reads"], [[k, c, n] for k in sorted(combo) for c, n in combo[k].most_common(15)])
W("trna_membership.tsv", ["pool", "NH|in_GTF_set|in_tRNAscan168_set", "reads"], [[k, c, n] for k in sorted(trm) for c, n in sorted(trm[k].items())])
W("length_histograms.tsv", ["pool", "class", "read_length(61=61+)", "reads"],
  [[k, *c.split("|"), n] for k in sorted(lenpool) for c, n in sorted(lenpool[k].items(), key=lambda x: (x[0].split("|")[0], int(x[0].split("|")[1])))])
print(open(f"{OUT}/class_summary.tsv").read()); print("CENSUS_ALL_GATES_PASS" if ok else "CENSUS_GATE_FAILURE"); sys.exit(0 if ok else 1)
