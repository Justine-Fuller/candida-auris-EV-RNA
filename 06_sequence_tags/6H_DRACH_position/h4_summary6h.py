#!/usr/bin/env python3
"""h4_summary6h.py -- Step 6H: H1 gate (sequences == Step 6E K1 PCG_CDS reads), synthetic gates, and the before/after table.
Usage: h4_summary6h.py h1gate <manifest> <H1 dir> <6E K1 dir> <outfile>
       h4_summary6h.py sim <G2 pos (SIM)> <G2 posgc (SIM)> <G2 pos (SIM nojitter)> <outfile>   (r2: null judged with position+GC
         matching and on GC-wobble-free fake data; the r1 position-only null on GC-wobbled data is reported as INFO)
       h4_summary6h.py real <6G G2_k5drach dir> <G2 gc dir> <G2 pos dir> <G2 posgc dir> <outdir>
Rule (Q4 a): the DRACH excess is "explained by position" if, after position matching, DRACH is no longer EV-higher at q < 0.05 in
>= 3 of 4 groups; effect sizes (observed EV / Cell proportion ratio, edgeR log2FC) are reported before and after."""
import sys, os, gzip, json, csv, collections
mode = sys.argv[1]
def drach(d, G, T="PCG_CDS"):
    return next(r for r in csv.DictReader(gzip.open(f"{d}/{G}.{T}.words.tsv.gz", "rt"), delimiter="\t") if r["word"] == "DRACH")
sig = lambda r, k: r[f"q_{k}"] not in ("NA", "") and float(r[f"q_{k}"]) < 0.05
if mode == "h1gate":
    MAN, H1, K1, OF = sys.argv[2:6]; ok = True; out = []
    for r in csv.DictReader(open(MAN), delimiter="\t"):
        L = r["library"]; h = json.load(gzip.open(f"{H1}/{L}.h1.json.gz", "rt")); k = json.load(gzip.open(f"{K1}/{L}.k1.json.gz", "rt"))
        agg = collections.Counter()
        for key, c in h["reads"].items(): agg[key.split("\t")[0]] += c
        same = agg == collections.Counter(k["types"]["PCG_CDS"]); ok &= same
        a = h["accounting"]; out.append(f"{L}\t{'PASS' if same else 'FAIL'}\t{a.get('PCG_CDS_reads', 0)}\t{a.get('class_known', 0)}\t{a.get('class_NA', 0)}")
    open(OF, "w").write("library\tsequences_equal_6E_K1\tPCG_reads\twith_position\tposition_NA\n" + "\n".join(out) + "\n")
    print(open(OF).read()); print("H1_GATE_PASS" if ok else "H1_GATE_FAIL"); sys.exit(0 if ok else 1)
if mode == "sim":
    D, DPG, DNJ, OF = sys.argv[2:6]; lines = []
    c = drach(D, "SIM_posconfound"); t = drach(D, "SIM_tag")
    ncall = lambda d, G: sum(sig(r, "primary") for r in csv.DictReader(gzip.open(f"{d}/{G}.PCG_CDS.words.tsv.gz", "rt"), delimiter="\t"))
    n_pos, n_pg, n_nj = ncall(D, "SIM_null"), ncall(DPG, "SIM_null"), ncall(DNJ, "SIM_null_nojitter")
    g = [("H_SIM_confound_called_without_position_matching", sig(c, "naive") and c["direction_naive"] == "higher_in_EV", f"q_naive={c['q_naive']}"),
         ("H_SIM_confound_not_called_after_position_matching", not sig(c, "primary"), f"q_matched={c['q_primary']} {c['direction_primary']}"),
         ("H_SIM_tag_called_after_position_matching", sig(t, "primary") and t["direction_primary"] == "higher_in_EV", f"q_matched={t['q_primary']}"),
         ("H_SIM_null_0_calls_position_and_GC_matched", n_pg == 0, f"calls {n_pg}"),
         ("H_SIM_null_nojitter_0_calls_position_matched", n_nj == 0, f"calls {n_nj}"),
         ("INFO_SIM_null_GC_wobble_position_only_calls", True, f"calls {n_pos} (GC difference built into the r1 fake data)")]
    open(OF, "w").write("gate\tresult\tdetail\n" + "".join(f"{n}\t{'PASS' if o else 'FAIL'}\t{d}\n" for n, o, d in g))
    print(open(OF).read()); ok = all(o for _, o, _ in g); print("H_SIM_ALL_PASS" if ok else "H_SIM_FAIL"); sys.exit(0 if ok else 1)
G6, GC, POS, PG, OUT = sys.argv[2:7]; os.makedirs(OUT, exist_ok=False)
GROUPS = ["B8441_untreated", "B8441_caspofungin", "MMC1_untreated", "MMC1_caspofungin"]; rows = []
for G in GROUPS:
    for label, d, k in (("6G_all_reads_GC_matched", G6, "primary"), ("6G_all_reads_no_GC_control", G6, "naive"),
                        ("6H_reads_with_position_no_matching", GC, "naive"), ("6H_GC_matched", GC, "primary"),
                        ("6H_position_matched", POS, "primary"), ("6H_position_and_GC_matched", PG, "primary")):
        r = drach(d, G); e, c = float(r[f"share_EV_{k}"]), float(r[f"share_Cell_{k}"])
        rows.append(dict(group=G, version=label, q=r[f"q_{k}"], direction=r[f"direction_{k}"], log2FC_edgeR=r[f"log2FC_{k}"],
                         proportion_EV=e, proportion_Cell=c, EV_over_Cell=round(e / c, 4) if c else "NA"))
with open(f"{OUT}/DRACH_before_after_position.tsv", "w") as w:
    k = list(rows[0]); w.write("\t".join(k) + "\n"); [w.write("\t".join(str(r[x]) for x in k) + "\n") for r in rows]
rule = collections.Counter((r["version"], r["direction"]) for r in rows if r["q"] not in ("NA", "") and float(r["q"]) < 0.05)
with open(f"{OUT}/DRACH_rule_ge3of4.tsv", "w") as w:
    w.write("version\tgroups_EV_higher\tgroups_Cell_higher\tcall_ge3of4\n")
    for v in dict.fromkeys(r["version"] for r in rows):
        a, b = rule[(v, "higher_in_EV")], rule[(v, "higher_in_Cell")]
        w.write(f"{v}\t{a}\t{b}\t{'higher_in_EV' if a >= 3 else ('higher_in_Cell' if b >= 3 else 'none')}\n")
print(open(f"{OUT}/DRACH_rule_ge3of4.tsv").read())
for r in rows: print(r["group"], r["version"], r["q"], r["direction"], r["EV_over_Cell"])
print("H4_DONE")
