#!/usr/bin/env python3
"""g9_power_table6g.py -- Step 6G DRACH gene-level power check: per dose set (g8) and version (GC / noGC), DRACH and GGAC q and the
observed fold of the mean DRACH rate (per kb, EV-enriched / comparison); next to the same fold in the REAL gene sets (g7 outputs).
Usage: g9_power_table6g.py <G8 dir> <G10 perm dir GC> <G10 perm dir noGC> <G7_genes_GC dir> <G7_genes_noGC dir> <outdir>"""
import sys, os, csv, gzip
G8, PGC, PNO, R7GC, R7NO, OUT = sys.argv[1:7]; os.makedirs(OUT, exist_ok=False)
sets = list(csv.DictReader(open(f"{G8}/sets_summary.tsv"), delimiter="\t"))
rows = []
for ver, D in (("GC", PGC), ("noGC", PNO)):
    for s in sets:
        f = f"{D}/{s['set']}.words.tsv.gz"
        if not os.path.exists(f): rows.append(dict(version=ver, set=s["set"], note="missing")); continue
        W = {r["word"]: r for r in csv.DictReader(gzip.open(f, "rt"), delimiter="\t")}
        d, g = W["DRACH"], W["GGAC"]; me, mc = float(d["mean_rate_per_kb_enriched"]), float(d["mean_rate_per_kb_comparison"])
        rows.append(dict(version=ver, set=s["set"], region=s["region"], fraction_planted=s["variant"].split("_")[1][1:] + "%",
                         sites_per_planted_gene=s["variant"].split("_d")[1], planted_genes=s["planted_genes"], enriched=s["enriched"],
                         DRACH_rate_enriched=round(me, 3), DRACH_rate_comparison=round(mc, 3), DRACH_fold=round(me / mc, 3) if mc else "NA",
                         DRACH_q=d["q"], DRACH_called=float(d["q"]) < 0.05 and d["direction"] == "higher_in_enriched",
                         GGAC_q=g["q"], GGAC_called=float(g["q"]) < 0.05 and g["direction"] == "higher_in_enriched"))
real = []
for ver, D in (("GC", R7GC), ("noGC", R7NO)):
    for r in csv.DictReader(open(f"{D}/GGAC_DRACH_gene_level.tsv"), delimiter="\t"):
        if r["word"] != "DRACH": continue
        me, mc = float(r["mean_rate_per_kb_enriched"]), float(r["mean_rate_per_kb_comparison"])
        real.append(dict(version=ver, set=r["set"], DRACH_rate_enriched=round(me, 3), DRACH_rate_comparison=round(mc, 3),
                         DRACH_fold=round(me / mc, 3) if mc else "NA", DRACH_q=r["q"], direction=r["direction"]))
def wr(p, R):
    k = list(dict.fromkeys(x for r in R for x in r))
    with open(p, "w") as w: w.write("\t".join(k) + "\n"); [w.write("\t".join(str(r.get(x, "")) for x in k) + "\n") for r in R]
wr(f"{OUT}/power_dose_table.tsv", rows); wr(f"{OUT}/real_DRACH_folds.tsv", real)
for r in rows: print(r)
for r in real: print(r)
print("G9_DONE")
