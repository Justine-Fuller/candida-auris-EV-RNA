#!/usr/bin/env python3
"""g3_reads_summary6g.py -- Step 6G (v2, 2026-10-14), read level: GGAC (from the 256 4-mer test) and DRACH (from the DRACH-partition
test), EV vs Cell sRNA reads, per group x RNA type, with GC control ("primary" = 6E GC-matched draw 1) and without ("naive").
GGAC is also ranked among the 256 4-mers by p (primary / naive). Rule (6G fixed method): q < 0.05, same direction, >= 3 of 4 groups.
Usage: g3_reads_summary6g.py real <G2 k4 dir> <G2 k5drach dir> <outdir>
       g3_reads_summary6g.py sim  <G2 k4 SIM dir> <G2 k5drach SIM dir> <outdir>
  sim gates: SIM_planted (6E synthetic, 3% of EV reads carry GGACTA = GGAC + DRACH site GGACT): GGAC and DRACH called higher in EV
  (primary); SIM_gcshift: 0 calls of any unit (primary); SIM_null: 0 calls (primary)."""
import sys, os, gzip, csv, collections
MODE, K4, K5, OUT = sys.argv[1:5]; os.makedirs(OUT, exist_ok=False)
GROUPS = ["B8441_untreated", "B8441_caspofungin", "MMC1_untreated", "MMC1_caspofungin"]; TYPES = ["PCG_CDS", "lncRNA", "tRNA", "5S_rRNA"]
fl = lambda x: float("nan") if x in ("NA", "") else float(x)
def rd(d, G, T):
    f = f"{d}/{G}.{T}.words.tsv.gz"
    return list(csv.DictReader(gzip.open(f, "rt"), delimiter="\t")) if os.path.exists(f) else None
def sig(r, k): return r[f"q_{k}"] not in ("NA", "") and float(r[f"q_{k}"]) < 0.05
def ncalls(rows, k): return sum(sig(r, k) for r in rows)
def rank(rows, word, k):
    ps = sorted((fl(r[f"p_{k}"]), r["word"]) for r in rows if r[f"p_{k}"] not in ("NA", ""))
    return next((i + 1 for i, (p, w) in enumerate(ps) if w == word), "NA")
fails = []
if MODE == "sim":
    lines = []
    for G, check in (("SIM_planted", "planted"), ("SIM_gcshift", "zero"), ("SIM_null", "zero")):
        a, b = rd(K4, G, "PCG_CDS"), rd(K5, G, "PCG_CDS")
        if check == "planted":
            g = next(r for r in a if r["word"] == "GGAC"); d = next(r for r in b if r["word"] == "DRACH")
            ok = sig(g, "primary") and g["direction_primary"] == "higher_in_EV" and sig(d, "primary") and d["direction_primary"] == "higher_in_EV"
            lines.append((f"G3_SIM_planted_GGAC_and_DRACH_called_EV", ok, f"GGAC q={g['q_primary']} rank {rank(a, 'GGAC', 'primary')}/256; DRACH q={d['q_primary']}"))
        else:
            n4, n5 = ncalls(a, "primary"), ncalls(b, "primary")
            lines.append((f"G3_{G}_0_calls_after_GC_matching", n4 == 0 and n5 == 0, f"k4 {n4}, DRACH-partition {n5}; naive k4 {ncalls(a, 'naive')}"))
    with open(f"{OUT}/sim_gates.tsv", "w") as w:
        w.write("gate\tresult\tdetail\n")
        for n, ok, d in lines: w.write(f"{n}\t{'PASS' if ok else 'FAIL'}\t{d}\n"); print(n, "PASS" if ok else "FAIL", d); ok or fails.append(n)
    print("G3_SIM_ALL_PASS" if not fails else f"G3_SIM_FAIL {fails}"); sys.exit(1 if fails else 0)
rows = []
for T in TYPES:
    for G in GROUPS:
        a, b = rd(K4, G, T), rd(K5, G, T)
        for unit, src in (("GGAC", a), ("DRACH", b)):
            r = next(x for x in src if x["word"] == unit)
            row = dict(type=T, group=G, unit=unit)
            for k in ("primary", "naive", "gcclass"):
                if f"q_{k}" not in r: continue
                row.update({f"q_{k}": r[f"q_{k}"], f"direction_{k}": r[f"direction_{k}"], f"log2FC_{k}": r[f"log2FC_{k}"],
                            f"proportion_EV_{k}": r[f"share_EV_{k}"], f"proportion_Cell_{k}": r[f"share_Cell_{k}"]})
            if unit == "GGAC":
                row["GGAC_rank_among_256_primary"] = rank(a, "GGAC", "primary"); row["GGAC_rank_among_256_naive"] = rank(a, "GGAC", "naive")
                row["calls_256_primary"] = ncalls(a, "primary"); row["calls_256_naive"] = ncalls(a, "naive")
            rows.append(row)
keys = list(dict.fromkeys(k for r in rows for k in r))
with open(f"{OUT}/GGAC_DRACH_read_level.tsv", "w") as w:
    w.write("\t".join(keys) + "\n"); [w.write("\t".join(str(r.get(k, "")) for k in keys) + "\n") for r in rows]
rule = []
for T in TYPES:
    for unit in ("GGAC", "DRACH"):
        for k in ("primary", "naive"):
            c = collections.Counter(r[f"direction_{k}"] for r in rows if r["type"] == T and r["unit"] == unit
                                    and r.get(f"q_{k}") not in (None, "", "NA") and float(r[f"q_{k}"]) < 0.05)
            call = next((d for d, n in c.items() if n >= 3), "none")
            rule.append((T, unit, "with_GC_control" if k == "primary" else "without_GC_control", c.get("higher_in_EV", 0), c.get("higher_in_Cell", 0), call))
with open(f"{OUT}/GGAC_DRACH_rule_ge3of4.tsv", "w") as w:
    w.write("type\tunit\tversion\tgroups_EV_higher\tgroups_Cell_higher\tcall_ge3of4\n"); [w.write("\t".join(map(str, x)) + "\n") for x in rule]
for x in rule: print(*x, sep="\t")
print("G3_DONE")
