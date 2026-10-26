#!/usr/bin/env python3
"""check_annotation100.py -- Option-1 100-nt branch (2026-10-08). Independent checks of the 100-nt bin annotation against
the frozen 150-nt B1 annotation (same models, same eligibility, same exonic positions; only the binning differs).
Usage: check_annotation100.py <bin100_outputs_dir> <frozen_B1_outputs_dir>
Coordinates 0-based half-open. Writes <bin100_outputs_dir>/../annotation_checks.tsv; exit 1 on any FAIL."""
import sys, os, csv, hashlib, collections
O100, OB1 = sys.argv[1], sys.argv[2]
BIN = 100
res = []
def gate(n, ok, d=""): res.append((n, "PASS" if ok else "FAIL", str(d)))
md5 = lambda p: hashlib.md5(open(p, "rb").read()).hexdigest()

# G8 eligibility outputs byte-identical (eligibility does not depend on bin size)
for f in ("lncRNA_eligibility_audit.tsv", "lncRNA_same_strand_eligible_pairs.tsv"):
    gate(f"G8_identical_to_B1_{f}", md5(f"{O100}/{f}") == md5(f"{OB1}/{f}"), md5(f"{O100}/{f}"))

def load(p):
    rows = list(csv.DictReader(open(p), delimiter="\t"))
    by = collections.defaultdict(list)
    for r in rows: by[r["model_id"]].append(r)
    return rows, by
r100, m100 = load(f"{O100}/bin100_master_annotation.tsv")
r150, m150 = load(f"{OB1}/bin150_master_annotation.tsv")
KEYS = ("rna_class", "model_id", "gene_id", "chrom", "strand", "role", "transcript_length", "rdna_flag", "GANNOT",
        "D4_mappability_F1_transcript", "D4_mappability_F2_transcript")
s100 = {tuple(v[0][k] for k in KEYS) for v in m100.values()}
s150 = {tuple(v[0][k] for k in KEYS) for v in m150.values()}
gate("G9_same_models_and_model_attributes_as_B1", s100 == s150, f"{len(s100)} models vs {len(s150)}")

def positions(rows):
    out = set()
    for r in rows:
        for blk in r["genomic_blocks_5to3"].split(";"):
            s, e = map(int, blk.split("-")); out.update(range(s, e))
    return out
bad10, bad11, maxnb = [], [], 0
for m, rows in m100.items():
    rows = sorted(rows, key=lambda r: int(r["bin_number"]))
    L = int(rows[0]["transcript_length"]); nb = -(-L // BIN); maxnb = max(maxnb, nb)
    ok = len(rows) == nb and all(int(r["n_bins"]) == nb for r in rows)
    for i, r in enumerate(rows, 1):
        a, b = (i - 1) * BIN, min(i * BIN, L)
        ok &= (int(r["bin_number"]) == i and int(r["tx_start0"]) == a and int(r["tx_end"]) == b and int(r["bin_length"]) == b - a
               and r["bin_id"] == f"{m}|bin{i:03d}" and (r["final_bin_shorter"] == "yes") == (i == nb and b - a < BIN))
        ok &= sum(e - s for s, e in (map(int, x.split("-")) for x in r["genomic_blocks_5to3"].split(";"))) == b - a
    if not ok: bad10.append(m)
    if positions(rows) != positions(m150[m]): bad11.append(m)
gate("G10_bins_consecutive_100nt_from_5prime_last_bin_true_length", not bad10, f"bad={bad10[:5]}")
gate("G11_same_exonic_positions_as_B1_per_model", not bad11, f"bad={bad11[:5]}")
gate("G12_bin_number_fits_3_digits", maxnb <= 999, f"max_n_bins={maxnb}")
gate("G13_bin_ids_unique", len({r["bin_id"] for r in r100}) == len(r100), len(r100))

summ = []
for cls in ("PCG_CDS", "lncRNA"):
    rows = [r for r in r100 if r["rna_class"] == cls]; ms = {r["model_id"] for r in rows}
    last = [int(r["bin_length"]) for r in rows if r["final_bin_shorter"] == "yes"]
    summ += [(f"{cls}_transcripts", len(ms)), (f"{cls}_bins", len(rows)), (f"{cls}_total_nt", sum(int(r["bin_length"]) for r in rows)),
             (f"{cls}_full_100nt_bins", sum(1 for r in rows if int(r["bin_length"]) == BIN)), (f"{cls}_shorter_final_bins", len(last)),
             (f"{cls}_shorter_final_bin_length_min_median_max", f"{min(last)}/{sorted(last)[len(last) // 2]}/{max(last)}" if last else "NA"),
             (f"{cls}_single_bin_transcripts(<=100nt)", sum(1 for m in ms if int(m100[m][0]["n_bins"]) == 1)),
             (f"{cls}_competitor_only_transcripts", len({r["model_id"] for r in rows if r["role"] == "competitor_only"})),
             (f"{cls}_bins_same_strand_overlap_other_model", sum(1 for r in rows if r["same_strand_binned_model_overlap"]))]
summ += [("all_transcripts", len(m100)), ("all_bins", len(r100))]
with open(f"{O100}/../annotation_checks.tsv", "w") as o:
    for k, v in summ: o.write(f"{k}\t{v}\n")
    for g in res: o.write("\t".join(g) + "\n")
print(open(f"{O100}/../annotation_checks.tsv").read())
sys.exit(0 if all(g[1] == "PASS" for g in res) else 1)
