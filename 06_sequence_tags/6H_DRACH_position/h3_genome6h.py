#!/usr/bin/env python3
"""h3_genome6h.py -- Step 6H Q3 a: genome-only check. DRACH density (sites per 100 nt, RNA strand) in each 100-nt CDS bin counted
from the 3' end (merged bin100 layout M1, all PCG_CDS models), pooled over genes per class; plus per gene: last bin vs the gene's
other bins (sign count and Wilcoxon-free summary: how many genes have a denser last bin).
Usage: h3_genome6h.py <bin100m merged annotation> <genome.fasta> <outdir>"""
import sys, os, csv, collections
AM, FA, OUT = sys.argv[1:4]; os.makedirs(OUT, exist_ok=False)
DRACH = {a + b + "AC" + c for a in "AGT" for b in "AG" for c in "ACT"}
G = {}; name = None; buf = []
for l in open(FA):
    if l.startswith(">"):
        if name: G[name] = "".join(buf)
        name = l[1:].split()[0]; buf = []
    else: buf.append(l.strip().upper())
G[name] = "".join(buf); RC = str.maketrans("ACGTN", "TGCAN")
tot = collections.defaultdict(lambda: [0, 0, 0]); per = collections.defaultdict(dict)
for r in csv.DictReader(open(AM), delimiter="\t"):
    if r["rna_class"] != "PCG_CDS" or r["role"] == "competitor_only": continue
    parts = []
    for blk in r["genomic_blocks_5to3"].split(";"):
        s, e = map(int, blk.split("-")); x = G[r["chrom"]][s:e]; parts.append(x.translate(RC)[::-1] if r["strand"] == "-" else x)
    s = "".join(parts); n = sum(s[i:i + 5] in DRACH for i in range(len(s) - 4)); b3 = int(r["n_bins"]) - int(r["bin_number"])
    t = tot[b3]; t[0] += n; t[1] += len(s); t[2] += 1; per[r["model_id"]][b3] = n * 100.0 / len(s)
with open(f"{OUT}/drach_density_by_bin_from_3prime_end.tsv", "w") as w:
    w.write("bins_from_3prime_end\tbins\tDRACH_sites\tnt\tDRACH_per_100nt\n")
    for b3 in sorted(tot): t = tot[b3]; w.write(f"{b3}\t{t[2]}\t{t[0]}\t{t[1]}\t{t[0] * 100 / t[1]:.4f}\n")
up = down = eq = 0
for m, d in per.items():
    if 0 not in d or len(d) < 2: continue
    rest = sum(v for k, v in d.items() if k != 0) / (len(d) - 1)
    up += d[0] > rest; down += d[0] < rest; eq += d[0] == rest
with open(f"{OUT}/last_bin_vs_rest_per_gene.tsv", "w") as w:
    w.write("genes_last_bin_denser\tgenes_last_bin_sparser\tgenes_equal\n" f"{up}\t{down}\t{eq}\n")
print(open(f"{OUT}/drach_density_by_bin_from_3prime_end.tsv").read()[:1500]); print("last bin denser / sparser / equal:", up, down, eq); print("H3_DONE")
