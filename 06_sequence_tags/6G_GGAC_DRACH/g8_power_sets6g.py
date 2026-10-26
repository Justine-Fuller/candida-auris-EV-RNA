#!/usr/bin/env python3
"""g8_power_sets6g.py -- Step 6G gene level, DRACH power check (Justine 2026-10-15 07:04 UTC, option b: "first measure how big an
effect the gene test could see, then decide").
Synthetic dose series on the real FX B8441_caspofungin sets (largest group: 66 / 150 enriched genes): GGACT (= one DRACH site + GGAC)
is written d times (non-overlapping, random positions) into the region of a fraction f of the EV-enriched genes; comparison genes
untouched. Doses d = 1, 2, 3, 5 sites per planted gene; fractions f = 0.3 (the original gate) and 1.0; regions UTR3 and whole mRNA.
A gene too short for d sites gets as many as fit. Each set is then tested by g6_gene_perm6g.R exactly like the real sets (GC and
noGC); g9 tabulates DRACH q and the observed DRACH rate fold (mean rate per kb, enriched / comparison) per dose, next to the fold
seen in the real sets.
Built with g5_sets6g.py's own functions (imported, unchanged).
Usage: g8_power_sets6g.py <F3 dir> <gene_sequences.tsv.gz> <outdir> <seed>"""
import sys, os, csv, gzip, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import g5_sets6g as g5
F3, SEQ, OUT, SEED = sys.argv[1:5]; SEED = int(SEED); os.makedirs(f"{OUT}/sets")
seqs = {r["model_id"]: r for r in csv.DictReader(gzip.open(SEQ, "rt"), delimiter="\t")}
T = list(csv.DictReader(gzip.open(f"{F3}/B8441_caspofungin.FX.gene_tests.tsv.gz", "rt"), delimiter="\t"))
summ = []
for reg in ("UTR3", "whole"):
    rows = [r for r in T if r["tested"] == "TRUE" and r["rna_class"] == "PCG_CDS" and int(seqs[r["model_id"]]["len_" + reg]) > 0]
    E = {r["model_id"] for r in rows if r["ev_enriched"] == "TRUE"}; C = {r["model_id"] for r in rows} - E
    seqd = {r["model_id"]: seqs[r["model_id"]]["seq_" + reg] for r in rows}
    info = {r["model_id"]: dict(tpm=r["cell_polyA_TPM_mean"], len=seqs[r["model_id"]]["len_whole"]) for r in rows}
    for f in (0.3, 1.0):
        for d in (1, 2, 3, 5):
            rng = random.Random(f"{SEED}-{reg}-{f}-{d}")
            chosen = sorted(rng.sample(sorted(E), round(f * len(E))))
            sd = dict(seqd); placed = 0
            for x in chosen:
                s = sd[x]; slots = list(range(0, len(s) - 4, 5)); rng.shuffle(slots)   # non-overlapping 5-nt slots
                for p in sorted(slots[:d]): s = s[:p] + "GGACT" + s[p + 5:]; placed += 1
                sd[x] = s
            tag = f"FX.B8441_caspofungin.PCG_CDS.{reg}.power_f{int(f * 100)}_d{d}"
            g5.write_set(f"{OUT}/sets/{tag}", E, C, sd, info)
            summ.append([tag, "FX", "B8441_caspofungin", "PCG_CDS", reg, f"power_f{int(f * 100)}_d{d}", len(E), len(C), len(chosen), placed])
with open(f"{OUT}/sets_summary.tsv", "w") as w:
    w.write("set\tfilter\tgroup\trna_class\tregion\tvariant\tenriched\tcomparison\tplanted_genes\tsites_placed\n")
    for s in summ: w.write("\t".join(map(str, s)) + "\n")
print(open(f"{OUT}/sets_summary.tsv").read()); print("G8_DONE")
