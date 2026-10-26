#!/usr/bin/env python3
"""
make_null_foregrounds.py  (v2 Step 6C (3) null, 2026-10-12; r2 2026-10-12)
r2: --primary_only is passed to codon_scorer.load_gene_codon_counts (must match the real-run setting).

All alternative codon sets for the 6 approved amino acids (one synonymous codon each, standard code):
Tyr 2 x Arg 6 x Ala 4 x His 2 x Phe 2 x Leu 6 = 1,152 sets; the real set is excluded -> 1,151.
Each set is scored with codon_scorer.score_genes (same code path as the real run); foreground =
top 10%. Checks: (a) the real set, scored here, gives exactly the real-run foreground file;
(b) every set scores the same genes (background identical), so the enrichment universe is fixed.

Outputs (in --outdir)
  null_sets.tsv            set_id, codon spec, number of codons shared with the real set
  null_foregrounds.tsv     set_id <TAB> comma-separated Ensembl gene IDs (foreground, rank order)
  null_generation_QC.txt
"""
import argparse
import itertools
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import codon_scorer as cs  # noqa: E402

# r3 (2026-10-13): the real set may hold several codons per amino acid (--real "Tyr:UAC;Arg:CGG,CGU;..."); the null =
# EVERY alternative set with the same number of codons per amino acid (itertools.combinations within each amino acid),
# i.e. the exhaustive r2 null generalised. With the r2 real set this enumerates exactly the r2 1,152 sets in the same order.
# --task/--ntasks: write only every ntasks-th alternative set (for parallel generation); files are concatenated afterwards.
AA_ORDER = ["Tyr", "Arg", "Ala", "His", "Phe", "Leu", "Thr"]
SYN = {
    "Tyr": ["UAU", "UAC"],
    "Arg": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "Ala": ["GCU", "GCC", "GCA", "GCG"],
    "His": ["CAU", "CAC"],
    "Phe": ["UUU", "UUC"],
    "Leu": ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
    "Thr": ["ACU", "ACC", "ACA", "ACG"],
}


def parse(spec_str):
    return {x.split(":")[0]: tuple(sorted(x.split(":")[1].split(","), key=lambda c: SYN[x.split(":")[0]].index(c))) for x in spec_str.split(";")}


def spec(d):
    return ";".join(f"{aa}:{','.join(d[aa])}" for aa in AA_ORDER if aa in d)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--real_foreground", required=True)
    ap.add_argument("--real_background", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--primary_only", action="store_true")
    ap.add_argument("--real", default="Tyr:UAC;Arg:CGG;Ala:GCU;His:CAC;Phe:UUC;Leu:UUA")
    ap.add_argument("--task", type=int, default=1)
    ap.add_argument("--ntasks", type=int, default=1)
    a = ap.parse_args()
    REAL = parse(a.real); AAS = [aa for aa in AA_ORDER if aa in REAL]
    os.makedirs(a.outdir, exist_ok=True)
    for aa in AA_ORDER:  # synonymous lists must match the standard code exactly
        assert sorted(SYN[aa]) == sorted(c for c, x in cs.CODON_TABLE.items() if x == aa), aa

    genes, st = cs.load_gene_codon_counts(a.fasta, primary_only=a.primary_only)
    real_fg = [x.strip() for x in open(a.real_foreground) if x.strip()]
    real_bg = set(x.strip() for x in open(a.real_background) if x.strip())

    rows = cs.score_genes(genes, {aa: set(REAL[aa]) for aa in AAS})
    fg = [r["gene_id"] for r in rows[:cs.n_foreground(rows)]]
    check_real_fg = fg == real_fg
    check_real_bg = set(r["gene_id"] for r in rows) == real_bg

    all_sets = [dict(zip(AAS, combo)) for combo in itertools.product(*[list(itertools.combinations(SYN[aa], len(REAL[aa]))) for aa in AAS])]
    n_total = len(all_sets); n_exp = 1
    for aa in AAS: n_exp *= math.comb(len(SYN[aa]), len(REAL[aa]))
    alt = [d for d in all_sets if d != REAL]
    width = max(4, len(str(len(alt))))
    mine = [(i, d) for i, d in enumerate(alt, 1) if (i - 1) % a.ntasks == a.task - 1]
    suf = "" if a.ntasks == 1 else f"_part{a.task:03d}"
    n_bg_mismatch = 0
    with open(os.path.join(a.outdir, f"null_sets{suf}.tsv"), "w") as fs, \
         open(os.path.join(a.outdir, f"null_foregrounds{suf}.tsv"), "w") as ff:
        if a.task == 1: fs.write("set_id\tcodons\tn_shared_with_real\tn_foreground\tn_scored\n")
        for i, d in mine:
            sid = f"S{i:0{width}d}"
            rows = cs.score_genes(genes, {aa: set(d[aa]) for aa in AAS})
            if set(r["gene_id"] for r in rows) != real_bg:
                n_bg_mismatch += 1
            nfg = cs.n_foreground(rows)
            shared = sum(len(set(d[aa]) & set(REAL[aa])) for aa in AAS)
            fs.write(f"{sid}\t{spec(d)}\t{shared}\t{nfg}\t{len(rows)}\n")
            ff.write(sid + "\t" + ",".join(r["gene_id"] for r in rows[:nfg]) + "\n")
    ok = check_real_fg and check_real_bg and n_total == n_exp and len(alt) == n_exp - 1 and n_bg_mismatch == 0
    with open(os.path.join(a.outdir, f"null_generation_QC{suf}.txt"), "w") as fh:
        fh.write(f"real set: {spec(REAL)}\n")
        fh.write(f"primary-assembly filter: {'ON' if a.primary_only else 'OFF'}\n")
        fh.write(f"codon sets enumerated: {n_total} (expected {n_exp})\n")
        fh.write(f"alternative sets: {len(alt)} (expected {n_exp - 1}); written by this task: {len(mine)} (task {a.task}/{a.ntasks})\n")
        fh.write(f"real set re-scored here == real-run foreground file: {check_real_fg}\n")
        fh.write(f"real set re-scored here == real-run background file: {check_real_bg}\n")
        fh.write(f"alternative sets whose scored-gene set differs from the real background: {n_bg_mismatch}\n")
        fh.write(f"NULL GENERATION: {'PASS' if ok else 'FAIL'}\n")
    print(open(os.path.join(a.outdir, f"null_generation_QC{suf}.txt")).read())
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
