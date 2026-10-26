#!/usr/bin/env python3
"""s6f_sets_kmers.py -- Step 6F (v2, 2026-10-13), stage F4: gene sets and 6-mer counts per gene, per filter x group x region.
Justine's choice "V2" (2026-10-13 15:19 UTC): no matching (plan Q2 b). Enriched = F3 ev_enriched (q < 0.05, log2 share ratio > 1);
comparison = every other TESTED gene of the same class in that group and filter (all non-enriched genes); abundance (Cell polyA
TPM), transcript length (len_whole) and region GC are covariates in the R regression. Filters: FX (primary) and FE (sensitivity).
Regions: PCG_CDS -> UTR5, CDS, UTR3, whole; lncRNA -> whole (only genes where the region exists).
6-mers: every overlapping 6-letter window on the RNA strand (sequences are already 5'->3'); windows with a non-ACGT letter are
skipped. Counts are written per gene (long format, nonzero); rates are formed in R (count * 1000 / region length).
Usage: s6f_sets_kmers.py real <F3 dir> <gene_sequences.tsv.gz> <outdir>
       s6f_sets_kmers.py sim  <F3 dir> <gene_sequences.tsv.gz> <outdir> <seed>
         sim: FX B8441_caspofungin PCG UTR3 only; "planted" = the rarest 6-mer in that set (ties: alphabetical) written once at a
         random position into the 3' UTR of 30% of the enriched genes; "gcshift" = in every enriched 3' UTR each A/T becomes
         G/C with probability 0.1 (no planted word).
       s6f_sets_kmers.py selftest"""
import sys, csv, gzip, random, itertools, collections, os, json
K = 6; WORDS = ["".join(p) for p in itertools.product("ACGT", repeat=K)]; WI = set(WORDS)
GROUPS = ["B8441_untreated", "B8441_caspofungin", "MMC1_untreated", "MMC1_caspofungin"]
REGIONS = [("PCG_CDS", "UTR5"), ("PCG_CDS", "CDS"), ("PCG_CDS", "UTR3"), ("PCG_CDS", "whole"), ("lncRNA", "whole")]
def kmer_counts(s):
    c = collections.Counter()
    for i in range(len(s) - K + 1):
        w = s[i:i + K]
        if w in WI: c[w] += 1
    return c
def gc(s): return (s.count("G") + s.count("C")) / len(s)
def selftest():
    ok = True
    def chk(n, cond): nonlocal ok; ok &= cond; print(n, "PASS" if cond else "FAIL")
    chk("overlapping windows", kmer_counts("ACGTACGTAC") == collections.Counter({"ACGTAC": 2, "CGTACG": 1, "GTACGT": 1, "TACGTA": 1}))
    chk("homopolymer overlap", kmer_counts("AAAAAAAA")["AAAAAA"] == 3)
    chk("N windows skipped", kmer_counts("ACGTNACGTAC") == collections.Counter({"ACGTAC": 1}))
    chk("shorter than k", kmer_counts("ACGTA") == collections.Counter())
    chk("lowercase not counted (sequences are upper-case)", kmer_counts("acgtac") == collections.Counter())
    chk("4096 words", len(WORDS) == 4096 and len(set(WORDS)) == 4096)
    chk("gc", abs(gc("GGCCAT") - 4 / 6) < 1e-12)
    print("F4_SELFTEST_PASS" if ok else "F4_SELFTEST_FAIL"); sys.exit(0 if ok else 1)
def write_set(d, E, C, seqd, info):
    os.makedirs(d)
    with open(f"{d}/genes.tsv", "w") as w:
        w.write("model_id\tset\tcell_polyA_TPM\tlen_whole\tregion_length\tregion_GC\n")
        for g in sorted(E) + sorted(C):
            w.write(f"{g}\t{'enriched' if g in E else 'comparison'}\t{info[g]['tpm']}\t{info[g]['len']}\t{len(seqd[g])}\t{gc(seqd[g]):.6f}\n")
    with gzip.open(f"{d}/kmer_counts.tsv.gz", "wt") as w:
        w.write("model_id\tword\tcount\n")
        for g in sorted(E) + sorted(C):
            for word, n in sorted(kmer_counts(seqd[g]).items()): w.write(f"{g}\t{word}\t{n}\n")
    for nm, gs in (("enriched", sorted(E)), ("comparison", sorted(C))):   # STREME input, RNA alphabet (T -> U)
        with open(f"{d}/{nm}.fa", "w") as w:
            for g in gs: w.write(f">{g}\n{seqd[g].replace('T', 'U')}\n")
def main():
    mode = sys.argv[1]
    if mode == "selftest": selftest()
    F3, SEQ, OUT = sys.argv[2:5]; os.makedirs(f"{OUT}/sets")
    seqs = {r["model_id"]: r for r in csv.DictReader(gzip.open(SEQ, "rt"), delimiter="\t")}
    jobs = [(f, g, c, r) for f in ("FX", "FE") for g in GROUPS for c, r in REGIONS] if mode == "real" else [("FX", "B8441_caspofungin", "PCG_CDS", "UTR3")]
    summ = []
    for flt, g, cls, reg in jobs:
        T = list(csv.DictReader(gzip.open(f"{F3}/{g}.{flt}.gene_tests.tsv.gz", "rt"), delimiter="\t"))
        rows = [r for r in T if r["tested"] == "TRUE" and r["rna_class"] == cls and int(seqs[r["model_id"]]["len_" + reg]) > 0]
        E = {r["model_id"] for r in rows if r["ev_enriched"] == "TRUE"}; C = {r["model_id"] for r in rows} - E
        seqd = {r["model_id"]: seqs[r["model_id"]]["seq_" + reg] for r in rows}
        info = {r["model_id"]: dict(tpm=r["cell_polyA_TPM_mean"], len=seqs[r["model_id"]]["len_whole"]) for r in rows}
        variants = [("", seqd)]
        if mode == "sim":
            rng = random.Random(int(sys.argv[5]))
            tot = collections.Counter(); [tot.update(kmer_counts(seqd[x])) for x in seqd]
            rare = min(WORDS, key=lambda w: (tot[w], w))
            planted = dict(seqd); chosen = sorted(rng.sample(sorted(E), round(0.3 * len(E))))
            for x in chosen:
                s = planted[x]; p = rng.randrange(0, len(s) - K + 1); planted[x] = s[:p] + rare + s[p + K:]
            shifted = {x: ("".join(rng.choice("GC") if ch in "AT" and rng.random() < 0.1 else ch for ch in s) if x in E else s) for x, s in seqd.items()}
            json.dump(dict(planted_word=rare, planted_word_total_count_before=tot[rare], planted_genes=chosen, n_enriched=len(E),
                           n_comparison=len(C)), open(f"{OUT}/sim_design.json", "w"), indent=1)
            open(f"{OUT}/sim_planted_word.txt", "w").write(rare + "\n")
            variants = [(".planted", planted), (".gcshift", shifted)]
        for vn, sd in variants:
            tag = f"{flt}.{g}.{cls}.{reg}{vn}"; write_set(f"{OUT}/sets/{tag}", E, C, sd, info)
            summ.append([tag, flt, g, cls, reg, vn.strip(".") or "real", len(E), len(C)])
    with open(f"{OUT}/sets_summary.tsv", "w") as w:
        w.write("set\tfilter\tgroup\trna_class\tregion\tvariant\tenriched\tcomparison\n")
        for s in summ: w.write("\t".join(map(str, s)) + "\n")
    print(open(f"{OUT}/sets_summary.tsv").read())
if __name__ == "__main__": main()
