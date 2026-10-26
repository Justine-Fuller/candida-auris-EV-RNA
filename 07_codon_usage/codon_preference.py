#!/usr/bin/env python3
"""
codon_preference.py  (v2 Step 6C (2), 2026-10-12; r3 2026-10-13: chosen families given by --chosen, default = r2 six)

For every tRNA family in the 43-family catalogue: decoded codon = exact reverse complement of
the anticodon (approved Q2 = a). Its share of synonymous codon usage and its rank (1 = most used)
in human (Kazusa table, standard code; same table as v1 script 19) and in C. auris (v1 genome-wide
table c_auris_codon_usage.tsv, CTG clade: CUG = Ser).

Outputs (in --outdir)
  hand_check_6_families.tsv   anticodon -> codon check for the 6 chosen families
  codon_preference_43.tsv     one row per family
  fisher_counts.tsv           2x2 counts (chosen 6 vs the other families) for the Fisher test in R
"""
import argparse
import csv
import os

# Kazusa Homo sapiens (gbpri), per thousand; copied unchanged from v1 19_c_auris_genomewide_codon_usage.py
HUMAN_USAGE = {
    "Phe": {"UUU": 17.6, "UUC": 20.3},
    "Leu": {"UUA": 7.7, "UUG": 12.9, "CUU": 13.2, "CUC": 19.6, "CUA": 7.2, "CUG": 39.6},
    "Ser": {"UCU": 15.2, "UCC": 17.7, "UCA": 12.2, "UCG": 4.4, "AGU": 12.1, "AGC": 19.5},
    "Tyr": {"UAU": 12.2, "UAC": 15.3},
    "Cys": {"UGU": 10.6, "UGC": 12.6},
    "Trp": {"UGG": 13.2},
    "Pro": {"CCU": 17.5, "CCC": 19.8, "CCA": 16.9, "CCG": 6.9},
    "His": {"CAU": 10.9, "CAC": 15.1},
    "Gln": {"CAA": 12.3, "CAG": 34.2},
    "Arg": {"CGU": 4.5, "CGC": 10.4, "CGA": 6.2, "CGG": 11.4, "AGA": 12.2, "AGG": 12.0},
    "Ile": {"AUU": 16.0, "AUC": 20.8, "AUA": 7.5},
    "Met": {"AUG": 22.0},
    "Thr": {"ACU": 13.1, "ACC": 18.9, "ACA": 15.1, "ACG": 6.1},
    "Asn": {"AAU": 17.0, "AAC": 19.1},
    "Lys": {"AAA": 24.4, "AAG": 31.9},
    "Val": {"GUU": 11.0, "GUC": 14.5, "GUA": 7.1, "GUG": 28.1},
    "Ala": {"GCU": 18.4, "GCC": 27.7, "GCA": 15.8, "GCG": 7.4},
    "Asp": {"GAU": 21.8, "GAC": 25.1},
    "Glu": {"GAA": 29.0, "GAG": 39.6},
    "Gly": {"GGU": 10.8, "GGC": 22.2, "GGA": 16.5, "GGG": 16.5},
}
HUMAN_AA = {c: aa for aa, d in HUMAN_USAGE.items() for c in d}

# approved 6 families and their hand-written codons (plan 00_PLAN_DRAFT_codon_bias_v2_2026-10-12.md)
CHOSEN = [("Tyr/GTA", "UAC"), ("Arg/CCG", "CGG"), ("Ala/AGC", "GCU"),
          ("His/GTG", "CAC"), ("Phe/GAA", "UUC"), ("Leu/TAA", "UUA")]
# r3 (Justine 2026-10-13, Q2 a after the CCA recount): + Ala/CGC (GCG), Arg/ACG (CGU), Thr/CGT (ACG)
CHOSEN_R3 = CHOSEN + [("Ala/CGC", "GCG"), ("Arg/ACG", "CGU"), ("Thr/CGT", "ACG")]

COMP = {"A": "U", "C": "G", "G": "C", "T": "A", "U": "A"}


def revcomp_rna(anticodon_dna):
    return "".join(COMP[b] for b in reversed(anticodon_dna.upper()))


def share_rank(counts, codon):
    """share (%) of codon within counts, rank (1 = largest; ties share the smallest rank)."""
    tot = sum(counts.values())
    v = counts[codon]
    rank = 1 + sum(1 for x in counts.values() if x > v)
    top = [c for c, x in counts.items() if x == max(counts.values())]
    return 100.0 * v / tot, rank, len(counts), ",".join(sorted(top))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--catalogue", required=True)
    ap.add_argument("--cauris", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--chosen", choices=["r2", "r3"], default="r2")
    a = ap.parse_args()
    global CHOSEN
    if a.chosen == "r3": CHOSEN = CHOSEN_R3
    os.makedirs(a.outdir, exist_ok=True)

    cauris = {}
    with open(a.cauris) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            cauris.setdefault(r["aa"], {})[r["codon"]] = int(r["cauris_count"])
    cauris_aa = {c: aa for aa, d in cauris.items() for c in d}

    fams = []
    with open(a.catalogue) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            fams.append(r["ctg_family"])
    assert len(fams) == 43, len(fams)

    chosen = dict(CHOSEN)
    with open(os.path.join(a.outdir, "hand_check_6_families.tsv"), "w") as fh:
        fh.write("family\tin_catalogue\tanticodon_5to3\tcodon_hand\tcodon_revcomp_computed\tmatch\t"
                 "aa_family\taa_codon_human_code\taa_codon_CTG_code\n")
        for fam, hand in CHOSEN:
            aa, ac = fam.split("/")
            comp = revcomp_rna(ac)
            fh.write(f"{fam}\t{'yes' if fam in fams else 'NO'}\t{ac}\t{hand}\t{comp}\t"
                     f"{'yes' if comp == hand else 'NO'}\t{aa}\t{HUMAN_AA.get(comp)}\t{cauris_aa.get(comp)}\n")
            assert fam in fams and comp == hand and HUMAN_AA[comp] == aa and cauris_aa[comp] == aa, fam

    rows = []
    for fam in fams:
        aa_lab, ac = fam.split("/")
        aa = "Met" if aa_lab == "iMet" else aa_lab
        codon = revcomp_rna(ac)
        assert cauris_aa[codon] == aa, (fam, codon, cauris_aa[codon])  # CTG code
        cs, cr, cn, ctop = share_rank(cauris[aa], codon)
        if HUMAN_AA[codon] == aa:
            hs, hr, hn, htop = share_rank(HUMAN_USAGE[aa], codon)
            hpref = "yes" if hr == 1 else "no"
        else:  # CUG: Ser in C. auris, Leu in human -> no human Ser share (as v1: NA, CTG-clade specific)
            hs = hr = hn = htop = None
            hpref = "NA"
        rows.append({
            "family": fam, "aa": aa, "anticodon": ac, "codon": codon, "chosen": "yes" if fam in chosen else "no",
            "human_aa_of_codon": HUMAN_AA[codon],
            "human_share_pct": "NA" if hs is None else f"{hs:.2f}", "human_rank": "NA" if hr is None else hr,
            "human_n_synonymous": "NA" if hn is None else hn, "human_preferred_codon": "NA" if htop is None else htop,
            "decodes_human_preferred": hpref,
            "cauris_share_pct": f"{cs:.2f}", "cauris_rank": cr, "cauris_n_synonymous": cn,
            "cauris_preferred_codon": ctop, "decodes_cauris_preferred": "yes" if cr == 1 else "no",
            "single_codon_aa": "yes" if cn == 1 else "no",
        })
    cols = list(rows[0].keys())
    with open(os.path.join(a.outdir, "codon_preference_43.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    with open(os.path.join(a.outdir, "fisher_counts.tsv"), "w") as fh:
        fh.write("species\tgroup\tn_families\tn_decode_preferred\tn_not\tn_NA\n")
        for sp, key in (("human", "decodes_human_preferred"), ("cauris", "decodes_cauris_preferred")):
            k = len(CHOSEN)
            for grp, sel in ((f"chosen{k}", lambda r: r["chosen"] == "yes"),
                             (f"other{43 - k}", lambda r: r["chosen"] == "no"),
                             ("all43", lambda r: True)):
                sub = [r for r in rows if sel(r)]
                y = sum(r[key] == "yes" for r in sub)
                n = sum(r[key] == "no" for r in sub)
                na = sum(r[key] == "NA" for r in sub)
                fh.write(f"{sp}\t{grp}\t{len(sub)}\t{y}\t{n}\t{na}\n")
    print("done")


if __name__ == "__main__":
    main()
