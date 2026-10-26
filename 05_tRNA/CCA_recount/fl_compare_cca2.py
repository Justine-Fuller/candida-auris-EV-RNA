#!/usr/bin/env python3
"""fl_compare_cca2.py -- CCA recount (2026-10-13): compare FL r3 (frozen) with FL r3 rebuilt on the recounted CCA reads.
Usage: fl_compare_cca2.py <FL_r3 unit_table.tsv.gz> <FL_r3_cca2_incl dir> <FL_r3_cca2_excl dir> <results_A_real dir> <out.tsv>
Rows: every unit whose full-length groups change (any compartment, TPM with and without rRNA) + a sanity count of units that
did not change. Codon list (plan Q1 b): families higher in EV (T2 option A, q_dataset < 0.05) in >= 3 of 4 groups AND
full-length in EV (TPM with rRNA, FL r3 rule) in >= 1 group -- old vs new."""
import sys, csv, gzip, collections
OLD, INC, EXC, A, OUT = sys.argv[1:6]
rd = lambda p: {(r["compartment"], r["unit"]): r for r in csv.DictReader(gzip.open(p, "rt"), delimiter="\t")}
O = rd(OLD); V = {"incl": rd(f"{INC}/unit_table.tsv.gz"), "excl": rd(f"{EXC}/unit_table.tsv.gz")}
cols = ["full_length_groups_TPM_with_rRNA", "full_length_groups_TPM_without_rRNA"]
lines = ["version\tcompartment\tunit\tclass\tcolumn\told\tnew\tCCA_old\tCCA_new"]; same = collections.Counter()
for v, D in V.items():
    assert set(D) == set(O)
    for k, r in O.items():
        ch = False
        for c in cols:
            if r[c] != D[k][c]:
                ch = True; lines.append(f"{v}\t{k[0]}\t{k[1]}\t{r['class']}\t{c}\t{r[c] or '-'}\t{D[k][c] or '-'}\t{r['CCA_reads_compartment']}\t{D[k]['CCA_reads_compartment']}")
        if not ch: same[(v, r["class"] != "tRNA_family")] += 1
G = ["B8441_untreated", "B8441_caspofungin", "MMC1_untreated", "MMC1_caspofungin"]
up = collections.Counter()
for g in G:
    for r in csv.DictReader(gzip.open(f"{A}/{g}.ev_vs_cell_bins.tsv.gz", "rt"), delimiter="\t"):
        if r["q_dataset"] not in ("NA", "") and float(r["q_dataset"]) < 0.05 and r["direction_observed"] == "higher_in_EV": up[r["bin_id"].split("|")[1]] += 1
up3 = sorted(f for f, n in up.items() if n >= 3)
fl = lambda D, f: D.get(("EV", f), {}).get("full_length_groups_TPM_with_rRNA", "") not in ("", "NA")
lines.append(f"#codon_Q1b\tEV-up>=3/4: {','.join(up3)}")
lines.append(f"#codon_Q1b\told: {','.join(f for f in up3 if fl(O, f))}")
for v, D in V.items(): lines.append(f"#codon_Q1b\tnew_{v}: {','.join(f for f in up3 if fl(D, f))}")
lines.append(f"#unchanged_units\t{dict(same)}")
open(OUT, "w").write("\n".join(lines) + "\n"); print("\n".join(lines))
