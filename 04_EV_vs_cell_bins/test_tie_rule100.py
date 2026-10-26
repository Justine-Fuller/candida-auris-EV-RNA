#!/usr/bin/env python3
"""Known-answer test for apply_tie_rule100.py (B4). Usage: test_tie_rule100.py <scriptdir> <workdir> (r2)
Toy: transcript X (3 bins), Y (2 bins). Strict counts X1=2, X2=3. Ties: two reads X1/X2, one read Y1/Y2.
Expected final: X1=4, X2=5, X3=0, Y1=1, Y2=1; X physical 7 (5 strict + 2 ties), summed 9; Y physical 1, summed 2.
r2: a tie between NON-neighbouring bins X1/X3 is left unassigned (X: transcript_unassigned_nonadjacent_tie_reads = 1, not
in physical reads or bin counts). Negative control: a tie naming a bin of a DIFFERENT transcript must make gate T1 FAIL."""
import sys, os, subprocess, gzip, json, csv
S, W = sys.argv[1], sys.argv[2]
H = "library\tread\tchrom\tstart0\trna_strand\tcigar\tread_length\tmodel_id\tn_bins_touched\tn_bins_at_max\tmax_overlap\tbin_overlaps_desc"
def setup(d, ties):
    os.makedirs(f"{d}/counts", exist_ok=True)
    with open(f"{d}/ann.tsv", "w") as o:
        o.write("bin_id\trna_class\tmodel_id\tgene_id\trole\trdna_flag\tbin_number\tn_bins\ttx_start0\ttx_end\tbin_length\tfinal_bin_shorter\n")
        for m, n, L in (("X", 3, 250), ("Y", 2, 150)):
            for i in range(1, n + 1):
                a, b = (i - 1) * 100, min(i * 100, L)
                o.write(f"{m}|bin{i:03d}\tPCG_CDS\t{m}\t{m}\tordinary\t\t{i}\t{n}\t{a}\t{b}\t{b - a}\t{'yes' if b - a < 100 else ''}\n")
    with open(f"{d}/man.tsv", "w") as o:
        o.write("library\tbam\tkind\tcompartment\tgroup\treplicate\nL1\tx.bam\tEV\tEV\tB8441_untreated\t1\n")
    tx = {"X": [5, 0], "Y": [0, 0]}
    for t in ties: tx[t[0]][1] += 1
    for sfx in ("", ".c1"):
        with gzip.open(f"{d}/counts/L1.bin_counts{sfx}.tsv.gz", "wt") as o: o.write("bin_id\tassigned_reads\nX|bin001\t2\nX|bin002\t3\n")
        with gzip.open(f"{d}/counts/L1.transcript_totals{sfx}.tsv.gz", "wt") as o:
            o.write("model_id\tassigned_reads\tsummed_bin_counts\tunresolved_tie_reads\n")
            for m, (a, t) in tx.items(): o.write(f"{m}\t{a}\t{a}\t{t}\n")
        with gzip.open(f"{d}/counts/L1.ties{sfx}.tsv.gz", "wt") as o:
            o.write(H + "\n")
            for k, (m, b1, b2) in enumerate(ties):
                m2 = "Y" if b2 == 99 else m; b2 = 2 if b2 == 99 else b2          # b2 == 99: bin of the other transcript
                o.write(f"L1\tq{k}\tc1\t0\t+\t20M\t20\t{m}\t2\t2\t10\t{m}|bin{b1:03d}:10;{m2}|bin{b2:03d}:10\n")
    acc = {"eligible_NH1_ge18": 5 + len(ties) + 4, "no_binned_model": 3, "ambiguous_model": 1, "competitor_only": 0,
           "assigned": 5, "tie_unresolved_2bin": len(ties), "tie_unresolved_3plus": 0}
    json.dump(acc, open(f"{d}/counts/L1.accounting.json", "w"))
    return subprocess.run([sys.executable, f"{S}/apply_tie_rule100.py", f"{d}/ann.tsv", f"{d}/counts", f"{d}/man.tsv", f"{d}/out"], capture_output=True, text=True)
ok = True
p = setup(f"{W}/pos", [("X", 1, 2), ("X", 1, 2), ("Y", 1, 2), ("X", 1, 3)])
if p.returncode: ok = False; print("positive case FAILED to pass gates\n", p.stdout[-2000:], p.stderr[-2000:])
else:
    t = {r["bin_id"]: r for r in csv.DictReader(gzip.open(f"{W}/pos/out/L1.bin100_final_table.tsv.gz", "rt"), delimiter="\t")}
    got = {b: int(r["bin_count"]) for b, r in t.items()}
    exp = {"X|bin001": 4, "X|bin002": 5, "X|bin003": 0, "Y|bin001": 1, "Y|bin002": 1}
    tx = {r["transcript_id"]: (int(r["transcript_physical_reads"]), int(r["transcript_tie_reads"]), int(r["transcript_summed_bin_counts"]),
                               int(r["transcript_unassigned_nonadjacent_tie_reads"])) for r in t.values()}
    na = open(f"{W}/pos/out/nonadjacent_ties_unassigned.tsv").read().splitlines()
    r1 = got == exp; r2 = tx == {"X": (7, 2, 9, 1), "Y": (1, 1, 2, 0)} and len(na) == 2 and "X|bin001:10;X|bin003:10" in na[1]
    print("positive: bin counts", r1, "transcripts", r2); ok &= r1 and r2
    if not (r1 and r2): print(" got", got, tx)
p = setup(f"{W}/neg", [("X", 1, 99)])
fails = [l.split("\t")[0] for l in p.stdout.splitlines() if "\tFAIL" in l]
neg = p.returncode == 1 and "T1_ties_two_bins_same_transcript_L1" in fails
print("negative control (tie naming another transcript's bin) rejected:", neg, fails); ok &= neg
print("ALL_TIE_RULE_TESTS_PASS" if ok else "TIE_RULE_TESTS_FAILED"); sys.exit(0 if ok else 1)
