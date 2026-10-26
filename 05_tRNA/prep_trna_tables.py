#!/usr/bin/env python3
"""prep_trna_tables.py -- v2 tRNA branch, step T2 (2026-10-09). Builds per-library count tables, in the format read by the
validated ev_vs_cell100.R (step E), for two EV-vs-Cell sRNA comparisons chosen by Justine (options A and C):
  A  family share: "transcript" = the pooled tRNA reads of a library (tRNA_pool), "bins" = the 43 anticodon families;
     bin_count = family reads; transcript_physical_reads = all family-assigned tRNA reads of that library.
  C  position within a family: "transcript" = an anticodon family, "bins" = consecutive 20-nt windows from the mature 5' end
     (window w = [20w, 20w + 20), the last window ending at the family's longest mature length). Each read (mature start/end from
     step T1) is counted once, in the window with the greatest overlap; an exact tie between two NEIGHBOURING windows adds +1 to
     both (Caro's tie rule, as for the 100-nt bins); ties across >= 3 windows are left unassigned and counted.
     transcript_physical_reads = reads placed in the family's windows (tie reads once).
Usage: prep_trna_tables.py <T1_counts_dir> <manifest (header)> <outdir>
Writes <outdir>/A/{ann.tsv,<lib>.bin100_final_table.tsv.gz} and <outdir>/C/{...}, plus <outdir>/prep_accounting.tsv and gates."""
import sys, os, json, csv, gzip, collections
CD, MAN, OUT = sys.argv[1:4]
man = list(csv.DictReader(open(MAN), delimiter="\t"))
J = {r["library"]: json.load(open(f"{CD}/{r['library']}.trna.json")) for r in man}
fams = sorted({F for o in J.values() for F in o["trna_families"]})
maxlen = collections.Counter()
for o in J.values():
    for k in o["positions"]:
        F, ml = k.split("|")[:2]; maxlen[F] = max(maxlen[F], int(ml))
W = {F: [(20 * w, min(20 * w + 20, maxlen[F])) for w in range((maxlen[F] + 19) // 20)] for F in fams if maxlen[F]}
HDR = ["transcript_id", "gene_id", "rna_class", "role", "bin_id", "bin_number", "n_bins", "bin_tx_start0", "bin_tx_end", "bin_length", "bin_count", "transcript_physical_reads"]
for d in ("A", "C"): os.makedirs(f"{OUT}/{d}", exist_ok=True)
with open(f"{OUT}/A/ann.tsv", "w") as o:
    o.write("bin_id\trole\tsame_strand_binned_model_overlap\n"); [o.write(f"tRNA_pool|{F}\tordinary\t\n") for F in fams]
with open(f"{OUT}/C/ann.tsv", "w") as o:
    o.write("bin_id\trole\tsame_strand_binned_model_overlap\n")
    for F in W: [o.write(f"{F}|w{i + 1:02d}\tordinary\t\n") for i in range(len(W[F]))]
acc_rows = []; gates = []; ok = True
for r in man:
    L = r["library"]; o = J[L]; tf = o["trna_families"]
    pool = sum(tf.get(F, {}).get("reads", 0) for F in fams)
    with gzip.open(f"{OUT}/A/{L}.bin100_final_table.tsv.gz", "wt") as f:
        f.write("\t".join(HDR) + "\n")
        for i, F in enumerate(fams, 1):
            f.write("\t".join(map(str, ["tRNA_pool", "tRNA_pool", "tRNA_family", "ordinary", f"tRNA_pool|{F}", i, len(fams), "NA", "NA", "NA", tf.get(F, {}).get("reads", 0), pool])) + "\n")
    wc = collections.defaultdict(collections.Counter); phys = collections.Counter(); ties2 = collections.Counter(); ties3 = 0
    for k, n in o["positions"].items():
        F, ml, s, e = k.split("|")[:4]; s, e = int(s), int(e)
        ov = [min(e, y) - max(s, x) for x, y in W[F]]; m = max(ov); top = [i for i, v in enumerate(ov) if v == m and v > 0]
        if len(top) == 1: wc[F][top[0]] += n; phys[F] += n
        elif len(top) == 2 and top[1] - top[0] == 1: wc[F][top[0]] += n; wc[F][top[1]] += n; phys[F] += n; ties2[F] += n
        else: ties3 += n
    with gzip.open(f"{OUT}/C/{L}.bin100_final_table.tsv.gz", "wt") as f:
        f.write("\t".join(HDR) + "\n")
        for F in W:
            for i, (x, y) in enumerate(W[F]):
                f.write("\t".join(map(str, [F, F, "tRNA_family", "ordinary", f"{F}|w{i + 1:02d}", i + 1, len(W[F]), x, y, y - x, wc[F][i], phys[F]])) + "\n")
    npos = sum(o["positions"].values())
    g1 = pool == sum(v["reads"] for v in tf.values())                                     # A: pool = all assigned reads
    g2 = all(sum(wc[F].values()) == phys[F] + ties2[F] for F in W)                       # C: window sum = physical + ties
    g3 = sum(phys.values()) + ties3 == npos                                              # C: every positioned read accounted for
    for gname, v in (("P1_pool_eq_assigned", g1), ("P2_windows_eq_physical_plus_ties", g2), ("P3_positioned_reads_accounted", g3)):
        gates.append((L, gname, "PASS" if v else "FAIL")); ok &= v
    acc_rows.append([L, r["compartment"], r["group"], pool, npos, sum(phys.values()), sum(ties2.values()), ties3,
                     round(100 * sum(ties2.values()) / max(1, sum(phys.values())), 2)])
def Wt(name, hdr, rows):
    with open(f"{OUT}/{name}", "w") as f: f.write("\t".join(hdr) + "\n"); [f.write("\t".join(map(str, x)) + "\n") for x in rows]
Wt("prep_accounting.tsv", ["library", "compartment", "group", "A_pool_reads", "C_positioned_reads", "C_physical_reads_in_windows", "C_neighbouring_tie_reads",
                           "C_ties_3plus_unassigned", "C_pct_tie_reads"], acc_rows)
Wt("prep_gates.tsv", ["library", "gate", "result"], gates)
Wt("windows.tsv", ["family", "longest_mature_length", "n_windows", "windows"], [[F, maxlen[F], len(W[F]), ";".join(f"{x}-{y}" for x, y in W[F])] for F in W])
print(open(f"{OUT}/prep_accounting.tsv").read()); print("PREP_ALL_GATES_PASS" if ok else "PREP_GATE_FAILURE"); sys.exit(0 if ok else 1)
