#!/usr/bin/env python3
"""apply_tie_rule100.py -- Option-1 100-nt branch, stage B4 (2026-10-08). Applies Caro's confirmed tie rule to the
validated B2 counts WITHOUT recounting BAMs. NO statistics.
TIE RULE (Caro, confirmed 2026-10-08; documented intentional exception to one-read-one-bin):
  a read whose greatest aligned overlap is shared exactly by two neighbouring 100-nt bins of the same transcript
  contributes +1 to BOTH bins. Reads touching > 1 same-strand transcript (ambiguous_model) stay unassigned.
  Consequence: per transcript, summed bin counts = physical reads + tie reads (>= physical reads).
r2 (2026-10-08, Justine/Caro decision after r1 gate T1 stopped on 7 reads): an exact tie between two NON-neighbouring bins
  of the same transcript (short reads split by STAR across an unannotated gap inside the transcript) is LEFT UNASSIGNED,
  like ambiguous_model reads: no bin count, not in the physical transcript total; listed in nonadjacent_ties_unassigned.tsv
  and counted per transcript as transcript_unassigned_nonadjacent_tie_reads.
Inputs: B2 per-library bin_counts (strict greatest-overlap reads), transcript_totals and ties files (C0 and C1).
Usage: apply_tie_rule100.py <bin100_master_annotation.tsv> <B2_counts_dir> <manifest_EV.tsv> <outdir>
Gates (exit 1 on failure):
  T1 every tie read: exactly 2 bins at the maximum, both in the annotation, both in the read's transcript; no 3+-bin ties
     (neighbouring -> +1 each; non-neighbouring -> unassigned, r2)
  T2 final counts built from C0 files == built from C1 files (bins and transcripts)
  T3 per transcript: physical = strict + ties; summed bin counts = strict + 2 x ties = physical + ties; and the sum of the
     final bin counts over the transcript's bins == summed bin counts
  T4 per library: eligible = no_binned_model + ambiguous_model + competitor_only + nonadjacent ties + physical reads; sum of
     final bin counts = assigned + 2 x neighbouring ties; ambiguous_model reads unchanged and absent from all counts
  T5 every final bin count >= its strict count; bins receiving +1 from ties == distinct bins named in the ties files
  T6 independent re-read of the written per-library tables reproduces T3 for every transcript"""
import sys, os, csv, gzip, json, collections
ANN, CD, MAN, OUT = sys.argv[1:5]
os.makedirs(OUT, exist_ok=True)
rows = lambda p: list(csv.DictReader(gzip.open(p, "rt"), delimiter="\t"))
gates = []
def gate(n, ok, d=""): gates.append((n, "PASS" if ok else "FAIL", str(d)))
ann = list(csv.DictReader(open(ANN), delimiter="\t"))
binfo = {b["bin_id"]: b for b in ann}
models = collections.OrderedDict()
for b in ann: models.setdefault(b["model_id"], b)
man = list(csv.DictReader(open(MAN), delimiter="\t"))
libs = [r["library"] for r in man]

def build(L, sfx):
    """final counts from one counter's files: strict bin counts + tie +1s; returns bins, transcripts, tie checks"""
    strict = {x["bin_id"]: int(x["assigned_reads"]) for x in rows(f"{CD}/{L}.bin_counts{sfx}.tsv.gz")}
    tt = {x["model_id"]: (int(x["assigned_reads"]), int(x["unresolved_tie_reads"])) for x in rows(f"{CD}/{L}.transcript_totals{sfx}.tsv.gz")}
    plus = collections.Counter(); bad = []; tie_tx = collections.Counter(); nonadj = []; na_tx = collections.Counter()
    for x in rows(f"{CD}/{L}.ties{sfx}.tsv.gz"):
        parts = [p.rsplit(":", 1) for p in x["bin_overlaps_desc"].split(";")]
        top = [b for b, v in parts if int(v) == int(x["max_overlap"])]
        ok = (int(x["n_bins_at_max"]) == 2 and len(top) == 2 and all(b in binfo for b in top)
              and all(binfo[b]["model_id"] == x["model_id"] for b in top))
        if not ok: bad.append(x["read"]); continue
        if abs(int(binfo[top[0]]["bin_number"]) - int(binfo[top[1]]["bin_number"])) != 1:   # r2: non-neighbouring -> unassigned
            nonadj.append(x); na_tx[x["model_id"]] += 1; continue
        for b in top: plus[b] += 1
        tie_tx[x["model_id"]] += 1
    return strict, tt, plus, tie_tx, bad, nonadj, na_tx

summary, mats, mtx = [], {}, {}; NONADJ = []
for r in man:
    L = r["library"]; strain, treat = r["group"].split("_", 1)
    acc = json.load(open(f"{CD}/{L}.accounting.json"))
    s0, t0, p0, tx0, bad0, na0, nt0 = build(L, "")
    s1, t1, p1, tx1, bad1, na1, nt1 = build(L, ".c1")
    gate(f"T1_ties_two_bins_same_transcript_{L}", not bad0 and not bad1 and acc.get("tie_unresolved_3plus", 0) == 0,
         f"bad={bad0[:3]} neighbouring_ties={sum(tx0.values())} nonadjacent_unassigned={len(na0)}")
    gate(f"T2_C0_eq_C1_final_{L}", (s0, t0, p0, tx0, na0, nt0) == (s1, t1, p1, tx1, na1, nt1))
    NONADJ.extend(na0)
    final = collections.Counter(s0); final.update(p0)
    # per-transcript physical vs summed
    sb = collections.Counter()
    for b, n in final.items(): sb[binfo[b]["model_id"]] += n
    t3 = []
    T = {}
    for m in set(t0) | set(sb):
        strict, all_ties = t0.get(m, (0, 0)); ties = all_ties - nt0[m]        # r2: neighbouring ties only
        phys, summed = strict + ties, strict + 2 * ties
        T[m] = (phys, strict, ties, summed, nt0[m])
        if ties != tx0[m] or sb[m] != summed or summed != phys + ties: t3.append(m)
    gate(f"T3_transcript_summed_eq_physical_plus_ties_{L}", not t3, f"bad={t3[:3]}")
    ties_n = acc.get("tie_unresolved_2bin", 0) + acc.get("tie_unresolved_3plus", 0) - len(na0)   # neighbouring ties
    phys_lib = acc["assigned"] + ties_n
    t4 = (acc["eligible_NH1_ge18"] == acc.get("no_binned_model", 0) + acc.get("ambiguous_model", 0) + acc.get("competitor_only", 0) + acc.get("excluded_3prime_overhang", 0) + len(na0) + phys_lib   # r3
          and sum(final.values()) == acc["assigned"] + 2 * ties_n == sum(v[3] for v in T.values())
          and sum(v[0] for v in T.values()) == phys_lib)
    gate(f"T4_library_accounting_{L}", t4, f"physical={phys_lib} summed_bins={sum(final.values())} ambiguous_unassigned={acc.get('ambiguous_model', 0)}")
    gate(f"T5_final_ge_strict_and_tie_bins_{L}", all(final[b] >= s0.get(b, 0) for b in final) and
         sum(1 for b in final if final[b] > s0.get(b, 0)) == len(p0))
    with gzip.open(f"{OUT}/{L}.bin100_final_table.tsv.gz", "wt") as o:
        o.write("library\tstrain\ttreatment\treplicate\ttranscript_id\tgene_id\trna_class\trole\trdna_flag\tbin_id\tbin_number\tn_bins"
                "\tbin_tx_start0\tbin_tx_end\tbin_length\tfinal_bin_shorter\tbin_reads_strict_greatest_overlap\tbin_tie_plus1\tbin_count"
                "\ttranscript_physical_reads\ttranscript_strict_assigned_reads\ttranscript_tie_reads\ttranscript_summed_bin_counts"
                "\ttranscript_unassigned_nonadjacent_tie_reads\n")
        for b in ann:
            m = b["model_id"]; ph, st, ti, su, nu = T.get(m, (0, 0, 0, 0, 0))
            o.write("\t".join(map(str, (L, strain, treat, r["replicate"], m, b["gene_id"], b["rna_class"], b["role"], b["rdna_flag"], b["bin_id"],
                                        b["bin_number"], b["n_bins"], b["tx_start0"], b["tx_end"], b["bin_length"], b["final_bin_shorter"],
                                        s0.get(b["bin_id"], 0), p0.get(b["bin_id"], 0), final.get(b["bin_id"], 0), ph, st, ti, su, nu))) + "\n")
    mats[L] = final; mtx[L] = T
    summary.append((L, r["group"], r["replicate"], acc["eligible_NH1_ge18"], phys_lib, acc["assigned"], ties_n, sum(final.values()),
                    sum(final.values()) - phys_lib, f"{100 * (sum(final.values()) - phys_lib) / phys_lib:.3f}", len(p0),
                    sum(1 for v in T.values() if v[2]), max(p0.values()) if p0 else 0, len(na0), acc.get("excluded_3prime_overhang", 0), acc.get("ambiguous_model", 0),
                    acc.get("no_binned_model", 0), acc.get("competitor_only", 0)))

# T6 independent re-read of the written tables
t6 = True
for L in libs:
    s = collections.Counter(); su = {}; ph = {}; ti = {}
    for x in csv.DictReader(gzip.open(f"{OUT}/{L}.bin100_final_table.tsv.gz", "rt"), delimiter="\t"):
        m = x["transcript_id"]; s[m] += int(x["bin_count"]); su[m] = int(x["transcript_summed_bin_counts"])
        ph[m] = int(x["transcript_physical_reads"]); ti[m] = int(x["transcript_tie_reads"])
        t6 &= int(x["bin_count"]) == int(x["bin_reads_strict_greatest_overlap"]) + int(x["bin_tie_plus1"])
    t6 &= all(s[m] == su[m] == ph[m] + ti[m] for m in su) and len(su) == len(models)
gate("T6_written_tables_reproduce_identities_all_libraries", t6)

W = lambda name, hdr, data: open(f"{OUT}/{name}", "w").write("\t".join(hdr) + "\n" + "".join("\t".join(map(str, x)) + "\n" for x in data))
W("tie_rule_summary.tsv", ["library", "group", "replicate", "eligible_reads", "physical_reads_in_bins", "strict_assigned_reads", "tie_reads",
                           "summed_bin_counts", "excess_summed_minus_physical", "pct_excess", "bins_receiving_tie_plus1",
                           "transcripts_with_ties", "max_tie_plus1_in_one_bin", "nonadjacent_tie_unassigned", "excluded_3prime_overhang", "ambiguous_model_unassigned", "no_binned_model",
                           "competitor_only"], summary)
with gzip.open(f"{OUT}/matrix_bin_count_final.tsv.gz", "wt") as o:
    o.write("bin_id\ttranscript_id\trna_class\tbin_number\tbin_length\t" + "\t".join(libs) + "\n")
    for b in ann:
        o.write("\t".join([b["bin_id"], b["model_id"], b["rna_class"], b["bin_number"], b["bin_length"]] + [str(mats[L].get(b["bin_id"], 0)) for L in libs]) + "\n")
with gzip.open(f"{OUT}/matrix_transcript_physical_vs_summed.tsv.gz", "wt") as o:
    o.write("transcript_id\trna_class\tn_bins\t" + "\t".join(f"{L}_physical_reads" for L in libs) + "\t" +
            "\t".join(f"{L}_tie_reads" for L in libs) + "\t" + "\t".join(f"{L}_summed_bin_counts" for L in libs) + "\n")
    for m, b in models.items():
        g = lambda L, i: str(mtx[L].get(m, (0, 0, 0, 0, 0))[i])
        o.write("\t".join([m, b["rna_class"], b["n_bins"]] + [g(L, 0) for L in libs] + [g(L, 2) for L in libs] + [g(L, 3) for L in libs]) + "\n")
TH = ["library", "read", "chrom", "start0", "rna_strand", "cigar", "read_length", "model_id", "n_bins_touched", "n_bins_at_max", "max_overlap", "bin_overlaps_desc"]
W("nonadjacent_ties_unassigned.tsv", TH, [[x[k] for k in TH] for x in NONADJ])
W("gates.tsv", ["gate", "result", "detail"], gates)
print(open(f"{OUT}/gates.tsv").read()); print(open(f"{OUT}/tie_rule_summary.tsv").read())
ok = all(g[1] == "PASS" for g in gates); print("B4_ALL_GATES_PASS" if ok else "B4_GATE_FAILURE"); sys.exit(0 if ok else 1)
