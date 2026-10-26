#!/usr/bin/env python3
"""build_tables100.py -- Option-1 100-nt branch, validation gates + per-replicate tables (2026-10-08). NO statistics.
Usage: build_tables100.py <bin100_master_annotation.tsv> <counts_dir> <manifest.tsv (header)> <outdir> [kind EV|sRNA]
r2 (2026-10-08): expected library kind is a parameter (default EV) so Cell sRNA can be tabulated; nothing else changed.
Gates (exit 1 on any failure):
  V1 C0 == C1 byte-identical for bin_counts, transcript_totals, ties, multibin_examples, overlap_patterns; accounting equal
  V2 accounting closes: eligible = no_binned_model + ambiguous_model + competitor_only + ties + assigned
  V3 every counted bin exists in the annotation; per transcript, sum of its bin counts (mapped through the annotation)
     == its assigned_reads; library total of bin counts == accounting 'assigned' (each assigned read counted once)
  V4 tie reads: rows in ties file == tie_unresolved_2bin + tie_unresolved_3plus == sum of transcript unresolved_tie_reads;
     no read name appears twice in a library's ties file; no tie read is in the assigned totals (by construction, checked
     as assigned + ties per transcript == assigned_reads + unresolved_tie_reads columns)
Outputs: per-library tables <lib>.bin100_table.tsv.gz (every bin of every eligible transcript, incl. zero-count bins);
  matrix_bin_assigned_reads.tsv.gz; matrix_transcript_totals.tsv.gz; accounting_summary.tsv; tie_summary.tsv;
  tie_examples.tsv; assignment_examples.tsv; overlap_pattern_summary.tsv; gates.tsv."""
import sys, os, csv, gzip, json, collections
ANN, CD, MAN, OUT = sys.argv[1:5]
KIND = sys.argv[5] if len(sys.argv) > 5 else "EV"   # r2: expected manifest kind (EV or sRNA)
os.makedirs(OUT, exist_ok=True)
rd = lambda p: gzip.open(p, "rt").read()
rows = lambda p: list(csv.DictReader(gzip.open(p, "rt"), delimiter="\t"))
gates = []
def gate(n, ok, d=""): gates.append((n, "PASS" if ok else "FAIL", str(d)))

ann = list(csv.DictReader(open(ANN), delimiter="\t"))
bin2model = {b["bin_id"]: b["model_id"] for b in ann}
models = collections.OrderedDict()
for b in ann: models.setdefault(b["model_id"], b)
man = list(csv.DictReader(open(MAN), delimiter="\t"))
libs = [r["library"] for r in man]
gate(f"V0_manifest_{KIND}_only_12_libraries", len(libs) == 12 and all(r["kind"] == KIND for r in man), libs)

FILES = ("bin_counts", "transcript_totals", "ties", "multibin_examples", "overlap_patterns") + (("excluded_3prime",) if os.path.exists(f"{CD}/{man[0]['library']}.excluded_3prime.tsv.gz") else ())   # r3
AKEYS = ("records", "primary_mapped", "eligible_NH1_ge18", "no_binned_model", "ambiguous_model", "competitor_only", "excluded_3prime_overhang",
         "tie_unresolved_2bin", "tie_unresolved_3plus", "assigned", "assigned_touching_2plus_bins")
acc_rows, tie_rows, tie_ex, asg_ex = [], [], [], []
patt = collections.Counter()
MB = {}; MT = {}
for r in man:
    L = r["library"]; strain, treat = r["group"].split("_", 1)
    a0 = json.load(open(f"{CD}/{L}.accounting.json")); a1 = json.load(open(f"{CD}/{L}.accounting.c1.json"))
    v1 = all(rd(f"{CD}/{L}.{n}.tsv.gz") == rd(f"{CD}/{L}.{n}.c1.tsv.gz") for n in FILES) and all(a0.get(k, 0) == a1.get(k, 0) for k in AKEYS)
    gate(f"V1_C0_eq_C1_{L}", v1)
    ties_n = a0.get("tie_unresolved_2bin", 0) + a0.get("tie_unresolved_3plus", 0)
    v2 = a0["eligible_NH1_ge18"] == sum(a0.get(k, 0) for k in ("no_binned_model", "ambiguous_model", "competitor_only", "excluded_3prime_overhang", "assigned")) + ties_n   # r3
    gate(f"V2_accounting_closes_{L}", v2 and a0["closes"])
    bc = {x["bin_id"]: int(x["assigned_reads"]) for x in rows(f"{CD}/{L}.bin_counts.tsv.gz")}
    tt = {x["model_id"]: (int(x["assigned_reads"]), int(x["unresolved_tie_reads"])) for x in rows(f"{CD}/{L}.transcript_totals.tsv.gz")}
    unknown = [b for b in bc if b not in bin2model]
    per = collections.Counter()
    for b, n in bc.items():
        if b in bin2model: per[bin2model[b]] += n
    mism = [m for m in set(per) | set(tt) if per[m] != tt.get(m, (0, 0))[0]]
    gate(f"V3_bins_sum_to_transcript_assigned_{L}", not unknown and not mism and sum(bc.values()) == a0["assigned"] == sum(v[0] for v in tt.values()),
         f"unknown={unknown[:3]} mismatched={mism[:3]} bin_sum={sum(bc.values())} assigned={a0['assigned']}")
    ti = rows(f"{CD}/{L}.ties.tsv.gz")
    names = [x["read"] for x in ti]
    gate(f"V4_ties_accounted_{L}", len(ti) == ties_n == sum(v[1] for v in tt.values()) and len(set(names)) == len(names),
         f"ties_file={len(ti)} accounting={ties_n} transcripts={sum(v[1] for v in tt.values())}")
    MB[L] = bc; MT[L] = tt
    # per-replicate full table
    with gzip.open(f"{OUT}/{L}.bin100_table.tsv.gz", "wt") as o:
        o.write("library\tstrain\ttreatment\treplicate\ttranscript_id\tgene_id\trna_class\trole\trdna_flag\tbin_id\tbin_number\tn_bins"
                "\tbin_tx_start0\tbin_tx_end\tbin_length\tfinal_bin_shorter\tbin_assigned_reads\ttranscript_assigned_reads"
                "\ttranscript_unresolved_tie_reads\n")
        for b in ann:
            m = b["model_id"]; ta, tu = tt.get(m, (0, 0))
            o.write("\t".join(map(str, (L, strain, treat, r["replicate"], m, b["gene_id"], b["rna_class"], b["role"], b["rdna_flag"],
                                        b["bin_id"], b["bin_number"], b["n_bins"], b["tx_start0"], b["tx_end"], b["bin_length"],
                                        b["final_bin_shorter"], bc.get(b["bin_id"], 0), ta, tu))) + "\n")
    # tie summaries
    cls = collections.Counter(models[x["model_id"]]["rna_class"] for x in ti)
    cig = collections.Counter(("N" in x["cigar"], "D" in x["cigar"], "S" in x["cigar"], "I" in x["cigar"]) for x in ti)
    rl = sorted(int(x["read_length"]) for x in ti); mo = sorted(int(x["max_overlap"]) for x in ti)
    denom = a0["assigned"] + ties_n
    tie_rows.append((L, r["group"], r["replicate"], a0["assigned"], ties_n, a0.get("tie_unresolved_2bin", 0), a0.get("tie_unresolved_3plus", 0),
                     f"{100 * ties_n / denom:.4f}" if denom else "NA", cls["PCG_CDS"], cls["lncRNA"],
                     sum(v for k, v in cig.items() if k[0]), sum(v for k, v in cig.items() if k[1]), sum(v for k, v in cig.items() if k[2]),
                     sum(v for k, v in cig.items() if k[3]),
                     f"{rl[0]}/{rl[len(rl) // 2]}/{rl[-1]}" if rl else "NA", f"{mo[0]}/{mo[len(mo) // 2]}/{mo[-1]}" if mo else "NA",
                     sum(1 for v in tt.values() if v[1])))
    tie_ex += ti[:5]
    ex = rows(f"{CD}/{L}.multibin_examples.tsv.gz"); asg_ex += ex[:3]
    for x in rows(f"{CD}/{L}.overlap_patterns.tsv.gz"): patt[(L, int(x["n_bins_touched"]), x["is_tie"], x["overlaps_desc"])] += int(x["reads"])
    acc_rows.append([L, r["group"], r["replicate"]] + [a0.get(k, 0) for k in AKEYS] +
                    [f"{100 * a0['assigned'] / a0['eligible_NH1_ge18']:.2f}", len(bc), sum(1 for v in tt.values() if v[0])])

W = lambda name, hdr, data: open(f"{OUT}/{name}", "w").write("\t".join(hdr) + "\n" + "".join("\t".join(map(str, x)) + "\n" for x in data))
W("accounting_summary.tsv", ["library", "group", "replicate", *AKEYS, "pct_eligible_assigned", "bins_with_reads", "transcripts_with_assigned_reads"], acc_rows)
W("tie_summary.tsv", ["library", "group", "replicate", "assigned_reads", "unresolved_tie_reads", "tie_2bin", "tie_3plus",
                      "pct_ties_of_assigned_plus_ties", "ties_PCG_CDS", "ties_lncRNA", "ties_cigar_has_N", "ties_cigar_has_D",
                      "ties_cigar_has_S", "ties_cigar_has_I", "tie_read_length_min/median/max", "tie_max_overlap_min/median/max",
                      "transcripts_with_ties"], tie_rows)
H = ["library", "read", "chrom", "start0", "rna_strand", "cigar", "read_length", "model_id", "n_bins_touched", "n_bins_at_max", "max_overlap", "bin_overlaps_desc"]
W("tie_examples.tsv", H, [[x[k] for k in H] for x in tie_ex])
W("assignment_examples.tsv", H + ["assigned_bin"], [[x[k] for k in H + ["assigned_bin"]] for x in asg_ex])
# overlap pattern summary pooled per library: share of assigned reads touching >= 2 bins, and the margin top - second
ps = []
for L in libs:
    res = [(n, t, v, c) for (l, n, t, v), c in patt.items() if l == L]
    two = sum(c for n, t, v, c in res if t == "False")
    marg = collections.Counter()
    for n, t, v, c in res:
        if t == "False":
            a, b = map(int, v.split(",")[:2]); d = a - b
            marg["1-5" if d <= 5 else "6-20" if d <= 20 else ">20"] += c
    ps.append((L, two, sum(c for n, t, v, c in res if t == "True"), sum(c for n, t, v, c in res if n >= 3),
               marg["1-5"], marg["6-20"], marg[">20"]))
W("overlap_pattern_summary.tsv", ["library", "assigned_reads_touching_2plus_bins", "tie_reads", "reads_touching_3plus_bins",
                                  "assigned_margin_top_minus_second_1to5nt", "margin_6to20nt", "margin_gt20nt"], ps)
# matrices
with gzip.open(f"{OUT}/matrix_bin_assigned_reads.tsv.gz", "wt") as o:
    o.write("bin_id\ttranscript_id\trna_class\tbin_number\tbin_length\t" + "\t".join(libs) + "\n")
    for b in ann:
        o.write("\t".join([b["bin_id"], b["model_id"], b["rna_class"], b["bin_number"], b["bin_length"]] + [str(MB[L].get(b["bin_id"], 0)) for L in libs]) + "\n")
with gzip.open(f"{OUT}/matrix_transcript_totals.tsv.gz", "wt") as o:
    o.write("transcript_id\trna_class\tn_bins\t" + "\t".join(f"{L}_assigned" for L in libs) + "\t" + "\t".join(f"{L}_unresolved_ties" for L in libs) + "\n")
    for m, b in models.items():
        o.write("\t".join([m, b["rna_class"], b["n_bins"]] + [str(MT[L].get(m, (0, 0))[0]) for L in libs] + [str(MT[L].get(m, (0, 0))[1]) for L in libs]) + "\n")
# V5 independent re-read of the written per-library tables: per transcript, sum of bin column == transcript column
v5 = True
for L in libs:
    s = collections.Counter(); t = {}
    for x in csv.DictReader(gzip.open(f"{OUT}/{L}.bin100_table.tsv.gz", "rt"), delimiter="\t"):
        s[x["transcript_id"]] += int(x["bin_assigned_reads"]); t[x["transcript_id"]] = int(x["transcript_assigned_reads"])
    v5 &= all(s[m] == t[m] for m in t) and len(t) == len(models)
gate("V5_written_tables_bin_sum_eq_transcript_total_all_libraries", v5)
W("gates.tsv", ["gate", "result", "detail"], gates)
print(open(f"{OUT}/gates.tsv").read()); print(open(f"{OUT}/accounting_summary.tsv").read()); print(open(f"{OUT}/tie_summary.tsv").read())
print(open(f"{OUT}/overlap_pattern_summary.tsv").read())
ok = all(g[1] == "PASS" for g in gates); print("B3_ALL_GATES_PASS" if ok else "B3_GATE_FAILURE"); sys.exit(0 if ok else 1)
