#!/usr/bin/env python3
"""cca2_compare.py -- CCA recount (2026-10-13): gates, old-vs-new CCA table, and T1-format JSON copies for the FL rebuild.
Usage: cca2_compare.py <R1 cca2 dir> <T1_all dir (frozen)> <manifest_all36.tsv> <outdir>
Gates: (G1) the old-rule fields reproduce the frozen T1_all JSON for every library and family; (G2) INVARIANT_MISMATCH = 0
in every library. Writes:
  cca_table_by_compartment.tsv : per compartment x family, reads reaching the 3' end and CCA reads, old vs new (+ C / CC,
                                 templated-ambiguous), % CCA of 3'-end reads;
  T1_all_cca2_incl/  : copies of the T1_all JSONs with trna_families[F].CCA_reads = new CCA reads INCLUDING templated-ambiguous;
  T1_all_cca2_excl/  : same, EXCLUDING templated-ambiguous (strict). Every other field is unchanged (input for fl_classify_r3.py)."""
import sys, os, json, csv, collections
R1, T1, MAN, OUT = sys.argv[1:5]; os.makedirs(OUT, exist_ok=False)
man = list(csv.DictReader(open(MAN), delimiter="\t"))
OLDF = ["reads", "reads_NH1", "reach_3prime", "clipped_at_3prime", "CCA_reads", "reads_bigN_excluded_from_CCA_check"]
bad = []; inv = 0; tab = collections.defaultdict(collections.Counter)
for v in ("incl", "excl"): os.makedirs(f"{OUT}/T1_all_cca2_{v}")
for r in man:
    L = r["library"]; n = json.load(open(f"{R1}/{L}.cca2.json")); o = json.load(open(f"{T1}/{L}.trna.json"))
    N, O = n["trna_families"], o["trna_families"]
    if set(N) != set(O): bad.append((L, "family set"))
    for F in O:
        for k in OLDF:
            if N.get(F, {}).get(k, 0) != O[F].get(k, 0): bad.append((L, F, k, N.get(F, {}).get(k, 0), O[F].get(k, 0)))
    inv += n["accounting"].get("INVARIANT_MISMATCH", 0)
    for F, d in N.items():
        t = tab[(r["compartment"], F)]
        for k in ("reach_3prime", "CCA_reads", "reach_3prime_new", "tail_CCA", "tail_CCA_non_templated", "tail_CCA_templated_ambiguous", "tail_CC", "tail_C"):
            t[k] += d.get(k, 0)
    for v, key in (("incl", "tail_CCA"), ("excl", "tail_CCA_non_templated")):
        o2 = json.loads(json.dumps(o))
        for F in o2["trna_families"]: o2["trna_families"][F]["CCA_reads"] = N.get(F, {}).get(key, 0)
        o2["cca_recount_note"] = f"CCA_reads replaced by cca2 {key} (2026-10-13); all other fields = frozen T1_all"
        json.dump(o2, open(f"{OUT}/T1_all_cca2_{v}/{L}.trna.json", "w"), sort_keys=True)
pct = lambda a, b: f"{100 * a / b:.1f}" if b else "NA"
with open(f"{OUT}/cca_table_by_compartment.tsv", "w") as w:
    w.write("compartment\tfamily\treach_3prime_old\tCCA_old\tpct_CCA_old\treach_3prime_new\tCCA_new\tCCA_new_non_templated\t"
            "CCA_new_templated_ambiguous\tpct_CCA_new\ttail_CC\ttail_C\n")
    for (c, F), t in sorted(tab.items()):
        w.write(f"{c}\t{F}\t{t['reach_3prime']}\t{t['CCA_reads']}\t{pct(t['CCA_reads'], t['reach_3prime'])}\t{t['reach_3prime_new']}\t"
                f"{t['tail_CCA']}\t{t['tail_CCA_non_templated']}\t{t['tail_CCA_templated_ambiguous']}\t{pct(t['tail_CCA'], t['reach_3prime_new'])}\t"
                f"{t['tail_CC']}\t{t['tail_C']}\n")
with open(f"{OUT}/gates.tsv", "w") as w:
    w.write(f"G1_old_rule_reproduces_T1_all\t{'PASS' if not bad else 'FAIL'}\t{len(bad)} differences {bad[:5]}\n")
    w.write(f"G2_invariant_no_past_end_new_equals_old\t{'PASS' if inv == 0 else 'FAIL'}\t{inv}\n")
print(open(f"{OUT}/gates.tsv").read()); print("CCA2_COMPARE_DONE" if not bad and inv == 0 else "CCA2_GATE_FAILURE")
