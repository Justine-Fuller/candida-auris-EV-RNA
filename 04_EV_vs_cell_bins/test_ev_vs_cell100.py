#!/usr/bin/env python3
"""Known-answer test for ev_vs_cell100.R (2026-10-09). Usage: test_ev_vs_cell100.py <scriptdir> <workdir> <Rscript>
Synthetic group "G_test", 3 Cell + 3 EV libraries, 6-bin transcripts, replicate overdispersion (gamma noise on bin weights):
  T001-T300 null (same share profile in EV and Cell); D001-D030 bin 3 weight x6 in EV only (true EV-higher bins);
  L001-L010 only 5 reads in every EV library (fail read rule); Z001 bins 5-6 never get reads (not_tested_zero);
  O001 only bin 1 has reads (bin 1 -> not_tested_lt2_bins, others zero); R001 competitor_only (excluded_rDNA);
  V001 bin 2 flagged same-strand overlap (excluded_overlap).
Checks: exact status counts; >= 27/30 planted bins called at q<0.05 with direction higher_in_EV; <= 5% of null-transcript bins
called at q<0.05; observed shares / fold enrichment / direction recomputed here == R output; categories for hand-flagged bins;
all R gates PASS."""
import sys, os, subprocess, gzip, csv, random
S, W, RS = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(W, exist_ok=True); rnd = random.Random(20261009)
CELL = ["C1", "C2", "C3"]; EV = ["E1", "E2", "E3"]; LIBS = CELL + EV; NB = 6
tx = {}
for i in range(1, 301): tx[f"T{i:03d}"] = "null"
for i in range(1, 31): tx[f"D{i:03d}"] = "du"
for i in range(1, 11): tx[f"L{i:03d}"] = "low"
for t in ("Z001", "O001", "R001", "V001"): tx[t] = t
def multinom(n, w):
    s = sum(w); out = [0] * len(w)
    for _ in range(n):
        r = rnd.random() * s; acc = 0
        for j, x in enumerate(w):
            acc += x
            if r < acc: out[j] += 1; break
    return out
counts = {}
for t, kind in tx.items():
    base = [0.5 + rnd.random() for _ in range(NB)]
    if kind == "Z001": base[4] = base[5] = 0
    if kind == "O001": base = [1, 0, 0, 0, 0, 0]
    for L in LIBS:
        w = [b * rnd.gammavariate(20, 1 / 20) for b in base]
        if kind == "du" and L in EV: w[2] *= 6
        n = 5 if (kind == "low" and L in EV) else rnd.randint(120, 200)
        counts[(t, L)] = multinom(n, w)
with open(f"{W}/ann.tsv", "w") as o:
    o.write("bin_id\trole\tsame_strand_binned_model_overlap\n")
    for t in tx:
        for b in range(1, NB + 1):
            o.write(f"{t}|bin{b:03d}\t{'competitor_only' if t == 'R001' else 'ordinary'}\t{'X' if (t == 'V001' and b == 2) else ''}\n")
for comp, libs in (("EV", EV), ("Cell", CELL)):
    os.makedirs(f"{W}/{comp}", exist_ok=True)
    with open(f"{W}/man_{comp}.tsv", "w") as o:
        o.write("library\tbam\tkind\tcompartment\tgroup\treplicate\n")
        for i, L in enumerate(libs, 1): o.write(f"{L}\tx\tsRNA\t{comp}\tG_test\t{i}\n")
    for L in libs:
        with gzip.open(f"{W}/{comp}/{L}.bin100_final_table.tsv.gz", "wt") as o:
            o.write("transcript_id\tgene_id\trna_class\trole\tbin_id\tbin_number\tn_bins\tbin_tx_start0\tbin_tx_end\tbin_length\tbin_count\ttranscript_physical_reads\n")
            for t in tx:
                c = counts[(t, L)]
                for b in range(1, NB + 1):
                    o.write(f"{t}\t{t}\tPCG_CDS\t{'competitor_only' if t == 'R001' else 'ordinary'}\t{t}|bin{b:03d}\t{b}\t{NB}\t{(b-1)*100}\t{b*100}\t100\t{c[b-1]}\t{sum(c)}\n")
def hot(path, ev):
    os.makedirs(path, exist_ok=True)
    with gzip.open(f"{path}/G_test.group_calls.tsv.gz", "wt") as o:
        o.write("bin_id\tgroup_call_2of3_q0.01\tgroup_call_2of3_q0.05\tgroup_call_2of3_q0.1\n")
        for t in tx:
            for b in range(1, NB + 1):
                f = (ev and t.startswith("D") and b == 3) or (not ev and t in [f"D{i:03d}" for i in range(1, 11)] and b == 3) or (t == "T001" and b == 1)
                o.write(f"{t}|bin{b:03d}\t" + "\t".join(["TRUE" if f else "FALSE"] * 3) + "\n")
hot(f"{W}/hotEV", True); hot(f"{W}/hotCell", False)
OUT = f"{W}/out"
p = subprocess.run([RS, f"{S}/ev_vs_cell100.R", f"{W}/EV", f"{W}/Cell", f"{W}/ann.tsv", f"{W}/man_EV.tsv", f"{W}/man_Cell.tsv",
                    f"{W}/hotEV", f"{W}/hotCell", OUT, "real", "1"], capture_output=True, text=True)
ok = p.returncode == 0
if not ok: print("R FAILED\n", p.stdout[-3000:], p.stderr[-3000:])
else:
    R = {r["bin_id"]: r for r in csv.DictReader(gzip.open(f"{OUT}/G_test.ev_vs_cell_bins.tsv.gz", "rt"), delimiter="\t")}
    st = {}
    for r in R.values(): st[r["status"]] = st.get(r["status"], 0) + 1
    exp_st = {"tested": 330 * 6 + 4 + 5, "not_tested_read_rule": 60, "not_tested_zero_in_all_6": 2 + 5, "not_tested_lt2_bins": 1,
              "excluded_rDNA_competitor_only": 6, "excluded_same_strand_overlap": 1}
    c1 = st == exp_st; print("status counts exact:", c1, "" if c1 else st); ok &= c1
    du = [R[f"D{i:03d}|bin003"] for i in range(1, 31)]
    hit = sum(1 for r in du if r["q_dataset"] != "NA" and float(r["q_dataset"]) < 0.05 and r["direction_observed"] == "higher_in_EV")
    c2 = hit >= 27; print(f"planted EV-higher bins called (q<0.05, higher_in_EV): {hit}/30", c2); ok &= c2
    nb = [r for r in R.values() if r["transcript_id"].startswith("T") and r["status"] == "tested"]
    fp = sum(1 for r in nb if float(r["q_dataset"]) < 0.05)
    c3 = fp <= 0.05 * len(nb); print(f"null-transcript bins called at q<0.05: {fp}/{len(nb)}", c3); ok &= c3
    bad = 0
    for bid, r in R.items():
        t, b = bid.split("|bin"); b = int(b) - 1
        sh = {L: counts[(t, L)][b] / sum(counts[(t, L)]) for L in LIBS if sum(counts[(t, L)]) > 0}
        mE = [sh[L] for L in EV if L in sh]; mC = [sh[L] for L in CELL if L in sh]
        me = sum(mE) / len(mE) if mE else None; mc = sum(mC) / len(mC) if mC else None
        d = None if me is None or mc is None else ("equal" if abs(me - mc) <= 1e-12 * max(me, mc) else "higher_in_EV" if me > mc else "higher_in_Cell")
        if (r["direction_observed"] if r["direction_observed"] != "NA" else None) != d: bad += 1
        if me is not None and abs(float(r["mean_share_EV"]) - me) > 1e-9: bad += 1
        if mc is not None and abs(float(r["mean_share_Cell"]) - mc) > 1e-9: bad += 1
    c4 = bad == 0; print("shares and direction recomputed independently match:", c4, bad); ok &= c4
    cat = lambda bid: R[bid]["category_q0.05"]
    exp_cat_ok = all(cat(f"D{i:03d}|bin003") == ("EV_hotspot__also_Cell_hotspot__EV_share_higher" if i <= 10 else "EV_hotspot__EV_share_higher")
                     for i in range(1, 31) if float(R[f"D{i:03d}|bin003"]["q_dataset"]) < 0.05)
    t1 = cat("T001|bin001"); t1_ok = t1 in ("EV_hotspot__also_Cell_hotspot__no_significant_difference", "EV_hotspot__also_Cell_hotspot__EV_share_higher",
                                            "EV_hotspot__Cell_share_higher")
    lo = cat("L001|bin001") == "no_hotspot__no_difference" and R["L001|bin001"]["status"] == "not_tested_read_rule"
    c5 = exp_cat_ok and t1_ok and lo; print("categories as expected:", c5, t1); ok &= c5
    g = list(csv.DictReader(open(f"{OUT}/gates.tsv"), delimiter="\t"))
    c6 = all(x["result"] == "PASS" for x in g); print("R gates all PASS:", c6, len(g)); ok &= c6
print("ALL_EV_VS_CELL_TESTS_PASS" if ok else "EV_VS_CELL_TESTS_FAILED"); sys.exit(0 if ok else 1)
