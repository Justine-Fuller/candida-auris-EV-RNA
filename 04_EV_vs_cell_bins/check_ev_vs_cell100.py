#!/usr/bin/env python3
"""check_ev_vs_cell100.py -- INDEPENDENT re-derivation (Python, no edgeR) of everything in ev_vs_cell100.R output that does not
depend on the edgeR fit (2026-10-09). Usage:
  check_ev_vs_cell100.py <EV_tables> <Cell_tables> <annotation> <EV_manifest> <Cell_manifest> <EV_hotspot_dir> <Cell_hotspot_dir> <R_outdir>
Checks per group: counts and transcript reads copied correctly; status re-derived (D5 exclusions; >= 10 physical reads in
>= 2 of 3 EV and >= 2 of 3 Cell; zero in all 6; < 2 remaining bins); mean shares (libraries with 0 transcript reads
excluded), fold enrichment and observed direction; BH q (dataset-wide and within transcript) recomputed from the R p-values;
hotspot flags copied from the group-call files; r2: shares equal within 1e-12 (relative) = "equal", FE 1; category assignment re-derived at 0.01/0.05/0.1. Exit 1 on any mismatch."""
import sys, os, gzip, csv, collections, math
EVD, CD, ANN, MEV, MC, HEV, HC, OUT = sys.argv[1:9]
rows = lambda p: list(csv.DictReader(gzip.open(p, "rt") if p.endswith(".gz") else open(p), delimiter="\t"))
ann = rows(ANN)
ex_r = {a["bin_id"] for a in ann if a["role"] == "competitor_only"}
ex_o = {a["bin_id"] for a in ann if a["same_strand_binned_model_overlap"] != "" and a["role"] != "competitor_only"}
mev, mc = rows(MEV), rows(MC)
def bh(ps):
    m = len(ps); o = sorted(range(m), key=lambda i: ps[i]); q = [0.0] * m; prev = 1.0
    for rank in range(m, 0, -1):
        i = o[rank - 1]; prev = min(prev, ps[i] * m / rank); q[i] = min(prev, 1.0)
    return q
res = []; allok = True
def chk(name, ok, d=""):
    global allok
    res.append((name, "PASS" if ok else "FAIL", str(d))); allok &= ok
for g in dict.fromkeys(r["group"] for r in mev):
    le = [r["library"] for r in mev if r["group"] == g]; lc = [r["library"] for r in mc if r["group"] == g]
    T = {L: rows(f"{(EVD if L in le else CD)}/{L}.bin100_final_table.tsv.gz") for L in lc + le}
    R = rows(f"{OUT}/{g}.ev_vs_cell_bins.tsv.gz")
    bins = [r["bin_id"] for r in T[lc[0]]]
    chk(f"X1_same_bin_order_{g}", [r["bin_id"] for r in R] == bins)
    cnt = {L: [float(r["bin_count"]) for r in T[L]] for L in T}; txr = {L: [float(r["transcript_physical_reads"]) for r in T[L]] for L in T}
    copy_ok = all(float(R[i][f"count_{L}"]) == cnt[L][i] and float(R[i][f"transcript_reads_{L}"]) == txr[L][i] for L in T for i in range(len(bins)))
    chk(f"X2_counts_copied_{g}", copy_ok)
    tid = [r["transcript_id"] for r in T[lc[0]]]
    st = []
    for i, b in enumerate(bins):
        if b in ex_r: s = "excluded_rDNA_competitor_only"
        elif b in ex_o: s = "excluded_same_strand_overlap"
        elif not (sum(txr[L][i] >= 10 for L in le) >= 2 and sum(txr[L][i] >= 10 for L in lc) >= 2): s = "not_tested_read_rule"
        elif sum(cnt[L][i] for L in T) == 0: s = "not_tested_zero_in_all_6"
        else: s = "candidate"
        st.append(s)
    nc = collections.Counter(tid[i] for i, s in enumerate(st) if s == "candidate")
    st = [("tested" if nc[tid[i]] >= 2 else "not_tested_lt2_bins") if s == "candidate" else s for i, s in enumerate(st)]
    chk(f"X3_status_rederived_{g}", st == [r["status"] for r in R], collections.Counter(st)["tested"])
    bad = 0; dirs = []
    for i, r in enumerate(R):
        sh = lambda libs: [cnt[L][i] / txr[L][i] for L in libs if txr[L][i] > 0]
        e, c = sh(le), sh(lc); me = sum(e) / len(e) if e else None; mcv = sum(c) / len(c) if c else None
        eq = me is not None and mcv is not None and abs(me - mcv) <= 1e-12 * max(me, mcv)     # r2: float ties are equal
        d = None if me is None or mcv is None else ("equal" if eq else "higher_in_EV" if me > mcv else "higher_in_Cell")
        dirs.append(d)
        fe = None if me is None or mcv is None else (1.0 if (eq and mcv > 0) else me / mcv if mcv > 0 else (math.inf if me > 0 else None))
        rv = lambda k: None if r[k] == "NA" else float(r[k])
        if (r["direction_observed"] if r["direction_observed"] != "NA" else None) != d: bad += 1
        for k, v in (("mean_share_EV", me), ("mean_share_Cell", mcv)):
            if (v is None) != (rv(k) is None) or (v is not None and abs(v - rv(k)) > 1e-9): bad += 1
        f2 = r["fold_enrichment_EV_over_Cell"]; f2 = None if f2 == "NA" else (math.inf if f2 == "Inf" else float(f2))
        if (fe is None) != (f2 is None) or (fe is not None and fe != f2 and abs(fe - f2) > 1e-9 * max(1, abs(fe))): bad += 1
    chk(f"X4_shares_FE_direction_{g}", bad == 0, bad)
    te = [i for i, s in enumerate(st) if s == "tested"]
    p = [float(R[i]["p"]) for i in te]; q = bh(p)
    chk(f"X5_BH_dataset_{g}", all(abs(q[j] - float(R[i]["q_dataset"])) < 1e-12 for j, i in enumerate(te)), len(te))
    byt = collections.defaultdict(list)
    for j, i in enumerate(te): byt[tid[i]].append(j)
    qw_ok = True
    for t, js in byt.items():
        qq = bh([p[j] for j in js])
        qw_ok &= all(abs(qq[k] - float(R[te[j]]["q_within_transcript"])) < 1e-12 for k, j in enumerate(js))
    chk(f"X6_BH_within_transcript_{g}", qw_ok)
    he = {r["bin_id"]: r for r in rows(f"{HEV}/{g}.group_calls.tsv.gz")}; hc = {r["bin_id"]: r for r in rows(f"{HC}/{g}.group_calls.tsv.gz")}
    cat_ok = True
    for lv in ("0.01", "0.05", "0.1"):
        for i, r in enumerate(R):
            evh = he[r["bin_id"]][f"group_call_2of3_q{lv}"] == "TRUE"; ceh = hc[r["bin_id"]][f"group_call_2of3_q{lv}"] == "TRUE"
            tested = st[i] == "tested"; sig = tested and float(r["q_dataset"]) < float(lv); d = dirs[i]
            if evh:
                c = ("EV_hotspot__comparison_not_tested" if not tested else
                     ("EV_hotspot__also_Cell_hotspot__EV_share_higher" if ceh else "EV_hotspot__EV_share_higher") if sig and d == "higher_in_EV" else
                     "EV_hotspot__Cell_share_higher" if sig and d == "higher_in_Cell" else
                     "EV_hotspot__also_Cell_hotspot__no_significant_difference" if ceh else "EV_hotspot__no_Cell_hotspot__no_significant_difference")
            else:
                c = ("Cell_hotspot__Cell_share_higher" if ceh and sig and d == "higher_in_Cell" else "Cell_hotspot__other" if ceh else
                     f"no_hotspot__share_{d}" if sig else "no_hotspot__no_difference")
            cat_ok &= (r[f"EV_hotspot_q{lv}"] == ("TRUE" if evh else "FALSE") and r[f"Cell_hotspot_q{lv}"] == ("TRUE" if ceh else "FALSE")
                       and r[f"category_q{lv}"] == c)
    chk(f"X7_hotspot_flags_and_categories_{g}", cat_ok)
with open(f"{OUT}/independent_check.tsv", "w") as o:
    o.write("check\tresult\tdetail\n"); [o.write("\t".join(x) + "\n") for x in res]
print(open(f"{OUT}/independent_check.tsv").read()); print("INDEPENDENT_CHECK_PASS" if allok else "INDEPENDENT_CHECK_FAIL"); sys.exit(0 if allok else 1)
