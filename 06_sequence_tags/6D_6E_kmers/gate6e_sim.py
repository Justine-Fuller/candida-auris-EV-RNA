#!/usr/bin/env python3
"""gate6e_sim.py -- Step 6E gates on the synthetic groups (2026-10-12). Usage: gate6e_sim.py <K3 dir of the SIM run>
G1 planted: GGACTA significant (q < 0.05), higher in EV, in the primary (GC-matched) test (and reported for the GC-class test).
G2 GC shift: the primary test calls NO word (q < 0.05); the naive test (no GC control) is reported to show the shift is real.
G3 null: the primary test calls no word."""
import sys, gzip, csv
K3 = sys.argv[1]; fails = []
def rd(G): return {r["word"]: r for r in csv.DictReader(gzip.open(f"{K3}/{G}.PCG_CDS.words.tsv.gz", "rt"), delimiter="\t")}
def n(R, k): return sum(1 for r in R.values() if r[f"q_{k}"] != "NA" and float(r[f"q_{k}"]) < 0.05)
P = rd("SIM_planted"); g = P["GGACTA"]
ok = float(g["q_primary"]) < 0.05 and g["direction_primary"] == "higher_in_EV"
rank = sorted(P.values(), key=lambda r: float(r["p_primary"])).index(g) + 1
print(("PASS" if ok else "FAIL"), "G1_planted_word_called", f"q_primary={g['q_primary']} rank={rank} q_gcclass={g['q_gcclass']} calls_primary={n(P, 'primary')}")
ok or fails.append("G1")
S = rd("SIM_gcshift"); ok = n(S, "primary") == 0
print(("PASS" if ok else "FAIL"), "G2_gc_shift_not_called_after_matching",
      f"primary={n(S, 'primary')} gcclass={n(S, 'gcclass')} naive={n(S, 'naive')} draw2={n(S, 'draw2')} draw3={n(S, 'draw3')}")
ok or fails.append("G2")
N = rd("SIM_null"); ok = n(N, "primary") == 0
print(("PASS" if ok else "FAIL"), "G3_null_no_calls", f"primary={n(N, 'primary')} gcclass={n(N, 'gcclass')} naive={n(N, 'naive')}")
ok or fails.append("G3")
print("ALL_6E_SIM_GATES_PASS" if not fails else f"6E_SIM_GATES_FAILED {fails}")
