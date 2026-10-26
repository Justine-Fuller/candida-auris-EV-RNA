#!/usr/bin/env python3
"""Step 8 (v2): prioritise EV hotspot regions (plan 00_PLAN_DRAFT_v2_hotspot_prioritisation_2026-10-13.md, md5 839b2104…,
approved by Justine 2026-10-13: Q1 a regions; Q2 F2 only; F3/F4/F5 = columns; Q3 a; Q4 a all; Q5 a table + FASTA).

Unit (Q1 a): EV hotspot region = neighbouring EV hotspot bins (N5 group call, >= 2/3 replicates, q < 0.05) of one transcript
merged, exactly as Step 5A.
Filter (F2 only): a region passes if it contains >= 1 EV-selective bin (EV hotspot AND E1 category_q0.05 ending
"EV_share_higher", i.e. q_dataset < 0.05 with the higher observed share in EVs) -- same bin rule as Step 5A.
Sort (r2; Justine 2026-10-13 A3 + B2): the region keeps its full hotspot extent (table + FASTA), but is ranked by the
EV-selective bin's own observed EV share / Cell share (E1 fold_enrichment_EV_over_Cell; if a region has > 1 EV-selective bin,
the largest). Ties (incl. the Inf group, Cell share 0) are broken by EV read depth = mean EV reads in the region (F5), then
region_id. The whole-region ratio (sum of bin counts / transcript reads; G1) and its strict no-tie-read version are kept as
columns with alternative ranks (r1 sorted by the whole-region ratio and stopped at the gate "region ratio > 1": 2 regions 0.99).
Direction and shares come from observed proportions only.
Columns (not filters): F3 recurrence (groups with an overlapping EV hotspot region), F4 Step 5B annotation, F5 depth.
Sequences: transcript strand, from the frozen genome + the frozen bin annotation's genomic blocks.

usage: s8_prioritise.py selftest <workdir>
       s8_prioritise.py run N5 E1 N4_EV N4_CELL ANNOT GENOME S5A_REGIONS S5B_ANNOT OUT
"""
import csv, gzip, json, math, os, sys, collections

GROUPS = ["B8441_untreated", "B8441_caspofungin", "MMC1_untreated", "MMC1_caspofungin"]
EXPECTED_SEL = {"B8441_untreated": 8, "B8441_caspofungin": 34, "MMC1_untreated": 6, "MMC1_caspofungin": 37}
ANN_COLS = ["V3_product", "eggNOG_Preferred_name", "eggNOG_Description", "eggNOG_COG_category", "eggNOG_KEGG_Pathway",
            "Candida_albicans_orthologs", "Saccharomyces_cerevisiae_orthologs", "SignalP_prediction", "DeepTMHMM_has_TM",
            "PredGPI_prediction", "n_GO_terms", "in_n_groups_EV_hotspot", "lncRNA_nearest_PCG", "lncRNA_nearest_PCG_distance_bp",
            "nearest_PCG_V3_product", "nearest_PCG_eggNOG_Description"]
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rd(path):
    f = gzip.open(path, "rt") if path.endswith(".gz") else open(path)
    return csv.DictReader(f, delimiter="\t")


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name: seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    if name: seqs[name] = "".join(buf)
    return seqs


def merge_runs(bins):
    """sorted bin numbers -> list of runs of neighbouring bins (numbers differ by 1), as Step 5A"""
    out, cur = [], []
    for b in sorted(bins):
        if cur and b == cur[-1] + 1: cur.append(b)
        else:
            if cur: out.append(cur)
            cur = [b]
    if cur: out.append(cur)
    return out


def mean_share(counts, totals):
    """G1: replicate share undefined when the transcript has 0 reads in that library -> left out of the mean"""
    v = [c / t for c, t in zip(counts, totals) if t > 0]
    return (sum(v) / len(v)) if v else float("nan"), len(v)


def ratio(ev, cell):
    if math.isnan(ev) or math.isnan(cell): return float("nan")
    if cell == 0: return float("inf") if ev > 0 else float("nan")
    return ev / cell


def seq_of(genome, chrom, strand, blocks):
    s = "".join(genome[chrom][a:b] for a, b in blocks)  # blocks are listed transcript 5'->3'
    if strand == "-":
        s = "".join(genome[chrom][a:b].translate(COMP)[::-1] for a, b in blocks)
    return s.upper()


def full(x):
    if isinstance(x, float) and math.isnan(x): return "NA"
    if isinstance(x, float) and math.isinf(x): return "Inf"
    return repr(x)


def parse_blocks(txt):
    return [tuple(int(x) for x in blk.split("-")) for blk in txt.split(";") if blk]


def fmt(x):
    if isinstance(x, float):
        if math.isnan(x): return "NA"
        if math.isinf(x): return "Inf"
        return f"{x:.6g}"
    return str(x)


def main(N5, E1, N4E, N4C, ANNOT, GENOME, S5A, S5B, OUT, check_5a_counts=True):
    os.makedirs(OUT, exist_ok=False)
    gates = []
    def gate(name, ok, detail=""):
        gates.append((name, "PASS" if ok else "FAIL", detail))
        print(f"GATE {name}: {'PASS' if ok else 'FAIL'} {detail}", flush=True)
        return ok

    ann = {r["bin_id"]: r for r in rd(ANNOT)}
    genome = read_fasta(GENOME)

    # ---- per group: bins, hotspot regions, EV-selective bins
    hot_regions = {}  # group -> list of (transcript, first, last)
    data = {}
    for G in GROUPS:
        n5 = {r["bin_id"]: r for r in rd(f"{N5}/{G}.group_calls.tsv.gz")}
        e1 = {r["bin_id"]: r for r in rd(f"{E1}/{G}.ev_vs_cell_bins.tsv.gz")}
        assert set(n5) == set(e1), G
        assert all((e1[b]["EV_hotspot_q0.05"] == "TRUE") == (n5[b]["group_call_2of3_q0.05"] == "TRUE") for b in n5), G
        hdr = list(next(iter(e1.values())).keys())
        cell_libs = [c[6:] for c in hdr if c.startswith("count_") and "_sRNA" in c]
        ev_libs = [c[6:] for c in hdr if c.startswith("count_") and "_EV" in c]
        assert len(cell_libs) == 3 and len(ev_libs) == 3, (G, cell_libs, ev_libs)
        # recompute E1's observed mean shares from counts (G1) and check they match exactly
        bad = 0
        for b, r in e1.items():
            for side, libs in (("EV", ev_libs), ("Cell", cell_libs)):
                m, _ = mean_share([int(r["count_" + l]) for l in libs], [int(r["transcript_reads_" + l]) for l in libs])
                ref = r[f"mean_share_{side}"]
                if ref == "NA": bad += not math.isnan(m)
                elif math.isnan(m) or abs(m - float(ref)) > 1e-9 * max(1, abs(m)): bad += 1
        gate(f"{G}: recomputed bin mean shares == E1 (G1)", bad == 0, f"mismatches={bad} bins={len(e1)}")
        hot = collections.defaultdict(list); sel = collections.defaultdict(set)
        for b, r in n5.items():
            if r["group_call_2of3_q0.05"] == "TRUE": hot[r["transcript_id"]].append(int(r["bin_number"]))
            if e1[b]["EV_hotspot_q0.05"] == "TRUE" and e1[b]["category_q0.05"].endswith("EV_share_higher"):
                sel[r["transcript_id"]].add(int(r["bin_number"]))
        byid = {(r["transcript_id"], int(r["bin_number"])): r for r in e1.values()}
        hot_regions[G] = [(t, run[0], run[-1]) for t, bl in hot.items() for run in merge_runs(bl)]
        sel_regions = [(t, run[0], run[-1]) for t, bl in sel.items() for run in merge_runs(bl)]
        data[G] = dict(hot=hot, sel=sel, byid=byid, ev_libs=ev_libs, cell_libs=cell_libs, sel_regions=sel_regions)

    # ---- gate: regions reproduce Step 5A exactly (both lists)
    s5a = collections.defaultdict(set)
    for r in rd(S5A):
        s5a[(r["group"], r["list"])].add((r["transcript_id"], int(r["first_bin"]), int(r["last_bin"])))
    for G in GROUPS:
        gate(f"{G}: EV hotspot regions == Step 5A", set(hot_regions[G]) == s5a[(G, "EV_hotspot")],
             f"n={len(hot_regions[G])} 5A={len(s5a[(G, 'EV_hotspot')])}")
        gate(f"{G}: EV-selective-bin regions == Step 5A", set(data[G]["sel_regions"]) == s5a[(G, "EV_selective")],
             f"n={len(data[G]['sel_regions'])} 5A={len(s5a[(G, 'EV_selective')])}")

    # ---- F2 filter: hotspot regions with >= 1 EV-selective bin
    passing = []
    for G in GROUPS:
        for (t, a, b) in hot_regions[G]:
            sb = sorted(x for x in data[G]["sel"].get(t, ()) if a <= x <= b)
            if sb: passing.append((G, t, a, b, sb))
    n_pass = collections.Counter(p[0] for p in passing)
    if check_5a_counts:
        gate("passing regions per group == 8/34/6/37", all(n_pass[G] == EXPECTED_SEL[G] for G in GROUPS),
             " ".join(f"{G}={n_pass[G]}" for G in GROUPS))
    gate("passing regions per group == Step 5A EV-selective region counts",
         all(n_pass[G] == len(s5a[(G, "EV_selective")]) for G in GROUPS), "")
    gate("every EV-selective bin lies in a passing region",
         sum(len(p[4]) for p in passing) == sum(len(s) for G in GROUPS for s in data[G]["sel"].values()), "")

    # ---- strict (no tie) counts from the N4 per-library tables, only for bins in passing regions
    need = collections.defaultdict(set)
    for G, t, a, b, sb in passing:
        for x in range(a, b + 1): need[G].add(f"{t}|bin{x:03d}")
    allneed = set().union(*need.values()) if need else set()
    strict, tie, bcount, tphys = {}, {}, {}, {}
    for d in (N4E, N4C):
        for fn in sorted(os.listdir(d)):
            if not fn.endswith(".bin100_final_table.tsv.gz"): continue
            lib = fn.split(".bin100_final_table")[0]
            for r in rd(f"{d}/{fn}"):
                if r["bin_id"] in allneed:
                    k = (lib, r["bin_id"])
                    strict[k] = int(r["bin_reads_strict_greatest_overlap"]); tie[k] = int(r["bin_tie_plus1"])
                    bcount[k] = int(r["bin_count"]); tphys[k] = int(r["transcript_physical_reads"])

    # ---- build rows
    rows, fasta, n4bad = [], [], 0
    hot_iv = {G: collections.defaultdict(list) for G in GROUPS}
    for G in GROUPS:
        for t, a, b in hot_regions[G]: hot_iv[G][t].append((a, b))
    pass_iv = collections.defaultdict(list)
    for G, t, a, b, sb in passing: pass_iv[t].append((G, a, b))
    s5b = {}
    for r in rd(S5B):
        if r["list"] == "EV_hotspot": s5b[(r["group"], r["gene_id"])] = r
    for G, t, a, b, sb in passing:
        D = data[G]; byid = D["byid"]
        bins = [byid[(t, x)] for x in range(a, b + 1)]
        r0 = bins[0]
        tot = {l: int(r0["transcript_reads_" + l]) for l in D["ev_libs"] + D["cell_libs"]}
        cnt = {l: sum(int(r["count_" + l]) for r in bins) for l in D["ev_libs"] + D["cell_libs"]}
        bid = [r["bin_id"] for r in bins]
        for l in D["ev_libs"] + D["cell_libs"]:
            for r in bins:
                k = (l, r["bin_id"])
                if k not in bcount or bcount[k] != int(r["count_" + l]) or tphys[k] != tot[l]: n4bad += 1
        st = {l: sum(strict[(l, x)] for x in bid) for l in tot}
        ti = {l: sum(tie[(l, x)] for x in bid) for l in tot}
        evs, nev = mean_share([cnt[l] for l in D["ev_libs"]], [tot[l] for l in D["ev_libs"]])
        ces, nce = mean_share([cnt[l] for l in D["cell_libs"]], [tot[l] for l in D["cell_libs"]])
        evs_s, _ = mean_share([st[l] for l in D["ev_libs"]], [tot[l] for l in D["ev_libs"]])
        ces_s, _ = mean_share([st[l] for l in D["cell_libs"]], [tot[l] for l in D["cell_libs"]])
        selr = [byid[(t, x)] for x in sb]
        sc = {l: sum(int(r["count_" + l]) for r in selr) for l in tot}
        sev, _ = mean_share([sc[l] for l in D["ev_libs"]], [tot[l] for l in D["ev_libs"]])
        sce, _ = mean_share([sc[l] for l in D["cell_libs"]], [tot[l] for l in D["cell_libs"]])
        sel_fe = max(float("inf") if r["fold_enrichment_EV_over_Cell"] == "Inf" else float(r["fold_enrichment_EV_over_Cell"]) for r in selr)
        sel_q = min(float(r["q_dataset"]) for r in selr)
        a0, b0 = int(r0["bin_tx_start0"]), int(bins[-1]["bin_tx_end"])
        A = ann[r0["bin_id"]]; chrom, strand = A["chrom"], A["strand"]
        blocks = []
        for r in bins:
            for blk in parse_blocks(ann[r["bin_id"]]["genomic_blocks_5to3"]):
                if blocks and ((strand == "+" and blocks[-1][1] == blk[0]) or (strand == "-" and blocks[-1][0] == blk[1])):
                    blocks[-1] = (blocks[-1][0], blk[1]) if strand == "+" else (blk[0], blocks[-1][1])
                else: blocks.append(blk)
        seq = seq_of(genome, chrom, strand, blocks)
        rec = sorted({g for g in GROUPS for (x, y) in hot_iv[g].get(t, []) if x <= b and a <= y})
        recp = sorted({g for (g, x, y) in pass_iv[t] if x <= b and a <= y})
        an = s5b.get((G, r0["gene_id"]))
        rid = f"{G}|{t}|bin{a:03d}-{b:03d}"
        row = dict(region_id=rid, group=G, transcript_id=t, gene_id=r0["gene_id"], rna_class=r0["rna_class"],
                   rdna_flag=A.get("rdna_flag", ""), chrom=chrom, strand=strand, first_bin=a, last_bin=b,
                   n_bins_in_region=b - a + 1, n_bins_transcript=int(r0["n_bins"]), region_tx_start0=a0, region_tx_end=b0,
                   region_length=b0 - a0, genomic_blocks_5to3=";".join(f"{x}-{y}" for x, y in blocks),
                   ev_selective_bins=";".join(map(str, sb)), n_ev_selective_bins=len(sb),
                   region_mean_share_EV=evs, region_mean_share_Cell=ces, n_EV_proportions=nev, n_Cell_proportions=nce,
                   region_FE_EV_over_Cell=ratio(evs, ces),
                   region_FE_strict_no_tie_reads=ratio(evs_s, ces_s),
                   region_tie_plus1_EV_total=sum(ti[l] for l in D["ev_libs"]),
                   region_tie_plus1_Cell_total=sum(ti[l] for l in D["cell_libs"]),
                   ev_selective_bin_max_FE=sel_fe, ev_selective_bin_min_q_dataset=sel_q,
                   ev_selective_bins_combined_FE=ratio(sev, sce),
                   ev_selective_bins_mean_EV_reads=sum(sc[l] for l in D["ev_libs"]) / 3,
                   F3_n_groups_with_EV_hotspot_region=len(rec), F3_groups_with_EV_hotspot_region=";".join(rec),
                   F3_n_groups_passing=len(recp), F3_groups_passing=";".join(recp),
                   F5_mean_EV_reads_in_region=sum(cnt[l] for l in D["ev_libs"]) / 3,
                   F5_mean_Cell_reads_in_region=sum(cnt[l] for l in D["cell_libs"]) / 3,
                   F5_mean_EV_transcript_reads=sum(tot[l] for l in D["ev_libs"]) / 3,
                   F5_mean_Cell_transcript_reads=sum(tot[l] for l in D["cell_libs"]) / 3,
                   GC_percent=100 * sum(c in "GC" for c in seq) / len(seq) if seq else float("nan"),
                   F4_in_step5B=an is not None, **{"F4_" + c: (an[c] if an and c in an else "") for c in ANN_COLS},
                   sequence=seq)
        row["per_library_counts"] = ";".join(f"{l}={cnt[l]}/{tot[l]}" for l in D["cell_libs"] + D["ev_libs"])
        rows.append(row)
    gate("N4 bin counts / transcript reads == E1 for every region bin", n4bad == 0, f"mismatches={n4bad}")
    gate("sequence length == region length (all)", all(len(r["sequence"]) == r["region_length"] for r in rows), f"n={len(rows)}")
    gate("sequences contain only ACGT", all(set(r["sequence"]) <= set("ACGT") for r in rows), "")
    gate("ranking measure (EV-selective bin FE) > 1 for every region", all(r["ev_selective_bin_max_FE"] > 1 for r in rows),
         f"min={min((r['ev_selective_bin_max_FE'] for r in rows), default=float('nan'))}")
    gate("every passing gene found in Step 5B annotation", all(r["F4_in_step5B"] for r in rows),
         f"missing={sum(not r['F4_in_step5B'] for r in rows)}")

    # ---- sort (Q3 a) + alternative ranks
    def key(f, depth="F5_mean_EV_reads_in_region"):
        return lambda r: (-(r[f] if not math.isnan(r[f]) else -1), -r[depth], r["region_id"])
    rows.sort(key=key("ev_selective_bin_max_FE"))
    for i, r in enumerate(rows, 1): r["rank_overall"] = i
    for G in GROUPS:
        for i, r in enumerate([r for r in rows if r["group"] == G], 1): r["rank_in_group"] = i
    for f, n in (("region_FE_EV_over_Cell", "alt_rank_whole_region_FE"), ("region_FE_strict_no_tie_reads", "alt_rank_whole_region_FE_strict"),
                 ("ev_selective_bins_combined_FE", "alt_rank_selective_bins_combined_FE")):
        for i, r in enumerate(sorted(rows, key=key(f)), 1): r[n] = i
    for i, r in enumerate(sorted(rows, key=key("ev_selective_bin_max_FE", "ev_selective_bins_mean_EV_reads")), 1): r["alt_rank_tiebreak_selective_bin_depth"] = i
    first = ["rank_overall", "rank_in_group", "ev_selective_bin_max_FE", "F5_mean_EV_reads_in_region", "alt_rank_whole_region_FE",
             "alt_rank_whole_region_FE_strict", "alt_rank_selective_bins_combined_FE", "alt_rank_tiebreak_selective_bin_depth"]
    cols = first + [c for c in rows[0] if c not in first] if rows else first
    with open(f"{OUT}/prioritised_regions_q0.05.tsv", "w") as w:
        w.write("\t".join(cols) + "\n")
        for r in rows: w.write("\t".join(fmt(r[c]) for c in cols) + "\n")
    with open(f"{OUT}/prioritised_regions_q0.05.fa", "w") as w:
        for r in rows:
            w.write(f">{r['region_id']} gene={r['gene_id']} class={r['rna_class']} {r['chrom']}:{r['genomic_blocks_5to3']}({r['strand']}) "
                    f"tx={r['region_tx_start0']}-{r['region_tx_end']} len={r['region_length']} rank={r['rank_overall']} "
                    f"selective_bins={r['ev_selective_bins']} selective_bin_FE={fmt(r['ev_selective_bin_max_FE'])}\n")
            for i in range(0, len(r["sequence"]), 60): w.write(r["sequence"][i:i + 60] + "\n")
    multi = [r for r in rows if r["n_bins_in_region"] > 1]
    summ = dict(passing_per_group={G: n_pass[G] for G in GROUPS}, total=len(rows),
                by_class={G: dict(collections.Counter(r["rna_class"] for r in rows if r["group"] == G)) for G in GROUPS},
                multi_bin_regions=len(multi),
                multi_bin_regions_with_internal_tie_reads=sum(1 for r in multi if r["region_tie_plus1_EV_total"] + r["region_tie_plus1_Cell_total"] > 0),
                recurrent_ge2_groups=sum(r["F3_n_groups_with_EV_hotspot_region"] >= 2 for r in rows),
                n_inf_selective_bin_FE=sum(math.isinf(r["ev_selective_bin_max_FE"]) for r in rows),
                n_whole_region_FE_le1=sum(not r["region_FE_EV_over_Cell"] > 1 for r in rows),
                whole_region_FE_le1=[(r["region_id"], fmt(r["region_FE_EV_over_Cell"])) for r in rows if not r["region_FE_EV_over_Cell"] > 1],
                regions_with_gt1_selective_bin=sum(r["n_ev_selective_bins"] > 1 for r in rows),
                rank_changes_if_combined_selective_FE=sum(r["rank_overall"] != r["alt_rank_selective_bins_combined_FE"] for r in rows),
                rank_changes_if_tiebreak_selective_bin_depth=sum(r["rank_overall"] != r["alt_rank_tiebreak_selective_bin_depth"] for r in rows),
                rank_changes_vs_whole_region_FE=sum(r["rank_overall"] != r["alt_rank_whole_region_FE"] for r in rows),
                top10=[(r["rank_overall"], r["region_id"], r["gene_id"], fmt(r["ev_selective_bin_max_FE"]), fmt(r["F5_mean_EV_reads_in_region"]),
                        fmt(r["region_FE_EV_over_Cell"]), r["F4_V3_product"], r["F4_eggNOG_Preferred_name"]) for r in rows[:10]])
    json.dump(summ, open(f"{OUT}/summary.json", "w"), indent=1)
    with open(f"{OUT}/gates.tsv", "w") as w:
        w.write("gate\tresult\tdetail\n")
        for g in gates: w.write("\t".join(g) + "\n")
    print(json.dumps(summ, indent=1))
    ok = all(g[1] == "PASS" for g in gates)
    print("S8_RUN_GATES_" + ("PASS" if ok else "FAIL"))
    return ok, rows


# ------------------------------------------------------------------ known-answer self-test
def selftest(W):
    os.makedirs(W, exist_ok=False)
    for d in ("N5", "E1", "N4E", "N4C"): os.makedirs(f"{W}/{d}")
    chrom = "chrT"
    gseq = "".join("ACGT"[(i * 7 + i // 3) % 4] for i in range(2000))
    open(f"{W}/g.fa", "w").write(f">{chrom} test\n" + "\n".join(gseq[i:i + 60] for i in range(0, 2000, 60)) + "\n")
    # T1: + strand, 5 bins, genome [100,600) one exon. T2: - strand spliced, exons [1000,1150) and [1300,1450): 300 nt, 3 bins
    t2 = {1: [(1350, 1450)], 2: [(1300, 1350), (1100, 1150)], 3: [(1000, 1100)]}  # 5'->3' blocks per bin (minus strand)
    annot = [["bin_id", "rna_class", "model_id", "gene_id", "chrom", "strand", "n_bins", "bin_number", "tx_start0", "tx_end",
              "genomic_blocks_5to3", "rdna_flag"]]
    for i in range(5):
        annot.append([f"T1-T|bin{i+1:03d}", "PCG_CDS", "T1-T", "T1", chrom, "+", 5, i + 1, 100 * i, 100 * i + 100, f"{200+100*i-100}-{200+100*i}", ""])
    for i in (1, 2, 3):
        annot.append([f"T2-T|bin{i:03d}", "lncRNA", "T2-T", "T2", chrom, "-", 3, i, 100 * (i - 1), 100 * i, ";".join(f"{a}-{b}" for a, b in t2[i]), ""])
    open(f"{W}/annot.tsv", "w").write("\n".join("\t".join(map(str, r)) for r in annot) + "\n")

    # planted per group: hotspot (H) / selective category (C) per bin; counts (cell libs then ev libs), transcript reads
    # G0: T1 bins 1,2 hotspot, bin 2 EV-selective -> region 1-2 PASS; T1 bin 4 hotspot Cell_share_higher -> FAIL
    #     T2 bins 2,3 hotspot; bin 2 not significant, bin 3 EV_share_higher (also_Cell) -> spliced region 2-3 PASS
    #     Cell library 1 has 0 T1 reads (G1); T2 bin 3 Cell share 0 -> selective-bin FE Inf
    # G1: T1 bin 2 hotspot, not significant -> FAIL (and recurrence for G0's region)
    # G2: T2 bin 3 hotspot + EV-selective, all counts x2 -> bin FE Inf like G0's T2 region but more EV reads (B2 tie-break)
    plan = {
        0: {("T1-T", 1): ("TRUE", "EV_hotspot__no_Cell_hotspot__no_significant_difference"),
            ("T1-T", 2): ("TRUE", "EV_hotspot__EV_share_higher"),
            ("T1-T", 4): ("TRUE", "EV_hotspot__Cell_share_higher"),
            ("T2-T", 2): ("TRUE", "EV_hotspot__no_Cell_hotspot__no_significant_difference"),
            ("T2-T", 3): ("TRUE", "EV_hotspot__also_Cell_hotspot__EV_share_higher")},
        1: {("T1-T", 2): ("TRUE", "EV_hotspot__no_Cell_hotspot__no_significant_difference")},
        2: {("T2-T", 3): ("TRUE", "EV_hotspot__no_Cell_hotspot__EV_share_higher")}, 3: {}}
    # counts: T1 transcript reads cell (0,10,10) ev (20,20,20); bin1 cell (0,1,1) ev (4,4,4); bin2 cell (0,1,1) ev (8,6,4)
    # T1 tie +1: bin1 and bin2 share one tie read in EV_1 (bin1 +1, bin2 +1; both in region) -> strict EV_1 bin1 3 bin2 7
    # T2 transcript reads cell (5,5,5) ev (10,10,10); bin3 cell 0 ev (5,5,5)
    cnt = {("T1-T", 1): ([0, 1, 1], [4, 4, 4]), ("T1-T", 2): ([0, 1, 1], [8, 6, 4]), ("T1-T", 3): ([0, 2, 2], [2, 2, 2]),
           ("T1-T", 4): ([0, 4, 4], [5, 5, 5]), ("T1-T", 5): ([0, 2, 2], [1, 3, 5]),
           ("T2-T", 1): ([3, 3, 3], [2, 2, 2]), ("T2-T", 2): ([2, 2, 2], [3, 3, 3]), ("T2-T", 3): ([0, 0, 0], [5, 5, 5])}
    tot = {"T1-T": ([0, 10, 10], [20, 20, 20]), "T2-T": ([5, 5, 5], [10, 10, 10])}
    tiep = {("S0_EV_1", "T1-T|bin001"): 1, ("S0_EV_1", "T1-T|bin002"): 1}
    n4rows = collections.defaultdict(list)
    regions5a = [["group", "list", "transcript_id", "first_bin", "last_bin"]]
    for gi, G in enumerate(GROUPS):
        cell = [f"S{gi}_sRNA_{i}" for i in (1, 2, 3)]; ev = [f"S{gi}_EV_{i}" for i in (1, 2, 3)]
        e1h = (["group", "transcript_id", "gene_id", "rna_class", "bin_id", "bin_number", "n_bins", "bin_tx_start0", "bin_tx_end"]
               + ["count_" + l for l in cell + ev] + ["transcript_reads_" + l for l in cell + ev]
               + ["mean_share_EV", "mean_share_Cell", "fold_enrichment_EV_over_Cell", "q_dataset", "EV_hotspot_q0.05", "category_q0.05"])
        E, N = [e1h], [["transcript_id", "rna_class", "bin_id", "bin_number", "n_bins", "group_call_2of3_q0.05", "n_replicates_tested"]]
        for (t, b), (cc, ec) in cnt.items():
            ct, et = tot[t]; m = 2 if gi == 2 else 1
            cc, ec, ct, et = [x * m for x in cc], [x * m for x in ec], [x * m for x in ct], [x * m for x in et]; nb = 5 if t == "T1-T" else 3; cls = "PCG_CDS" if t == "T1-T" else "lncRNA"
            h, cat = plan[gi].get((t, b), ("FALSE", "no_hotspot__no_difference"))
            evm, _ = mean_share(ec, et); cem, _ = mean_share(cc, ct)
            fe = ratio(evm, cem)
            E.append([G, t, t[:2], cls, f"{t}|bin{b:03d}", b, nb, 100 * (b - 1), 100 * b] + cc + ec + ct + et
                     + [full(evm), full(cem), full(fe), "0.01", h, cat])
            N.append([t, cls, f"{t}|bin{b:03d}", b, nb, h, 3])
            for l, c, tt in zip(cell + ev, cc + ec, ct + et):
                tp = tiep.get((l, f"{t}|bin{b:03d}"), 0)
                n4rows[l].append([f"{t}|bin{b:03d}", c - tp, tp, c, tt])
        open(f"{W}/E1/{G}.ev_vs_cell_bins.tsv", "w").write("\n".join("\t".join(map(str, r)) for r in E) + "\n")
        open(f"{W}/N5/{G}.group_calls.tsv", "w").write("\n".join("\t".join(map(str, r)) for r in N) + "\n")
        hot = collections.defaultdict(list); sel = collections.defaultdict(list)
        for (t, b), (h, cat) in plan[gi].items():
            if h == "TRUE": hot[t].append(b)
            if h == "TRUE" and cat.endswith("EV_share_higher"): sel[t].append(b)
        for L, dd in (("EV_hotspot", hot), ("EV_selective", sel)):
            for t, bl in dd.items():
                for run in merge_runs(bl): regions5a.append([G, L, t, run[0], run[-1]])
    for f in os.listdir(f"{W}/E1"): os.system(f"gzip '{W}/E1/{f}'")
    for f in os.listdir(f"{W}/N5"): os.system(f"gzip '{W}/N5/{f}'")
    for l, rr in n4rows.items():
        d = "N4C" if "_sRNA" in l else "N4E"
        with gzip.open(f"{W}/{d}/{l}.bin100_final_table.tsv.gz", "wt") as w:
            w.write("bin_id\tbin_reads_strict_greatest_overlap\tbin_tie_plus1\tbin_count\ttranscript_physical_reads\n")
            for r in rr: w.write("\t".join(map(str, r)) + "\n")
    open(f"{W}/s5a.tsv", "w").write("\n".join("\t".join(map(str, r)) for r in regions5a) + "\n")
    s5b = [["group", "list", "gene_id"] + ANN_COLS]
    for G in GROUPS:
        for g in ("T1", "T2"): s5b.append([G, "EV_hotspot", g] + [f"{g}_{c}" for c in ANN_COLS])
    open(f"{W}/s5b.tsv", "w").write("\n".join("\t".join(r) for r in s5b) + "\n")

    ok, rows = main(f"{W}/N5", f"{W}/E1", f"{W}/N4E", f"{W}/N4C", f"{W}/annot.tsv", f"{W}/g.fa", f"{W}/s5a.tsv", f"{W}/s5b.tsv",
                    f"{W}/out", check_5a_counts=False)
    R = {r["region_id"]: r for r in rows}
    exp_ids = {"B8441_untreated|T1-T|bin001-002", "B8441_untreated|T2-T|bin002-003", "MMC1_untreated|T2-T|bin003-003"}
    good = ok
    def chk(name, cond):
        nonlocal good
        print(f"SELFTEST {name}: {'ok' if cond else 'FAIL'}"); good &= bool(cond)
    chk("exactly the planted passing regions (hotspot-only, Cell-higher, other-group regions fail)", set(R) == exp_ids)
    r1 = R.get("B8441_untreated|T1-T|bin001-002"); r2 = R.get("B8441_untreated|T2-T|bin002-003")
    if r1 and r2:
        # region T1 1-2: EV shares 12/20, 10/20, 8/20 -> mean 0.5; Cell (lib1 excluded, G1) 2/10, 2/10 -> 0.2; FE 2.5
        chk("region share EV (summed bins) = 0.5", abs(r1["region_mean_share_EV"] - 0.5) < 1e-12)
        chk("G1: Cell share over 2 libraries = 0.2", abs(r1["region_mean_share_Cell"] - 0.2) < 1e-12 and r1["n_Cell_proportions"] == 2)
        chk("region FE = 2.5", abs(r1["region_FE_EV_over_Cell"] - 2.5) < 1e-12)
        # strict: EV_1 12-2=10 -> (10/20+10/20+8/20)/3 = 0.466667; Cell 0.2 -> 2.33333
        chk("strict FE (internal tie removed) = 2.333333", abs(r1["region_FE_strict_no_tie_reads"] - 7 / 3) < 1e-9)
        chk("EV-selective bins listed = 2", r1["ev_selective_bins"] == "2" and r1["n_ev_selective_bins"] == 1)
        chk("F3 recurrence = 2 groups (B8441_untreated;B8441_caspofungin)",
            r1["F3_n_groups_with_EV_hotspot_region"] == 2 and r1["F3_n_groups_passing"] == 1)
        chk("plus-strand sequence = genome[100:300]", r1["sequence"] == gseq[100:300])
        r3 = R.get("MMC1_untreated|T2-T|bin003-003")
        # A3: rank by the EV-selective bin. T1 bin 2: EV (8,6,4)/20 = 0.3, Cell (lib1 G1-excluded) 1/10, 1/10 = 0.1 -> 3
        chk("T1 selective-bin FE = 3", abs(r1["ev_selective_bin_max_FE"] - 3) < 1e-12)
        # whole-region ratios kept as columns: T2 2-3 EV 0.8 / Cell 0.4 = 2; alt rank by whole region puts T1 (2.5) above T2 (2)
        chk("T2 whole-region FE = 2 (column)", abs(r2["region_FE_EV_over_Cell"] - 2) < 1e-12 and r1["alt_rank_whole_region_FE"] < r2["alt_rank_whole_region_FE"])
        # B2: two Inf regions; G2 (MMC1_untreated) has 10 EV reads vs G0's 8 -> ranks 1 and 2 although its ID sorts later; T1 rank 3
        chk("Inf tie broken by EV depth (MMC1 10 reads > B8441 8 reads)", r3 is not None and math.isinf(r3["ev_selective_bin_max_FE"])
            and r3["rank_overall"] == 1 and r2["rank_overall"] == 2 and r1["rank_overall"] == 3 and r3["F5_mean_EV_reads_in_region"] == 10)
        chk("rank_in_group restarts per group", r3["rank_in_group"] == 1 and r2["rank_in_group"] == 1 and r1["rank_in_group"] == 2)
        rc = lambda s: s.translate(COMP)[::-1]
        chk("minus-strand spliced sequence = rc(g[1300:1350]) + rc(g[1000:1150])", r2["sequence"] == rc(gseq[1300:1350]) + rc(gseq[1000:1150])
            and r2["genomic_blocks_5to3"] == "1300-1350;1000-1150" and r2["region_length"] == 200)
        chk("F5 depth EV region reads = 8", r2["F5_mean_EV_reads_in_region"] == 8)
        chk("ratio(): Cell 0 -> Inf; both 0 -> NA", math.isinf(ratio(0.3, 0)) and math.isnan(ratio(0, 0)))
        chk("F4 annotation joined", r1["F4_V3_product"] == "T1_V3_product")
    # a separate check of the minus-strand spliced path: bins 1-2 of T2 = revcomp(1350-1450) + revcomp(1300-1350) + revcomp(1100-1150)
    seq = seq_of({"c": gseq}, "c", "-", [(1300, 1450), (1100, 1150)])
    chk("spliced minus-strand join (merged blocks)", seq == (gseq[1300:1450].translate(COMP)[::-1] + gseq[1100:1150].translate(COMP)[::-1]))
    chk("merge_runs", merge_runs([2, 3, 4, 7, 9]) == [[2, 3, 4], [7], [9]])
    print("S8_SELFTEST_" + ("PASS" if good else "FAIL"))
    return good


if __name__ == "__main__":
    if sys.argv[1] == "selftest":
        sys.exit(0 if selftest(sys.argv[2]) else 1)
    elif sys.argv[1] == "run":
        ok, _ = main(*sys.argv[2:11])
        sys.exit(0 if ok else 1)
