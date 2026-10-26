#!/usr/bin/env python3
"""tpm_abundance_dump.py -- 2026-10-12, READ-ONLY derivative of fl_classify_r3.py: identical inputs and TPM (with rRNA)
computation; additionally writes per-library TPM per unit plus the rDNA-array unit (census rDNA reads / 8,130 nt), for the
abundance-weighted pies (Justine 2026-10-12: TPM, cells = Cell sRNA). Classification outputs are written too but not used.
--- fl_classify_r3 docstring ---
fl_classify_r3.py -- FL r3 (2026-10-09, v2_rRNA branch). Derived from frozen fl_classify.py r2 (md5 52260e9a...).
r3 CHANGES ONLY (Justine 2026-10-09: all RNAs with meaningful biological labels; ncRNA rule = PCG rule; no 18S/5.8S/25S):
  + 26 PCG units (GTF transcripts without CDS lines; NCBI V3 = mRNA; V3 CDS frame); reads and 100-nt bin aligned bases from
    fl_extra_count.py (same read rules as the r2 inputs, verified by gate_extra_counts.py);
  + one pooled 5S rRNA unit (13 V3 loci, 118 nt; multimappers pooled like a tRNA family); coverage = PCG rule on its unmerged
    100-nt bins (0-100, 100-118); reads = 5S reads from fl_extra_count.py;
  TPM: without_rRNA denominator = PCG + lncRNA + tRNA families (incl. the 26); the 5S unit's without_rRNA TPM = its RPK /
    (that denominator + its own RPK). with_rRNA denominator = all units incl. 5S + the census rDNA-array unit (the census 5S
    unit of r2 is replaced by the 5S unit, so 5S is not counted twice).
Extra usage args (after the r2 args, Cell_polyA_counts_dir required): <flextra_dir> <new_bins100.tsv>
--- r2 docstring follows ---
fl_classify.py -- full-length vs fragment-only classification (2026-10-09; criteria registered 2026-10-09, third entry).
Usage: fl_classify.py <T1_counts_dir> <manifest (header)> <bin100_master_annotation.tsv (unmerged)> <trna168.bed> <census_dir>
                      <EV_tables_dir (B4)> <Cell_tables_dir (C4)> <outdir> [Cell_polyA_counts_dir (v2_bin150 R2_counts)]
r2 (2026-10-09): Cell poly(A) libraries supported; their transcript physical reads come from the v2_bin150 R2 counter
(transcript_totals: physical_reads; NH = 1, >= 18 nt, strand-flipped), their bin coverage from trna_count.py run with kind polyA.
Units: protein-coding CDS and lncRNA transcripts (not the rDNA-associated LNCRNA.1946/.1057, which are rRNA), and tRNA anticodon
families (tRNAscan 168 / 43).
Per library:
  reads  - PCG/lncRNA: physical reads per transcript (NH = 1, >= 18 nt; B4/C4 tables); tRNA: family-assigned reads (T1);
           rRNA (only for the 'with rRNA' TPM): read-census reads in the rDNA array (one unit, LNCRNA.1946 + .1057 model length
           8,130 nt) and in 5S rRNA (one unit, 118 nt).
  TPM    - RPK = reads / (length in kb); TPM = RPK / sum(RPK) x 1e6, computed twice: without rRNA (PCG + lncRNA + tRNA families)
           and with rRNA (plus the two rRNA units). Length: transcript length; tRNA family = median mature length of its loci.
  cover  - PCG/lncRNA: every unmerged 100-nt bin has mean depth >= 3x (aligned bases / bin length; final bin at its true length);
           tRNA family: mean depth over the mature tRNA >= 3x (aligned bases / family length).
Per strain x treatment group (3 replicates):
  detected     - >= 1 read in >= 2 of 3 replicates;
  full-length  - TPM > 100 AND coverage criterion in >= 2 of 3 replicates; tRNA families additionally need >= 1 CCA read in
                 the group's libraries (read-level CCA evidence, v1 Step 6C follow-up (1) method).
Per compartment: detected in >= 1 group; full-length if full-length in >= 1 group; otherwise 'fragment-only'."""
import sys, os, json, csv, gzip, collections
T1, MAN, ANN, TB, CEN, EVT, CT, OUT, PAT, FXD, NEWB = sys.argv[1:12]; os.makedirs(OUT, exist_ok=True)
RDNA_LEN, R5S_LEN = 8130, 118
man = list(csv.DictReader(open(MAN), delimiter="\t"))
bins = collections.defaultdict(list); info = {}
for r in csv.DictReader(open(ANN), delimiter="\t"):
    if r["role"] == "competitor_only": continue
    bins[r["model_id"]].append((r["bin_id"], int(r["bin_length"]))); info[r["model_id"]] = (r["rna_class"], int(r["transcript_length"]))
ml = collections.defaultdict(list)
for line in open(TB):
    if line.startswith("#"): continue
    x = line.rstrip("\n").split("\t"); s, e = int(x[1]), int(x[2]); i = x[8]
    ml[x[6]].append((e - s) - (0 if i == "NA" else int(i.split("-")[1]) - int(i.split("-")[0])))
famlen = {F: sorted(v)[(len(v) - 1) // 2] for F, v in ml.items()}
for F in famlen: info[F] = ("tRNA_family", famlen[F])
NEW = set()
for r in csv.DictReader(open(NEWB), delimiter="\t"):
    bins[r["model_id"]].append((r["bin_id"], int(r["bin_length"]))); info[r["model_id"]] = ("PCG_CDS", int(r["transcript_length"])); NEW.add(r["model_id"])
assert len(NEW) == 26
U5 = "5S_rRNA"; info[U5] = ("rRNA_5S", R5S_LEN); bins[U5] = [("5S_rRNA|bin001", 100), ("5S_rRNA|bin002", R5S_LEN - 100)]
units = sorted(info)
per = {}; gates = []; ok = True; DUMP = []
for r in man:
    L = r["library"]; o = json.load(open(f"{T1}/{L}.trna.json"))
    txr = {}
    if r["compartment"] == "Cell_polyA":
        for x in csv.DictReader(gzip.open(f"{PAT}/{L}.transcript_totals.tsv.gz", "rt"), delimiter="\t"): txr[x["model_id"]] = int(x["physical_reads"])
    else:
        tabdir = EVT if r["compartment"] == "EV" else CT
        for x in csv.DictReader(gzip.open(f"{tabdir}/{L}.bin100_final_table.tsv.gz", "rt"), delimiter="\t"): txr.setdefault(x["transcript_id"], int(x["transcript_physical_reads"]))
    cen = json.load(open(f"{CEN}/{L}.census.json"))["primary_class"]
    fx = json.load(open(f"{FXD}/{L}.flextra.json")); assert fx["library"] == L
    rd = {u: (o["trna_families"].get(u, {}).get("reads", 0) if info[u][0] == "tRNA_family" else txr.get(u, 0)) for u in units}
    for u in NEW: rd[u] = fx["physical_reads"].get(u, 0)
    rd[U5] = fx["rRNA_5S"]["reads"]
    rpk = {u: rd[u] / (info[u][1] / 1000) for u in units}
    rr = {"rDNA": sum(v for k, v in cen.items() if k.startswith("rDNA_lncRNA|")) / (RDNA_LEN / 1000)}
    s0 = sum(v for u, v in rpk.items() if u != U5); s1 = s0 + rpk[U5] + sum(rr.values())
    tpm0 = {u: 1e6 * rpk[u] / s0 for u in units if u != U5}; tpm0[U5] = 1e6 * rpk[U5] / (s0 + rpk[U5])
    tpm1 = {u: 1e6 * rpk[u] / s1 for u in units}
    bb = dict(o["bin_aligned_bases"]); bb.update({k: v for k, v in fx["bin_aligned_bases"].items() if k.split("|")[0] in NEW})
    bb.update(fx["rRNA_5S"]["bin_aligned_bases"]); cov = {}; mind = {}
    for u in units:
        if info[u][0] == "tRNA_family":
            d = o["trna_families"].get(u, {}).get("bases", 0) / famlen[u]; cov[u] = d >= 3; mind[u] = d
        else:
            ds = [bb.get(b, 0) / bl for b, bl in bins[u]]; cov[u] = min(ds) >= 3; mind[u] = min(ds)
    cca = {u: o["trna_families"].get(u, {}).get("CCA_reads", 0) for u in units if info[u][0] == "tRNA_family"}
    per[L] = dict(r=r, rd=rd, tpm0=tpm0, tpm1=tpm1, cov=cov, mind=mind, cca=cca)
    DUMP.append((L, r["compartment"], {**tpm1, "rDNA_array_rRNA": 1e6 * rr["rDNA"] / s1}))
    g1 = abs(sum(v for u, v in tpm0.items() if u != U5) - 1e6) < 1e-3 and abs(sum(tpm1.values()) + 1e6 * sum(rr.values()) / s1 - 1e6) < 1e-3
    g2 = all(u in txr or rd[u] == 0 for u in units if info[u][0] not in ("tRNA_family", "rRNA_5S") and u not in NEW) and not (set(txr) & NEW) and set(txr) - set(units) <= {"LNCRNA.1946", "LNCRNA.1057"}
    if r["compartment"] == "Cell_polyA": g2 = set(txr) - set(units) <= {"LNCRNA.1946", "LNCRNA.1057"}   # R2 lists only transcripts with reads
    for gname, v in (("F1_TPM_sums_to_1e6_both_versions", g1), ("F2_transcripts_match_tables", g2)):
        gates.append((L, gname, "PASS" if v else "FAIL")); ok &= v
groups = collections.defaultdict(list)
for L, p in per.items(): groups[(p["r"]["compartment"], p["r"]["group"])].append(L)
rows = []; summ = collections.Counter(); gsum = collections.Counter()
for comp in sorted({c for c, _ in groups}):
    for u in units:
        det_g, fl0_g, fl1_g = [], [], []
        for (c, g), libs in sorted(groups.items()):
            if c != comp: continue
            assert len(libs) == 3
            det = sum(per[L]["rd"][u] >= 1 for L in libs) >= 2
            ccaok = True if info[u][0] != "tRNA_family" else sum(per[L]["cca"][u] for L in libs) >= 1
            f0 = sum(per[L]["tpm0"][u] > 100 and per[L]["cov"][u] for L in libs) >= 2 and ccaok
            f1 = sum(per[L]["tpm1"][u] > 100 and per[L]["cov"][u] for L in libs) >= 2 and ccaok
            if det: det_g.append(g)
            if f0: fl0_g.append(g)
            if f1: fl1_g.append(g)
            for v, f in (("without_rRNA", f0), ("with_rRNA", f1)):
                if det: gsum[(comp, g, v, info[u][0], "full_length" if f else "fragment_only")] += 1
        cls = info[u][0]
        if det_g or fl0_g or fl1_g:
            for v, fl in (("without_rRNA", fl0_g), ("with_rRNA", fl1_g)):
                summ[(comp, v, cls, "full_length" if fl else "fragment_only")] += 1
        libs_c = [L for (c, g), ls in groups.items() if c == comp for L in ls]
        rows.append([comp, u, cls, info[u][1], ",".join(det_g), ",".join(fl0_g), ",".join(fl1_g),
                     round(max(per[L]["tpm0"][u] for L in libs_c), 2), round(max(per[L]["tpm1"][u] for L in libs_c), 2),
                     round(max(per[L]["mind"][u] for L in libs_c), 3), sum(per[L]["cca"][u] for L in libs_c) if cls == "tRNA_family" else "NA"])
def Wt(name, hdr, data, gz=False):
    f = gzip.open(f"{OUT}/{name}", "wt") if gz else open(f"{OUT}/{name}", "w")
    f.write("\t".join(hdr) + "\n"); [f.write("\t".join(map(str, x)) + "\n") for x in data]; f.close()
Wt("unit_table.tsv.gz", ["compartment", "unit", "class", "length", "detected_groups", "full_length_groups_TPM_without_rRNA",
                         "full_length_groups_TPM_with_rRNA", "max_TPM_without_rRNA", "max_TPM_with_rRNA",
                         "max_min_bin_or_mean_depth", "CCA_reads_compartment"], rows, gz=True)
Wt("summary.tsv", ["compartment", "TPM_version", "class", "status", "units"], [[*k, v] for k, v in sorted(summ.items())])
Wt("group_summary.tsv", ["compartment", "group", "TPM_version", "class", "status", "units"], [[*k, v] for k, v in sorted(gsum.items())])
Wt("gates.tsv", ["library", "gate", "result"], gates)
with gzip.open(f"{OUT}/tpm_with_rRNA_per_library.tsv.gz", "wt") as w:
    w.write("library\tcompartment\tunit\tclass\tTPM_with_rRNA\n")
    for L, c, t in DUMP:
        assert abs(sum(t.values()) - 1e6) < 1e-3
        for u, v in sorted(t.items()): w.write(f"{L}\t{c}\t{u}\t{info[u][0] if u in info else 'rRNA_rDNA_array'}\t{v}\n")
print(open(f"{OUT}/summary.tsv").read()); print("FL_ALL_GATES_PASS" if ok else "FL_GATE_FAILURE"); sys.exit(0 if ok else 1)
