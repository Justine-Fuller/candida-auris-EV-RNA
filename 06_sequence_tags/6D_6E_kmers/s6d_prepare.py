#!/usr/bin/env python3
"""s6d_prepare.py -- Step 6D (v2, 2026-10-12): inputs for the annotation table and the enrichment of the FL r3 EV full-length PCGs.
Approved (Justine 2026-10-12): Q1 a = enrichment for B8441 caspofungin only (70 genes; MMC1 caspofungin's 4 are listed, not
tested); Q2 a = background = protein-coding genes detected in B8441 caspofungin EVs (>= 1 physical read in >= 2 of the 3 EV
libraries); method = s5b_enrich.R r5 (standard + read-depth-corrected), BH q < 0.05.
READ-ONLY inputs. Usage: s6d_prepare.py <FL r3 unit_table.tsv.gz> <tpm_with_rRNA_per_library.tsv.gz> <bin100 B4_tie_rule dir>
                                         <FL3_counts dir> <PCG_functional_annotation.tsv> <V3 GFF.gz> <outdir>
Reads per gene and library: transcript physical reads (NH = 1, >= 18 nt) from the bin100 B4 EV tables (ordinary models), and
fl_extra_count.py physical_reads for the models absent from B4 (the 26 FL r3 CDS models) (the same inputs FL r3 used).
Covariate for the read-depth-corrected test (s5b_enrich.R reads it from the column named bins_tested_ge2rep; the script is used
unchanged): two versions are written, side by side --
  prepare_TPM/   mean EV TPM (with rRNA, as FL r3) over the 3 B8441 caspofungin EV libraries;
  prepare_reads/ mean EV physical reads over the same libraries.
Outputs: prepare_*/enrichment_universe.tsv, term2gene_GO.tsv, term2gene_KEGG.tsv (GO / KEGG exactly as s5b_prepare.py r2);
fl_genes.tsv (the EV full-length PCGs with groups, reads, TPM and annotation columns incl. eggNOG GO ids)."""
import sys, os, gzip, csv, re, json, collections
UT, TPM, B4, FX, INT, V3, OUT = sys.argv[1:8]; os.makedirs(OUT, exist_ok=False)
G = "B8441_caspofungin"; LIBS = ["B8441_EV_CASP_1", "B8441_EV_CASP_2", "B8441_EV_CASP_3"]
SPLIT = lambda v: [t.strip() for t in re.split(r"[,;]", v or "") if t.strip()]
gene = lambda t: t[:-2] if t.endswith("-T") else t
fl = {}
for r in csv.DictReader(gzip.open(UT, "rt"), delimiter="\t"):
    if r["compartment"] == "EV" and r["class"] == "PCG_CDS" and r["full_length_groups_TPM_with_rRNA"] not in ("", "NA"):
        fl[r["unit"]] = r["full_length_groups_TPM_with_rRNA"]
reads = collections.defaultdict(dict)
for L in LIBS:
    seen = set(); allb4 = set()
    for r in csv.DictReader(gzip.open(f"{B4}/{L}.bin100_final_table.tsv.gz", "rt"), delimiter="\t"):
        allb4.add(r["transcript_id"])
        if r["rna_class"] != "PCG_CDS" or r["role"] != "ordinary" or r["transcript_id"] in seen: continue
        seen.add(r["transcript_id"]); reads[r["transcript_id"]][L] = int(r["transcript_physical_reads"])
    fx = json.load(open(f"{FX}/{L}.flextra.json"))["physical_reads"]
    for m, n in fx.items():                      # flextra holds all models; only the 26 FL r3 models are not in the B4 tables
        if m in allb4: continue
        reads[m][L] = int(n)
tpm = collections.defaultdict(dict)
for r in csv.DictReader(gzip.open(TPM, "rt"), delimiter="\t"):
    if r["library"] in LIBS and r["class"] == "PCG_CDS": tpm[r["unit"]][r["library"]] = float(r["TPM_with_rRNA"])
units = sorted(reads)
det = [u for u in units if sum(reads[u].get(L, 0) > 0 for L in LIBS) >= 2]
hits = {u for u, gs in fl.items() if G in gs.split(",")}
checks = {"fl_genes_total": len(fl), "fl_genes_B8441_caspofungin": len(hits), "pcg_units_with_reads_table": len(units),
          "universe_detected_ge2of3": len(det), "hits_in_universe": len(hits & set(det)),
          "universe_units_missing_TPM": sum(1 for u in det if len(tpm.get(u, {})) != 3),
          "tpm_zero_iff_reads_zero": all((tpm[u][L] > 0) == (reads[u].get(L, 0) > 0) for u in det for L in LIBS)}
print(json.dumps(checks)); assert checks["hits_in_universe"] == len(hits) and checks["universe_units_missing_TPM"] == 0 and checks["tpm_zero_iff_reads_zero"]
ann = {r["gene_id"]: r for r in csv.DictReader(open(INT), delimiter="\t")}
gos = {g: [t for t in SPLIT(a.get("eggNOG_GOs")) if t.startswith("GO:")] for g, a in ann.items()}
kegg = {g: sorted({t for t in SPLIT(a.get("eggNOG_KEGG_Pathway")) if t.startswith("ko")}) for g, a in ann.items()}
assert all(re.fullmatch(r"GO:\d{7}", t) for ts in gos.values() for t in ts) and all(re.fullmatch(r"ko\d{5}", t) for ts in kegg.values() for t in ts)
for name, cov in (("prepare_TPM", lambda u: sum(tpm[u][L] for L in LIBS) / 3), ("prepare_reads", lambda u: sum(reads[u].get(L, 0) for L in LIBS) / 3)):
    os.makedirs(f"{OUT}/{name}")
    with open(f"{OUT}/{name}/enrichment_universe.tsv", "w") as w:
        w.write("group\tgene_id\tbins_tested_ge2rep\thit\n")
        for u in det: w.write(f"{G}\t{gene(u)}\t{cov(u)!r}\t{int(u in hits)}\n")
    with open(f"{OUT}/{name}/term2gene_GO.tsv", "w") as w:
        w.write("term\tgene_id\n"); [w.write(f"{t}\t{g}\n") for g, ts in sorted(gos.items()) for t in sorted(set(ts))]
    with open(f"{OUT}/{name}/term2gene_KEGG.tsv", "w") as w:
        w.write("term\tgene_id\n"); [w.write(f"{t}\t{g}\n") for g, ts in sorted(kegg.items()) for t in ts]
v3 = {}
with gzip.open(V3, "rt") as f:
    for l in f:
        x = l.rstrip("\n").split("\t")
        if len(x) < 9 or x[2] in ("gene", "exon", "CDS", "region"): continue
        p = re.search(r"product=([^;]+)", x[8]); dec = lambda s: s.replace("%2C", ",").replace("%3B", ";").replace("%20", " ")
        if not p: continue
        m = re.search(r"B8441v2:([A-Za-z0-9_]+?)T0", x[8])
        if m: v3[m.group(1)] = dec(p.group(1))
        lt = re.search(r"locus_tag=([A-Za-z0-9_]+)", x[8])
        if lt: v3.setdefault(lt.group(1), dec(p.group(1)))
AC = ["eggNOG_Preferred_name", "eggNOG_Description", "eggNOG_COG_category", "eggNOG_KEGG_Pathway", "Candida_albicans_orthologs",
      "Saccharomyces_cerevisiae_orthologs", "SignalP_prediction", "DeepTMHMM_has_TM", "PredGPI_prediction", "eggNOG_GOs"]
with open(f"{OUT}/fl_genes.tsv", "w") as w:
    w.write("\t".join(["gene_id", "full_length_groups", "in_integrated_annotation", "V3_product"] + [f"reads_{L}" for L in LIBS] +
                      [f"TPM_{L}" for L in LIBS] + ["mean_TPM_B8441_caspofungin_EV"] + AC) + "\n")
    for u in sorted(fl):
        g = gene(u); a = ann.get(g, {})
        w.write("\t".join([g, fl[u], str(g in ann), v3.get(g, "")] + [str(reads[u].get(L, 0)) for L in LIBS] +
                          [f"{tpm.get(u, {}).get(L, 0):.2f}" for L in LIBS] + [f"{sum(tpm.get(u, {}).get(L, 0) for L in LIBS) / 3:.2f}"] +
                          [a.get(c, "") for c in AC]) + "\n")
print("fl genes without integrated annotation:", [gene(u) for u in fl if gene(u) not in ann], "without V3 product:", [gene(u) for u in fl if gene(u) not in v3])
print("S6D_PREPARE_DONE")
