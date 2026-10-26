# v2 of build_strain_DE_annotation.py (2026-10-14): same DE join, plus whole-proteome (WP_) InterProScan,
# SignalP, DeepTMHMM and PredGPI columns from the 2026-10-14 whole-proteome run (all 5,417 proteins,
# same tool settings as the old up/down runs). The old up/down-only columns are kept unchanged.
# Also writes the old-vs-new agreement check for genes that were in the old up/down runs.
# Raw tool outputs are parsed with the readers in build_pcg_annotation_table.py (the script that built
# the v1 annotation table), so the WP columns are built the same way as the old ones.
#
# usage: python3 build_strain_DE_annotation_v2.py <build_pcg_annotation_table.py> <DE csv> <v1 annotation tsv>
#                <whole-proteome results dir> <out dir>
import sys, importlib.util, logging, glob, os
import pandas as pd

BUILDER, EXP, ANN, WP, OUT = sys.argv[1:6]
spec = importlib.util.spec_from_file_location("bpat", BUILDER); B = importlib.util.module_from_spec(spec)
spec.loader.exec_module(B)
log = logging.getLogger("wp"); logging.basicConfig(level=logging.WARNING)
J = os.path.join

# ---------- whole-proteome tool outputs -> one row per gene ----------
tested = []
for line in open(J(WP, "input/proteome_5417.fasta")):
    if line.startswith(">"): tested.append(B.normalize_gene_id(line[1:].split()[0]))
assert len(tested) == len(set(tested)) == 5417
wp = pd.DataFrame({"gene_id": tested})

ipr = B._read_interpro_tsv(B.Path(J(WP, "InterProScan/wp_interproscan.tsv")), log)
ipr["gene_id"] = B.normalize_id_series(ipr.protein_accession.tolist(), "InterProScan", log)
ipr = ipr.groupby("gene_id", sort=False).apply(lambda g: pd.Series({
    "WP_InterPro_accessions": B.unique_join(g.interpro_accession),
    "WP_InterPro_descriptions": B.unique_join(g.interpro_description),
    "WP_InterPro_databases": B.unique_join(g.analysis),
    "WP_InterPro_signature_accessions": B.unique_join(g.signature_accession),
    "WP_InterPro_signature_descriptions": B.unique_join(g.signature_description),
    "WP_InterPro_GOs": B.split_and_collect(g.go_terms, split_on="|"),
    "WP_InterPro_pathways": B.split_and_collect(g.pathways, split_on="|")})).reset_index()

sp = B._read_signalp_file(B.Path(J(WP, "SignalP/prediction_results.txt")), log)
sp["gene_id"] = B.normalize_id_series(sp.protein_accession.tolist(), "SignalP", log)
assert not sp.gene_id.duplicated().any()
sp = sp.rename(columns={"prediction": "WP_SignalP_prediction", "prob_sp": "WP_SignalP_SP_probability",
                        "cs_position": "WP_SignalP_cleavage_site"})[
    ["gene_id", "WP_SignalP_prediction", "WP_SignalP_SP_probability", "WP_SignalP_cleavage_site"]]

chunks = sorted(glob.glob(J(WP, "DeepTMHMM/chunk*/biolib_results")))
# 55 chunks of <=99 proteins (biolib anonymous size limit). If some chunks are still missing, genes in them are
# marked WP_DeepTMHMM_tested = False (interim run); the final run must have all 55.
print(f"DeepTMHMM chunks with results: {len(chunks)} / 55" + ("" if len(chunks) == 55 else "  <-- INTERIM, incomplete"))
gff = pd.concat([B._parse_tmrs_gff3(B.Path(J(c, "TMRs.gff3")), log) for c in chunks], ignore_index=True)
cls = pd.concat([B._parse_3line_classes(B.Path(J(c, "predicted_topologies.3line")), log) for c in chunks],
                ignore_index=True)
tm = gff.merge(cls, on="protein_accession", how="outer", validate="1:1")
assert tm.DeepTMHMM_topology_class.notna().all() and tm.DeepTMHMM_TM_count.notna().all()
assert len(chunks) < 55 or len(tm) == 5417, len(tm)
tm["gene_id"] = B.normalize_id_series(tm.protein_accession.tolist(), "DeepTMHMM", log)
assert not tm.gene_id.duplicated().any()
tm["WP_DeepTMHMM_has_TM"] = tm.DeepTMHMM_TM_count > 0
tm["DeepTMHMM_TM_count"] = tm.DeepTMHMM_TM_count.apply(B.clean_integer_display)
tm = tm.rename(columns={"DeepTMHMM_TM_count": "WP_DeepTMHMM_TM_count", "DeepTMHMM_TM_regions": "WP_DeepTMHMM_TM_regions",
                        "DeepTMHMM_topology_class": "WP_DeepTMHMM_topology_class"})[
    ["gene_id", "WP_DeepTMHMM_has_TM", "WP_DeepTMHMM_TM_count", "WP_DeepTMHMM_TM_regions", "WP_DeepTMHMM_topology_class"]]

gpi = B._parse_predgpi_json(B.Path(J(WP, "PredGPI/wp_predgpi.json")), log)
gpi["gene_id"] = B.normalize_id_series(gpi.protein_accession.tolist(), "PredGPI", log)
assert not gpi.gene_id.duplicated().any()
gpi["PredGPI_anchor_position"] = gpi.PredGPI_anchor_position.apply(B.clean_integer_display)
gpi = gpi.rename(columns={c: "WP_" + c for c in ["PredGPI_prediction", "PredGPI_anchor_position", "PredGPI_score"]})[
    ["gene_id", "WP_PredGPI_prediction", "WP_PredGPI_anchor_position", "WP_PredGPI_score"]]

for name, df, hit in [("InterProScan", ipr, None), ("SignalP", sp, "WP_SignalP_prediction"),
                      ("DeepTMHMM", tm, None), ("PredGPI", gpi, None)]:
    wp[f"WP_{name}_tested"] = wp.gene_id.isin(df.gene_id) if name != "InterProScan" else True
    wp = wp.merge(df, on="gene_id", how="left", validate="1:1")
wp["WP_InterPro_has_hit"] = wp.WP_InterPro_accessions.notna() | wp.WP_InterPro_signature_accessions.notna()
wp = wp.rename(columns={"WP_InterProScan_tested": "WP_InterPro_tested"})
wp["WP_GPI_without_signal_peptide"] = wp.WP_PredGPI_prediction.eq("GPI-anchored") & wp.WP_SignalP_prediction.ne("SP")
# GPI rule chosen by Justine (2026-10-14): GPI-anchored only if PredGPI says GPI-anchored AND a signal peptide is
# predicted by either tool (SignalP = SP, or DeepTMHMM topology class SP / SP+TM). Raw PredGPI column kept.
wp["WP_GPI_final"] = wp.WP_PredGPI_prediction.eq("GPI-anchored") & (
    wp.WP_SignalP_prediction.eq("SP") | wp.WP_DeepTMHMM_topology_class.isin(["SP", "SP+TM"]))
cov = pd.DataFrame([
    dict(tool="InterProScan", genes_tested=5417, genes_with_hit=int(wp.WP_InterPro_has_hit.sum()),
         hit_definition="any signature hit (incl. MobiDB-lite, Coils)"),
    dict(tool="InterProScan", genes_tested=5417, genes_with_hit=int(wp.WP_InterPro_accessions.notna().sum()),
         hit_definition="hit with an InterPro (IPR) entry"),
    dict(tool="SignalP", genes_tested=int(wp.WP_SignalP_tested.sum()), genes_with_hit=int(wp.WP_SignalP_prediction.eq("SP").sum()),
         hit_definition="prediction SP"),
    dict(tool="DeepTMHMM", genes_tested=int(wp.WP_DeepTMHMM_tested.sum()), genes_with_hit=int(wp.WP_DeepTMHMM_has_TM.eq(True).sum()),
         hit_definition=">=1 TM region"),
    dict(tool="DeepTMHMM", genes_tested=int(wp.WP_DeepTMHMM_tested.sum()),
         genes_with_hit=int(wp.WP_DeepTMHMM_topology_class.isin(["SP", "SP+TM"]).sum()), hit_definition="topology class SP or SP+TM"),
    dict(tool="PredGPI", genes_tested=int(wp.WP_PredGPI_tested.sum()), genes_with_hit=int(wp.WP_PredGPI_prediction.eq("GPI-anchored").sum()),
         hit_definition="GPI-anchored"),
    dict(tool="GPI_final", genes_tested=5417, genes_with_hit=int(wp.WP_GPI_final.sum()),
         hit_definition="PredGPI GPI-anchored AND (SignalP SP OR DeepTMHMM SP/SP+TM)"),
    dict(tool="PredGPI+SignalP", genes_tested=5417, genes_with_hit=int(wp.WP_GPI_without_signal_peptide.sum()),
         hit_definition="GPI-anchored but SignalP not SP")])
cov.to_csv(J(OUT, "whole_proteome_coverage.csv"), index=False)
wp.to_csv(J(OUT, "whole_proteome_annotation_WP.tsv"), sep="\t", index=False)

# ---------- agreement with the old up/down runs ----------
ann = pd.read_csv(ANN, sep="\t", dtype=str)
a = ann.merge(wp, on="gene_id", how="left", validate="1:1")
def s(x): return set() if pd.isna(x) else set(str(x).split(";"))
rows = []
def agree(tool, field, mask, old, new, how="equal"):
    g = a[mask]
    if how == "set":
        eq = [s(o) == s(n) for o, n in zip(g[old], g[new])]
    else:
        eq = [(pd.isna(o) and pd.isna(n)) or (not pd.isna(o) and not pd.isna(n) and str(o) == str(n)) for o, n in zip(g[old], g[new])]
    rows.append(dict(tool=tool, field=field, genes_compared=len(g), agree=sum(eq), pct_agree=round(100 * sum(eq) / max(len(g), 1), 2)))
    return g[[not e for e in eq]][["gene_id", old, new]].assign(tool=tool, field=field)
old_ipr = a.InterPro_tested.eq("True")   # v1: genes with >=1 InterProScan row in the old up/down runs
dis = [agree("InterProScan", "set of InterPro (IPR) accessions", old_ipr, "InterPro_accessions", "WP_InterPro_accessions", "set"),
       agree("InterProScan", "set of member-database signatures", old_ipr, "InterPro_signature_accessions", "WP_InterPro_signature_accessions", "set"),
       agree("InterProScan", "set of InterPro GO terms", old_ipr, "InterPro_GOs", "WP_InterPro_GOs", "set"),
       agree("SignalP", "prediction (SP/OTHER)", a.SignalP_tested.eq("True"), "SignalP_prediction", "WP_SignalP_prediction"),
       agree("SignalP", "cleavage site", a.SignalP_tested.eq("True"), "SignalP_cleavage_site", "WP_SignalP_cleavage_site"),
       agree("DeepTMHMM", "topology class", a.DeepTMHMM_tested.eq("True"), "DeepTMHMM_topology_class", "WP_DeepTMHMM_topology_class"),
       agree("DeepTMHMM", "TM count", a.DeepTMHMM_tested.eq("True"), "DeepTMHMM_TM_count", "WP_DeepTMHMM_TM_count"),
       agree("PredGPI", "prediction (GPI/not)", a.PredGPI_tested.eq("True"), "PredGPI_prediction", "WP_PredGPI_prediction"),
       agree("PredGPI", "anchor position", a.PredGPI_tested.eq("True"), "PredGPI_anchor_position", "WP_PredGPI_anchor_position")]
pd.DataFrame(rows).to_csv(J(OUT, "agreement_old_updown_vs_whole_proteome.csv"), index=False)
pd.concat([d.rename(columns={d.columns[1]: "old", d.columns[2]: "new"}) for d in dis]).to_csv(
    J(OUT, "agreement_disagreements.tsv"), sep="\t", index=False)
print(pd.DataFrame(rows).to_string(index=False))

# ---------- strain DE table (same logic as v1) + WP columns ----------
de = pd.read_csv(EXP); de = de[de.Feature.str.startswith("B9J08_")]
def call(lfc, p):
    return pd.Series(["up" if (pp < 0.05 and l > 0) else "down" if (pp < 0.05 and l < 0) else "ns"
                      for l, pp in zip(lfc, p.fillna(1))], index=lfc.index)
de["B8441_DE"] = call(de.LFC_B8441, de.padj_B8441); de["MMC1_DE"] = call(de.LFC_MMC1, de.padj_MMC1)
def cat(b, m):
    if b == "ns" and m == "ns": return "not DE"
    if b == m: return "shared " + b
    if b == "ns": return "MMC1 only " + m
    if m == "ns": return "B8441 only " + b
    return "opposite (B8441 %s, MMC1 %s)" % (b, m)
de["DE_category"] = [cat(b, m) for b, m in zip(de.B8441_DE, de.MMC1_DE)]
cols = ["Feature", "LFC_B8441", "padj_B8441", "B8441_DE", "LFC_MMC1", "padj_MMC1", "MMC1_DE", "DE_category"]
m = de[cols].merge(ann, left_on="Feature", right_on="gene_id", how="left", indicator=True)
m["in_FungiDB_proteome"] = m._merge.eq("both"); m = m.drop(columns="_merge")
m = m.merge(wp, on="gene_id", how="left", validate="m:1")
for c in ["WP_InterPro_tested", "WP_SignalP_tested", "WP_DeepTMHMM_tested", "WP_PredGPI_tested"]:
    m[c] = m[c].eq(True)
m.to_csv(J(OUT, "PCG_annotation_by_strain.tsv"), sep="\t", index=False)
rows = []
for k, g in m.groupby("DE_category"):
    rows.append(dict(DE_category=k, genes=len(g), in_proteome=int(g.in_FungiDB_proteome.sum()),
        eggNOG_description=int(g.eggNOG_Description.notna().sum()),
        v1_tested_by_SignalP_TMHMM_GPI=int(g.SignalP_tested.eq("True").sum()),
        WP_tested_all_4_tools=int((g.WP_InterPro_tested & g.WP_SignalP_tested & g.WP_DeepTMHMM_tested & g.WP_PredGPI_tested).sum()),
        InterPro_IPR_hit=int(g.WP_InterPro_accessions.notna().sum()),
        signal_peptide=int(g.WP_SignalP_prediction.eq("SP").sum()),
        TM=int(g.WP_DeepTMHMM_has_TM.eq(True).sum()),
        GPI=int(g.WP_PredGPI_prediction.eq("GPI-anchored").sum()),
        GPI_without_signal_peptide=int(g.WP_GPI_without_signal_peptide.eq(True).sum()),
        GPI_final=int(g.WP_GPI_final.eq(True).sum())))
sm = pd.DataFrame(rows).sort_values("genes", ascending=False)
sm.to_csv(J(OUT, "summary_by_DE_category.csv"), index=False); print(); print(sm.to_string(index=False))
print("\nDE-table PCGs:", len(m), "| in proteome:", int(m.in_FungiDB_proteome.sum()),
      "| in proteome but missing WP results:", int((m.in_FungiDB_proteome & ~m.WP_SignalP_tested).sum()))
