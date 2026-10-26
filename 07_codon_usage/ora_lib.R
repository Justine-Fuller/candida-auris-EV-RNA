# ora_lib.R  (v2 Step 6C (3), 2026-10-12)
# Fast over-representation for the codon-set null. Uses the SAME function clusterProfiler 4.20
# calls internally (enrichGO/enrichKEGG -> enrichit::ora_gson -> enrichit::ora) on gene sets built
# exactly as ora_gson builds them, so per-term p-values are identical by construction; the gate in
# run_enrichment_real.R checks this against enrichGO/enrichKEGG output for the real codon set.

build_gene_sets <- function(gson, universe, minGSSize = 10, maxGSSize = 500) {
  gsid2gene <- gson@gsid2gene
  extID <- intersect(unique(gsid2gene$gene), universe)          # as ora_gson (force_universe FALSE)
  geneSets <- split(gsid2gene$gene, gsid2gene$gsid)
  geneSets <- lapply(geneSets, intersect, extID)
  sz <- vapply(geneSets, length, integer(1))
  geneSets <- geneSets[minGSSize <= sz & sz <= maxGSSize]
  list(geneSets = geneSets, extID = extID)
}

ora_terms <- function(gene, gs) {
  gene <- as.character(unique(gene))
  res <- enrichit::ora(gene, gs$geneSets, universe = gs$extID)
  res <- res[res$Count > 0, , drop = FALSE]                     # as ora_gson
  res$p.adjust <- p.adjust(res$pvalue, method = "BH")
  res[, c("ID", "Count", "pvalue", "p.adjust")]
}

# foreground (Ensembl IDs) -> gene IDs used by each analysis, exactly as v1 / enrichGO do it
fg_for_go <- function(ens) clusterProfiler:::map_to_entrezid(ens, fromType = "ENSEMBL", OrgDb = org.Hs.eg.db)
fg_for_kegg <- function(ens) {
  m <- suppressWarnings(suppressMessages(
    clusterProfiler::bitr(ens, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)))
  m$ENTREZID
}

THEME_REGEX <- "MHC|HLA|antigen processing and presentation|allograft|graft-versus-host|autoimmune thyroid|type I diabetes|epiderm|keratin|skin|cornified"
