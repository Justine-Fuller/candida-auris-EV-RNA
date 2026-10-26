#!/usr/bin/env Rscript
# run_enrichment_real.R  (v2 Step 6C (3), 2026-10-12)
# Based on v1 21_human_codon_gene_clusterProfiler.R: enrichGO (BP/MF/CC, keyType ENSEMBL) and
# enrichKEGG (bitr -> ENTREZ), BH, universe = all scored genes, default minGSSize 10 / maxGSSize 500.
# Changes vs v1 (technical only):
#  - KEGG data downloaded once (prepare_KEGG) and saved, then passed to enrichKEGG as a GSON object,
#    so the real run and the null use the identical KEGG snapshot.
#  - significant = BH p.adjust < 0.05 (approved spec), taken from the unfiltered table.
#  - gate: the fast null ORA (ora_lib.R) must reproduce every clusterProfiler p-value to 1e-10.
#  - builds the planted term (top 50 real-foreground genes, BP universe) and the tracked-term list.
# r2 (2026-10-12): two optional args, no change to the tests themselves:
#  - arg 5 = saved KEGG GSON (.rds) to reuse instead of downloading, so r2 uses the r1 KEGG snapshot;
#  - arg 6 = r1 real-run dir; its significant terms are added to the tracked list (comparison only).
# Args: <scoring_dir> <outdir> <v1_enrichment_output_dir> <libdir> [kegg_rds] [r1_run_dir]
args <- commandArgs(trailingOnly = TRUE)
sdir <- args[1]; odir <- args[2]; v1dir <- args[3]; libdir <- args[4]
kegg_rds <- if (length(args) >= 5) args[5] else NA
r1dir <- if (length(args) >= 6) args[6] else NA
suppressMessages({ library(clusterProfiler); library(org.Hs.eg.db) })
source(file.path(libdir, "ora_lib.R"))
dir.create(odir, showWarnings = FALSE, recursive = TRUE)
set.seed(20261012)

fg <- readLines(file.path(sdir, "foreground_top10pct.txt"))
bg <- readLines(file.path(sdir, "background_all_scored.txt"))
cat("Foreground:", length(fg), " Background:", length(bg), "\n")

sink(file.path(odir, "session_and_versions.txt"))
print(sessionInfo()); print(clusterProfiler:::map_to_entrezid)
sink()

gate <- data.frame()
gs_all <- list(); tables <- list()

# ---- GO ----
univ_go <- fg_for_go(bg)
gene_go <- fg_for_go(fg)
for (ont in c("BP", "MF", "CC")) {
  ego <- enrichGO(gene = fg, universe = bg, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = ont,
                  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
  df <- as.data.frame(ego)
  df$ontology <- ont
  tables[[ont]] <- df
  gson <- clusterProfiler:::get_GO_data(org.Hs.eg.db, ont, "ENTREZID")
  saveRDS(gson, file.path(odir, sprintf("gson_GO_%s.rds", ont)))
  gs <- build_gene_sets(gson, univ_go)
  gs_all[[ont]] <- gs
  mine <- ora_terms(gene_go, gs)
  m <- merge(df[, c("ID", "pvalue", "p.adjust", "Count")], mine, by = "ID", all = TRUE, suffixes = c(".cp", ".fast"))
  gate <- rbind(gate, data.frame(ontology = ont, n_terms_cp = nrow(df), n_terms_fast = nrow(mine),
    same_term_ids = setequal(df$ID, mine$ID), same_counts = isTRUE(all(m$Count.cp == m$Count.fast)),
    max_abs_diff_p = max(abs(m$pvalue.cp - m$pvalue.fast)),
    max_abs_diff_padj = max(abs(m$p.adjust.cp - m$p.adjust.fast)),
    n_universe = length(gs$extID), n_terms_tested_size_10_500 = length(gs$geneSets)))
}

# ---- KEGG ----
KEGG_DATA <- if (!is.na(kegg_rds)) readRDS(kegg_rds) else clusterProfiler:::prepare_KEGG("hsa", "KEGG", "kegg")
if (!is.na(kegg_rds)) cat("KEGG snapshot reused from", kegg_rds, "\n")
saveRDS(KEGG_DATA, file.path(odir, "gson_KEGG_hsa.rds"))
gene_kegg <- fg_for_kegg(fg); univ_kegg <- fg_for_kegg(bg)
cat("Foreground mapped to ENTREZ:", length(gene_kegg), "of", length(fg), "\n")
cat("Background mapped to ENTREZ:", length(univ_kegg), "of", length(bg), "\n")
ek <- enrichKEGG(gene = gene_kegg, organism = KEGG_DATA, universe = univ_kegg,
                 pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1)
dfk <- as.data.frame(setReadable(ek, OrgDb = org.Hs.eg.db, keyType = "ENTREZID"))
dfk$ontology <- "KEGG"
tables[["KEGG"]] <- dfk
gs <- build_gene_sets(KEGG_DATA, univ_kegg); gs_all[["KEGG"]] <- gs
mine <- ora_terms(gene_kegg, gs)
m <- merge(dfk[, c("ID", "pvalue", "p.adjust", "Count")], mine, by = "ID", all = TRUE, suffixes = c(".cp", ".fast"))
gate <- rbind(gate, data.frame(ontology = "KEGG", n_terms_cp = nrow(dfk), n_terms_fast = nrow(mine),
  same_term_ids = setequal(dfk$ID, mine$ID), same_counts = isTRUE(all(m$Count.cp == m$Count.fast)),
  max_abs_diff_p = max(abs(m$pvalue.cp - m$pvalue.fast)),
  max_abs_diff_padj = max(abs(m$p.adjust.cp - m$p.adjust.fast)),
  n_universe = length(gs$extID), n_terms_tested_size_10_500 = length(gs$geneSets)))
saveRDS(gs_all, file.path(odir, "gene_sets_frozen.rds"))
saveRDS(list(go_universe = univ_go, kegg_universe = univ_kegg), file.path(odir, "universes.rds"))

gate$PASS <- gate$same_term_ids & gate$same_counts & gate$max_abs_diff_p <= 1e-10 & gate$max_abs_diff_padj <= 1e-10
gate$PASS[is.na(gate$PASS)] <- FALSE
write.table(gate, file.path(odir, "gate_fast_ora_vs_clusterProfiler.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(gate)

# ---- write tables ----
v1sig <- list(); r1sig <- list()
for (ont in names(tables)) {
  df <- tables[[ont]]
  tag <- if (ont == "KEGG") "KEGG" else paste0("GO_", ont)
  write.table(df, file.path(odir, paste0(tag, "_all.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  sig <- df[df$p.adjust < 0.05, ]
  write.table(sig, file.path(odir, paste0(tag, "_significant.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  v1f <- file.path(v1dir, paste0(tag, "_significant.tsv"))
  v1 <- read.delim(v1f, quote = "", comment.char = "")
  v1sig[[ont]] <- v1$ID
  if (!is.na(r1dir)) {
    r1 <- read.delim(file.path(r1dir, paste0(tag, "_significant.tsv")), quote = "", comment.char = "")
    r1sig[[ont]] <- r1$ID
  }
  cat(sprintf("%s: %d terms tested with Count>0; %d with BH p.adjust < 0.05 (v1-style filter p.adjust<0.05 & qvalue<0.2: %d)\n",
              tag, nrow(df), nrow(sig), sum(df$p.adjust < 0.05 & df$qvalue < 0.2, na.rm = TRUE)))
}

# ---- tracked terms: real significant + v1 significant + HLA/MHC/skin theme terms ----
tracked <- do.call(rbind, lapply(names(tables), function(ont) {
  df <- tables[[ont]]
  ids <- unique(c(df$ID[df$p.adjust < 0.05], v1sig[[ont]],
                  df$ID[grepl(THEME_REGEX, df$Description, ignore.case = TRUE)], r1sig[[ont]]))
  allterms <- names(gs_all[[ont]]$geneSets)
  data.frame(ontology = ont, ID = ids, in_tested_sets = ids %in% allterms,
             real_sig = ids %in% df$ID[df$p.adjust < 0.05], v1_sig = ids %in% v1sig[[ont]],
             theme = ids %in% df$ID[grepl(THEME_REGEX, df$Description, ignore.case = TRUE)],
             r1_sig = ids %in% r1sig[[ont]])
}))
write.table(tracked, file.path(odir, "tracked_terms.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("Tracked terms:", nrow(tracked), " (not in tested sets:", sum(!tracked$in_tested_sets), ")\n")

# ---- planted term: top 50 real-foreground genes (rank order) that are in the BP universe ----
m <- suppressWarnings(suppressMessages(bitr(fg, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)))
m <- m[order(match(m$ENSEMBL, fg)), ]
pl <- unique(m$ENTREZID[m$ENTREZID %in% gs_all[["BP"]]$extID])[1:50]
planted <- list(geneSets = list(PLANTED_top50 = pl), extID = gs_all[["BP"]]$extID)
saveRDS(planted, file.path(odir, "planted_term.rds"))
pr <- ora_terms(gene_go, planted)
writeLines(c(sprintf("planted term: %d genes (top-ranked real foreground genes, ENTREZ, in BP universe)", length(pl)),
             sprintf("real count in foreground: %d; real p = %.6g", pr$Count, pr$pvalue)),
           file.path(odir, "planted_term_real.txt"))
cat("Planted real p:", pr$pvalue, "\n")
cat("ALL GATES (fast ORA):", if (all(gate$PASS)) "PASS" else "FAIL", "\n")
