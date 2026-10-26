## ============================================================================
## clusterProfiler analysis for Candida auris DE genes
## ORA (GO + KEGG) and GSEA, using custom eggNOG-derived term mappings
## because C. auris has no OrgDb / Bioconductor annotation package.
## ============================================================================
##
## WHAT THIS SCRIPT DOES
## 1. Builds custom GO and KEGG "TERM2GENE" maps from your full-proteome
##    eggNOG annotation (this is what replaces org.Xx.eg.db for a
##    non-model organism).
## 2. Defines a proper statistical universe/background = every
##    protein-coding gene that was actually TESTED in your DESeq2 run
##    AND has at least one functional annotation.
## 3. Runs ORA (over-representation analysis) for GO and KEGG, separately
##    for your up- and down-regulated gene sets.
## 4. Runs GSEA (gene set enrichment analysis) for GO and KEGG using the
##    full ranked gene list (DESeq2 Wald statistic), not just the
##    significant genes.
## 5. Saves all results as CSV + a dotplot PNG for each analysis.
##
## EDIT THE FILE PATHS IN THE "USER SETTINGS" BLOCK BELOW, THEN SOURCE/RUN.
## ============================================================================


## ---- 0. Install packages (run once, comment out afterwards) --------------
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler", "enrichplot", "GO.db"))
# install.packages(c("dplyr", "readr", "tidyr", "stringr", "ggplot2"))


## ---- 1. Libraries ----------------------------------------------------------
library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)

## clusterProfiler loads AnnotationDbi in the background, and
## AnnotationDbi::select() masks dplyr::select() (same for filter/rename/
## mutate in some setups). Force these back to the dplyr versions so the
## tidyverse pipes below work as written.
select   <- dplyr::select
filter   <- dplyr::filter
rename   <- dplyr::rename
mutate   <- dplyr::mutate
arrange  <- dplyr::arrange
distinct <- dplyr::distinct

## Optional: GO.db lets us turn "GO:0006355" into a readable name like
## "regulation of DNA-templated transcription". If it's not installed,
## the script falls back to using the GO ID itself as the label.
has_godb <- requireNamespace("GO.db", quietly = TRUE)


## ---- 2. USER SETTINGS: edit these paths -----------------------------------
full_proteome_emapper <- "/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/PCG annotation/full_proteome_annotation.emapper.annotations"
up_gene_file          <- "/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline/Gene_IDs_Upregulated_renamed_Annotated.txt"
down_gene_file        <- "/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline/Gene_IDs_Downregulated_renamed_Annotated.txt"
deseq2_full_results    <- "/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline/03_mRNA_Cells_DESeq2_Results.csv"
output_dir            <- "/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/PCG annotation/clusterProfiler_results"

dir.create(output_dir, showWarnings = FALSE)


## ---- 3. Helper: strip transcript/protein suffix off eggNOG "query" IDs ----
## e.g. "B9J08_000002-t37_1-p1" -> "B9J08_000002"
base_id <- function(x) sub("-.*", "", x)


## ---- 4. Load full-proteome eggNOG annotation (background) ----------------
bg_raw <- read_tsv(
  full_proteome_emapper,
  comment = "##",
  na = "-",
  show_col_types = FALSE
)
# The header line itself starts with "#query" - strip the leading "#"
colnames(bg_raw)[1] <- "query"

bg <- bg_raw %>%
  mutate(Gene_ID = base_id(query)) %>%
  distinct(Gene_ID, .keep_all = TRUE)

cat("Full proteome annotation loaded:", nrow(bg), "genes\n")


## ---- 5. Build GO TERM2GENE / TERM2NAME ------------------------------------
term2gene_go <- bg %>%
  filter(!is.na(GOs)) %>%
  select(Gene_ID, GOs) %>%
  separate_rows(GOs, sep = ",") %>%
  filter(GOs != "") %>%
  transmute(term = GOs, gene = Gene_ID) %>%
  distinct()

if (has_godb) {
  library(GO.db)
  go_terms <- unique(term2gene_go$term)
  go_names <- suppressMessages(
    AnnotationDbi::select(GO.db, keys = go_terms, columns = "TERM", keytype = "GOID")
  )
  term2name_go <- go_names %>%
    transmute(term = GOID, name = ifelse(is.na(TERM), GOID, TERM)) %>%
    distinct(term, .keep_all = TRUE)
} else {
  message("GO.db not installed - using raw GO IDs as labels. ",
          "Install GO.db for readable term names.")
  term2name_go <- term2gene_go %>% distinct(term) %>% mutate(name = term)
}

cat("GO TERM2GENE built:", nrow(term2gene_go), "gene-term pairs,",
    length(unique(term2gene_go$term)), "unique GO terms\n")


## ---- 6. Build KEGG TERM2GENE / TERM2NAME ----------------------------------
## We use the "map"-prefixed pathway IDs (generic KEGG pathway maps,
## e.g. map00010) rather than the "ko"-prefixed ones, since map IDs are
## the standard identifiers used in KEGG pathway browsers/figures.
term2gene_kegg <- bg %>%
  filter(!is.na(KEGG_Pathway)) %>%
  select(Gene_ID, KEGG_Pathway) %>%
  separate_rows(KEGG_Pathway, sep = ",") %>%
  filter(str_starts(KEGG_Pathway, "map")) %>%
  transmute(term = KEGG_Pathway, gene = Gene_ID) %>%
  distinct()

## Try to fetch human-readable KEGG pathway names via KEGGREST (needs
## internet access from wherever you run this). Falls back to raw IDs
## if unavailable (e.g. offline HPC compute node).
term2name_kegg <- tryCatch({
  if (!requireNamespace("KEGGREST", quietly = TRUE)) stop("KEGGREST not installed")
  library(KEGGREST)
  path_list <- keggList("pathway")             # named vector: "path:map00010" -> "Glycolysis..."
  tibble(term = sub("^path:", "", names(path_list)), name = unname(path_list)) %>%
    filter(term %in% unique(term2gene_kegg$term))
}, error = function(e) {
  message("Could not fetch KEGG pathway names (offline or KEGGREST missing). ",
          "Using raw pathway IDs as labels instead.")
  term2gene_kegg %>% distinct(term) %>% mutate(name = term)
})

cat("KEGG TERM2GENE built:", nrow(term2gene_kegg), "gene-pathway pairs,",
    length(unique(term2gene_kegg$term)), "unique pathways\n")


## ---- 7. Load your DE gene lists and full DESeq2 results -------------------
up_genes   <- readLines(up_gene_file)
down_genes <- readLines(down_gene_file)

deseq_full <- read_csv(deseq2_full_results, show_col_types = FALSE) %>%
  filter(str_starts(Feature, "B9J08_"))   # keep protein-coding genes only

cat("Up genes:", length(up_genes), " | Down genes:", length(down_genes),
    " | DESeq2 protein-coding features tested:", nrow(deseq_full), "\n")


## ---- 8. Define the statistical universe -----------------------------------
## Universe = genes that were BOTH tested in DESeq2 AND have at least one
## functional annotation (GO or KEGG). This is the correct background for
## ORA - using only your DE gene lists as their own background would
## inflate significance.
universe_go   <- intersect(deseq_full$Feature, unique(term2gene_go$gene))
universe_kegg <- intersect(deseq_full$Feature, unique(term2gene_kegg$gene))

cat("GO universe size:", length(universe_go), "\n")
cat("KEGG universe size:", length(universe_kegg), "\n")


## ============================================================================
## PART A: OVER-REPRESENTATION ANALYSIS (ORA)
## ============================================================================

run_ora <- function(genes, term2gene, term2name, universe, label) {
  enricher(
    gene          = genes,
    universe      = universe,
    TERM2GENE     = term2gene,
    TERM2NAME     = term2name,
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    minGSSize     = 3,
    maxGSSize     = 500
  )
}

## --- GO ORA ---
ora_go_up   <- run_ora(up_genes,   term2gene_go, term2name_go, universe_go, "GO_Up")
ora_go_down <- run_ora(down_genes, term2gene_go, term2name_go, universe_go, "GO_Down")

## --- KEGG ORA ---
ora_kegg_up   <- run_ora(up_genes,   term2gene_kegg, term2name_kegg, universe_kegg, "KEGG_Up")
ora_kegg_down <- run_ora(down_genes, term2gene_kegg, term2name_kegg, universe_kegg, "KEGG_Down")

## --- Save ORA results ---
save_ora <- function(result, name) {
  if (is.null(result) || nrow(as.data.frame(result)) == 0) {
    message("No significant terms found for: ", name)
    return(invisible(NULL))
  }
  write_csv(as.data.frame(result), file.path(output_dir, paste0(name, ".csv")))
  p <- dotplot(result, showCategory = 20) + ggtitle(name)
  ggsave(file.path(output_dir, paste0(name, "_dotplot.png")), p,
         width = 9, height = 7, dpi = 300)
}

save_barplot <- function(result, name) {
  if (is.null(result) || nrow(as.data.frame(result)) == 0) return(invisible(NULL))
  p <- barplot(result, showCategory = 15) + ggtitle(name)
  ggsave(file.path(output_dir, paste0(name, "_barplot.png")), p,
         width = 9, height = 7, dpi = 300)
}

save_ora(ora_go_up,     "ORA_GO_Upregulated")
save_barplot(ora_go_up, "ORA_GO_Upregulated")

save_ora(ora_go_down,   "ORA_GO_Downregulated")
save_barplot(ora_go_down, "ORA_GO_Downregulated")

save_ora(ora_kegg_up,   "ORA_KEGG_Upregulated")
save_barplot(ora_kegg_up, "ORA_KEGG_Upregulated")

save_ora(ora_kegg_down, "ORA_KEGG_Downregulated")
save_barplot(ora_kegg_down, "ORA_KEGG_Downregulated")

## --- Combined Up vs Down comparison plot (GO) ---
cc_go <- compareCluster(
  geneClusters = list(Upregulated = up_genes, Downregulated = down_genes),
  fun          = "enricher",
  universe     = universe_go,
  TERM2GENE    = term2gene_go,
  TERM2NAME    = term2name_go,
  pvalueCutoff = 0.05,
  minGSSize    = 3
)
if (!is.null(cc_go) && nrow(as.data.frame(cc_go)) > 0) {
  write_csv(as.data.frame(cc_go), file.path(output_dir, "ORA_GO_UpVsDown_comparison.csv"))
  p <- dotplot(cc_go, showCategory = 15) + ggtitle("GO ORA: Up vs Down")
  ggsave(file.path(output_dir, "ORA_GO_UpVsDown_dotplot.png"), p,
         width = 10, height = 8, dpi = 300)
}


## ============================================================================
## PART B: GENE SET ENRICHMENT ANALYSIS (GSEA)
## ============================================================================
## GSEA needs a ranked list of ALL tested genes (not just significant ones),
## ranked by a continuous metric. We use the DESeq2 Wald statistic ("stat")
## rather than raw log2FoldChange, since "stat" already accounts for each
## gene's variance/confidence and is the standard recommendation for
## ranking genes going into GSEA.

gene_rank <- deseq_full %>%
  filter(!is.na(stat)) %>%
  arrange(desc(stat))

gene_list <- gene_rank$stat
names(gene_list) <- gene_rank$Feature

cat("GSEA ranked gene list length:", length(gene_list), "\n")

run_gsea <- function(term2gene, term2name, label) {
  GSEA(
    geneList     = gene_list,
    TERM2GENE    = term2gene,
    TERM2NAME    = term2name,
    pvalueCutoff = 0.05,
    minGSSize    = 3,
    maxGSSize    = 500,
    seed         = TRUE
  )
}

gsea_go   <- run_gsea(term2gene_go,   term2name_go,   "GSEA_GO")
gsea_kegg <- run_gsea(term2gene_kegg, term2name_kegg, "GSEA_KEGG")

save_gsea <- function(result, name) {
  if (is.null(result) || nrow(as.data.frame(result)) == 0) {
    message("No significant terms found for: ", name)
    return(invisible(NULL))
  }
  write_csv(as.data.frame(result), file.path(output_dir, paste0(name, ".csv")))
  p <- dotplot(result, showCategory = 20) + ggtitle(name)
  ggsave(file.path(output_dir, paste0(name, "_dotplot.png")), p,
         width = 9, height = 7, dpi = 300)
}

save_gsea(gsea_go,   "GSEA_GO")
save_gsea(gsea_kegg, "GSEA_KEGG")


## ============================================================================
## DONE
## ============================================================================
cat("\n=====================================================\n")
cat("clusterProfiler analysis complete.\n")
cat("Results saved to:", normalizePath(output_dir), "\n")
cat("=====================================================\n")
