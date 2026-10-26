## =============================================================================
## Cellular mRNA DESeq2 re-run on STRANDED counts, and comparison with the old run
##
## Input : the 12 x 3 count files from 01_align_and_count.sbatch, copied from
##         HPC 08_Cell_mRNA_corrected_2026-09-25/03_stranded_2026-10-14/counts/
##         into COUNTS_DIR below
## Output: new B8441 / MMC1 DESeq2 results (stranded, and new-unstranded as a control)
##         plus an old-vs-new comparison, all written to OUT_DIR
##
## Three versions are compared so the cause of any change is clear:
##   old            = historical counts (old STAR settings, strand ignored)
##   new_unstranded = new alignment, strand ignored   -> effect of realignment only
##   new_stranded   = new alignment, strand used (-S) -> the corrected result
##
## Base R + dplyr only (stringr/stringi avoided, as in the step 7-9 scripts).
## =============================================================================

library(DESeq2)
library(dplyr)

BASE <- "/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV"
RUN_DIR    <- file.path(BASE, "Differential_expression_analysis/03_Cellular_DE_analysis/mRNA_DE_analysis/STRANDED_recount_2026-10-14")
COUNTS_DIR <- file.path(RUN_DIR, "counts")
OUT_DIR    <- file.path(RUN_DIR, "results")
OLD_COUNTS <- file.path(BASE, "Differential_expression_analysis/02_Input_data/count_matrix_transcripts.rds")
OLD_META   <- file.path(BASE, "Differential_expression_analysis/02_Input_data/sample_metadata.rds")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

SAMPLES <- c("B8441_mRNA_1", "B8441_mRNA_2", "B8441_mRNA_3",
             "B8441_mRNA_CASP_1", "B8441_mRNA_CASP_2", "B8441_mRNA_CASP_3",
             "MMC1_mRNA_1", "MMC1_mRNA_2", "MMC1_mRNA_3",
             "MMC1_mRNA_CASP_1", "MMC1_mRNA_CASP_2", "MMC1_mRNA_CASP_3")

meta <- data.frame(
  row.names = SAMPLES,
  Treatment = factor(ifelse(grepl("CASP", SAMPLES), "Treated", "Untreated"),
                     levels = c("Untreated", "Treated")),
  Strain    = factor(sub("_.*", "", SAMPLES), levels = c("B8441", "MMC1"))
)

## -----------------------------------------------------------------------------
## 1. Build a transcript-level count matrix from the bedtools files
##    bedtools coverage -counts output: chr, start, end, name, score, strand, count
##    Only *_transcript rows are used (same convention as the FINAL Cell sRNA run).
##    Rows are matched by NAME, never by position.
## -----------------------------------------------------------------------------
read_counts <- function(sample, suffix) {
  f <- file.path(COUNTS_DIR, paste0(sample, ".", suffix, ".txt"))
  x <- read.delim(f, header = FALSE, stringsAsFactors = FALSE,
                  col.names = c("chr", "start", "end", "name", "score", "strand", "count"))
  x <- x[grepl("_transcript$", x$name), ]
  # A few names (MSTRG.*) occur twice; sum them by name, as in the FINAL sRNA pipeline
  x %>%
    group_by(name) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    mutate(sample = sample)
}

build_matrix <- function(suffix) {
  long <- bind_rows(lapply(SAMPLES, read_counts, suffix = suffix))
  wide <- reshape(as.data.frame(long), idvar = "name", timevar = "sample", v.names = "count", direction = "wide")
  rownames(wide) <- sub("_transcript$", "", wide$name)
  m <- as.matrix(wide[, paste0("count.", SAMPLES)])
  colnames(m) <- SAMPLES
  m[is.na(m)] <- 0
  m
}

counts_stranded   <- build_matrix("stranded_S")
counts_unstranded <- build_matrix("unstranded")
counts_samestrand <- build_matrix("samestrand_s")

## Strand sanity check per sample (protein-coding genes): should be well above 0.9
pcg_rows <- grepl("^B9J08_", rownames(counts_stranded))
strand_check <- data.frame(
  sample = SAMPLES,
  opposite_strand_reads = colSums(counts_stranded[pcg_rows, ]),
  same_strand_reads     = colSums(counts_samestrand[rownames(counts_stranded)[pcg_rows], ])
)
strand_check$fraction_correct_strand <- with(strand_check,
  opposite_strand_reads / (opposite_strand_reads + same_strand_reads))
print(strand_check)
write.csv(strand_check, file.path(OUT_DIR, "strand_check_per_sample.csv"), row.names = FALSE)
if (any(strand_check$fraction_correct_strand < 0.9)) {
  warning("At least one library has <90% of PCG reads on the expected strand -- check before trusting results")
}

## Old matrix (historical counts, 36 libraries); keep the 12 cellular mRNA columns
old_counts <- as.matrix(readRDS(OLD_COUNTS))
rownames(old_counts) <- sub("_transcript$", "", rownames(old_counts))
old_counts <- old_counts[, SAMPLES]

## -----------------------------------------------------------------------------
## 2. DESeq2: same model as the original cellular mRNA analysis
##    ~ Treatment * Strain, Untreated and B8441 as reference levels.
##    Filter: >=10 reads in >=3 of the 12 libraries. (The old run filtered across all
##    36 libraries and then dropped zero-count rows; that can't be repeated here
##    because only the 12 mRNA libraries were recounted. The filter is applied the
##    same way to all three versions below, so the comparison stays fair.)
## -----------------------------------------------------------------------------
run_deseq <- function(counts, label) {
  counts <- round(counts)  # bedtools counts are already integers; this just enforces it
  keep <- rowSums(counts >= 10) >= 3
  cat(label, ": features kept", sum(keep), "of", nrow(counts), "\n")
  dds <- DESeqDataSetFromMatrix(counts[keep, ], colData = meta, design = ~ Treatment * Strain)
  dds <- DESeq(dds)
  b <- results(dds, contrast = c("Treatment", "Treated", "Untreated"), alpha = 0.05)
  m <- results(dds, contrast = list(c("Treatment_Treated_vs_Untreated",
                                      "TreatmentTreated.StrainMMC1")), alpha = 0.05)
  tidy <- function(r) {
    d <- as.data.frame(r)
    d$Feature <- rownames(d)
    d[, c("Feature", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  }
  out <- list(B8441 = tidy(b), MMC1 = tidy(m), dds = dds)
  write.csv(out$B8441, file.path(OUT_DIR, paste0(label, "_B8441_results.csv")), row.names = FALSE)
  write.csv(out$MMC1,  file.path(OUT_DIR, paste0(label, "_MMC1_results.csv")),  row.names = FALSE)
  saveRDS(dds, file.path(OUT_DIR, paste0(label, "_dds.rds")))
  out
}

res_old        <- run_deseq(old_counts,        "old")
res_unstranded <- run_deseq(counts_unstranded, "new_unstranded")
res_stranded   <- run_deseq(counts_stranded,   "new_stranded")

## -----------------------------------------------------------------------------
## 3. Compare: for each strain and gene type, how much do the results move?
## -----------------------------------------------------------------------------
compare_two <- function(a, b, name_a, name_b, strain) {
  j <- inner_join(a, b, by = "Feature", suffix = c(".a", ".b"))
  j$type <- ifelse(grepl("^LNCRNA", j$Feature), "lncRNA",
                   ifelse(grepl("^B9J08_", j$Feature), "PCG", "other"))
  sig_a <- !is.na(j$padj.a) & j$padj.a < 0.05
  sig_b <- !is.na(j$padj.b) & j$padj.b < 0.05
  j %>%
    mutate(sig_a = sig_a, sig_b = sig_b) %>%
    group_by(type) %>%
    summarise(
      comparison        = paste(name_a, "vs", name_b),
      strain            = strain,
      n_in_both         = n(),
      log2FC_pearson    = cor(log2FoldChange.a, log2FoldChange.b, use = "complete.obs"),
      log2FC_spearman   = cor(log2FoldChange.a, log2FoldChange.b, use = "complete.obs", method = "spearman"),
      sig_in_a          = sum(sig_a),
      sig_in_b          = sum(sig_b),
      sig_in_both       = sum(sig_a & sig_b),
      lost_significance = sum(sig_a & !sig_b),
      gained_significance = sum(!sig_a & sig_b),
      direction_flipped_while_sig_in_both =
        sum(sig_a & sig_b & sign(log2FoldChange.a) != sign(log2FoldChange.b)),
      .groups = "drop"
    )
}

summary_tbl <- bind_rows(
  compare_two(res_old$B8441, res_unstranded$B8441, "old", "new_unstranded", "B8441"),
  compare_two(res_old$MMC1,  res_unstranded$MMC1,  "old", "new_unstranded", "MMC1"),
  compare_two(res_unstranded$B8441, res_stranded$B8441, "new_unstranded", "new_stranded", "B8441"),
  compare_two(res_unstranded$MMC1,  res_stranded$MMC1,  "new_unstranded", "new_stranded", "MMC1"),
  compare_two(res_old$B8441, res_stranded$B8441, "old", "new_stranded", "B8441"),
  compare_two(res_old$MMC1,  res_stranded$MMC1,  "old", "new_stranded", "MMC1")
)
print(as.data.frame(summary_tbl))
write.csv(summary_tbl, file.path(OUT_DIR, "old_vs_new_summary.csv"), row.names = FALSE)

## Per-gene table old vs new_stranded, flagging every gene whose call changed.
## Use this to check whether any gene named in another Results section changed.
per_gene <- function(strain) {
  full_join(res_old[[strain]][, c("Feature", "log2FoldChange", "padj")],
            res_stranded[[strain]][, c("Feature", "log2FoldChange", "padj")],
            by = "Feature", suffix = c("_old", "_stranded")) %>%
    mutate(
      strain = strain,
      sig_old = !is.na(padj_old) & padj_old < 0.05,
      sig_stranded = !is.na(padj_stranded) & padj_stranded < 0.05,
      call_old = ifelse(sig_old, ifelse(log2FoldChange_old > 0, "up", "down"), "ns"),
      call_stranded = ifelse(sig_stranded, ifelse(log2FoldChange_stranded > 0, "up", "down"), "ns"),
      call_changed = call_old != call_stranded
    )
}
genes <- bind_rows(per_gene("B8441"), per_gene("MMC1"))
write.csv(genes, file.path(OUT_DIR, "per_gene_old_vs_stranded.csv"), row.names = FALSE)
cat("\nPCGs whose up/down/ns call changed (old -> stranded):\n")
print(with(genes[grepl("^B9J08_", genes$Feature), ], table(strain, call_changed)))

saveRDS(list(old = res_old[c("B8441", "MMC1")],
             new_unstranded = res_unstranded[c("B8441", "MMC1")],
             new_stranded = res_stranded[c("B8441", "MMC1")]),
        file.path(OUT_DIR, "all_results.rds"))
cat("\nDone. Results in", OUT_DIR, "\n")
