## =============================================================================
## lncRNA-PCG neighbour analysis, redone on STRANDED counts
## ("lncRNAs transcribed from neighbouring positions of PCGs", steps 4-9)
##
## Changes from the original run (see reports/lncRNA_PCG_neighbour_check_2026-10-14.md):
##   1. Expression comes from the stranded DESeq2 run (02_DESeq2_stranded_and_compare.R),
##      so antisense overlapping pairs no longer share reads.
##   2. Partner choice is no longer arbitrary: an lncRNA that overlaps a PCG is paired
##      with the PCG it overlaps by the MOST bases; a non-overlapping lncRNA keeps its
##      nearest PCG. Ties are counted and reported.
##   3. "Overlapping" now means at least 1 shared base. Ranges that only touch
##      end-to-end go into "<=1 kb".
##   4. A chance baseline is added: PCG partners are shuffled 1,000 times to show how
##      often random pairs look "concordant".
## Classification and ranking rules are unchanged from steps 8-9.
## Base R + dplyr + Bioconductor only.
## =============================================================================

library(rtracklayer)
library(GenomicRanges)
library(dplyr)

BASE    <- "/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV"
RUN_DIR <- file.path(BASE, "Differential_expression_analysis/03_Cellular_DE_analysis/mRNA_DE_analysis/STRANDED_recount_2026-10-14")
GTF     <- file.path(BASE, "Differential_expression_analysis/01_Reference_files/final.gtf")
OLD_CLASSIFIED <- file.path(BASE, "PCG_annotation/Cellular_RNA_analysis/04_lncRNA_PCG_analysis/lncRNA_PCG_classified_by_strain.csv")
OUT_DIR <- file.path(RUN_DIR, "lncRNA_PCG")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
SIG <- 0.05
set.seed(20261014)

## -----------------------------------------------------------------------------
## Steps 1-3: lncRNA and PCG ranges (unchanged)
## -----------------------------------------------------------------------------
gtf <- import(GTF)
tx  <- gtf[gtf$type == "transcript"]
tx  <- tx[!duplicated(tx$gene_id)]
lnc <- tx[grepl("^LNCRNA", tx$gene_id)]
pcg <- tx[grepl("^B9J08_", tx$gene_id)]
cat("lncRNAs:", length(lnc), " PCGs:", length(pcg), "\n")

## -----------------------------------------------------------------------------
## Step 4 (revised): choose each lncRNA's partner PCG
## -----------------------------------------------------------------------------
ov <- findOverlaps(lnc, pcg, ignore.strand = TRUE)
ov_df <- data.frame(
  q = queryHits(ov), s = subjectHits(ov),
  overlap_bp = width(pintersect(lnc[queryHits(ov)], pcg[subjectHits(ov)], ignore.strand = TRUE))
)
n_overlapping_pcgs <- tabulate(ov_df$q, nbins = length(lnc))

# Largest overlap wins; a tie is broken by the PCG's ID so it is reproducible
ov_best <- ov_df %>%
  mutate(pcg_id = pcg$gene_id[s]) %>%
  group_by(q) %>%
  mutate(n_tied_best = sum(overlap_bp == max(overlap_bp))) %>%
  arrange(desc(overlap_bp), pcg_id, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

pairs <- data.frame(
  lncRNA = lnc$gene_id,
  n_overlapping_PCGs = n_overlapping_pcgs,
  PCG = NA_character_, distance_bp = NA_integer_, overlap_bp = 0L,
  n_tied_best = NA_integer_, stringsAsFactors = FALSE
)
pairs$PCG[ov_best$q]         <- ov_best$pcg_id
pairs$distance_bp[ov_best$q] <- 0L
pairs$overlap_bp[ov_best$q]  <- ov_best$overlap_bp
pairs$n_tied_best[ov_best$q] <- ov_best$n_tied_best

# lncRNAs that overlap nothing: nearest PCG
no_ov <- which(n_overlapping_pcgs == 0)
nn <- distanceToNearest(lnc[no_ov], pcg, ignore.strand = TRUE)
pairs$PCG[no_ov[queryHits(nn)]]         <- pcg$gene_id[subjectHits(nn)]
pairs$distance_bp[no_ov[queryHits(nn)]] <- mcols(nn)$distance

# Orientation and how much of the lncRNA lies inside its partner
pidx <- match(pairs$PCG, pcg$gene_id)
pairs$orientation <- ifelse(as.character(strand(lnc)) == as.character(strand(pcg[pidx])),
                            "same strand", "antisense")
pairs$lncRNA_length <- width(lnc)
pairs$fraction_lncRNA_in_PCG <- pairs$overlap_bp / pairs$lncRNA_length
pairs$distance_category <- ifelse(
  pairs$overlap_bp > 0, "Overlapping",
  as.character(cut(pairs$distance_bp, breaks = c(-1, 1000, 5000, 10000, Inf),
                   labels = c("<=1 kb", "1-5 kb", "5-10 kb", ">10 kb")))
)

cat("\nlncRNAs overlapping >1 PCG:", sum(n_overlapping_pcgs > 1),
    "\nOf those, ties for largest overlap:", sum(pairs$n_tied_best > 1, na.rm = TRUE),
    "\nUnique PCGs used:", length(unique(pairs$PCG)), "\n")
print(table(pairs$distance_category))

## -----------------------------------------------------------------------------
## Steps 6-7: join with the STRANDED DESeq2 results
## -----------------------------------------------------------------------------
res <- readRDS(file.path(RUN_DIR, "results", "all_results.rds"))$new_stranded

add_strain <- function(df, r, strain) {
  r <- r[, c("Feature", "baseMean", "log2FoldChange", "padj")]
  l <- setNames(r, c("lncRNA", paste0(strain, "_lncRNA_", c("baseMean", "log2FC", "padj"))))
  p <- setNames(r, c("PCG",    paste0(strain, "_PCG_",    c("baseMean", "log2FC", "padj"))))
  df %>% left_join(l, by = "lncRNA") %>% left_join(p, by = "PCG")
}
expr <- pairs %>% add_strain(res$B8441, "B8441") %>% add_strain(res$MMC1, "MMC1")
stopifnot(nrow(expr) == length(lnc), !any(duplicated(expr$lncRNA)))

## -----------------------------------------------------------------------------
## Step 8: classify (same rules as before, plus "PCG not tested", which can now
## happen because the filter differs from the old run)
## -----------------------------------------------------------------------------
classify <- function(lb, ll, lp, pb, pl, pp) {
  ls <- !is.na(lp) & lp < SIG
  ps <- !is.na(pp) & pp < SIG
  out <- ifelse(is.na(lb), "lncRNA not tested",
         ifelse(is.na(pb), "PCG not tested",
         ifelse(ls & ps & ll > 0 & pl > 0, "Co-upregulated",
         ifelse(ls & ps & ll < 0 & pl < 0, "Co-downregulated",
         ifelse(ls & ps, "Opposite regulation",
         ifelse(ls | ps, "One significant", "Neither significant"))))))
  out
}
for (s in c("B8441", "MMC1")) {
  expr[[paste0(s, "_relationship")]] <- classify(
    expr[[paste0(s, "_lncRNA_baseMean")]], expr[[paste0(s, "_lncRNA_log2FC")]], expr[[paste0(s, "_lncRNA_padj")]],
    expr[[paste0(s, "_PCG_baseMean")]],    expr[[paste0(s, "_PCG_log2FC")]],    expr[[paste0(s, "_PCG_padj")]])
}
conc <- c("Co-upregulated", "Co-downregulated")
expr$cross_strain_relationship <- with(expr, case_when(
  B8441_relationship %in% c("lncRNA not tested", "PCG not tested") |
    MMC1_relationship %in% c("lncRNA not tested", "PCG not tested") ~ "Not tested",
  B8441_relationship == "Co-upregulated"   & MMC1_relationship == "Co-upregulated"   ~ "Same co-upregulated in both strains",
  B8441_relationship == "Co-downregulated" & MMC1_relationship == "Co-downregulated" ~ "Same co-downregulated in both strains",
  B8441_relationship %in% conc & MMC1_relationship %in% conc ~ "Co-regulated in both but different direction",
  B8441_relationship %in% conc ~ "Co-regulated in B8441 only",
  MMC1_relationship  %in% conc ~ "Co-regulated in MMC1 only",
  TRUE ~ "Not co-regulated in either"
))
write.csv(expr, file.path(OUT_DIR, "lncRNA_PCG_classified_by_strain_STRANDED.csv"), row.names = FALSE)

rel_levels <- c("Co-upregulated", "Co-downregulated", "Opposite regulation", "One significant",
                "Neither significant", "lncRNA not tested", "PCG not tested")
rel_tbl <- data.frame(Relationship = rel_levels,
  B8441 = as.integer(table(factor(expr$B8441_relationship, rel_levels))),
  MMC1  = as.integer(table(factor(expr$MMC1_relationship,  rel_levels))))
rel_tbl$B8441_pct <- round(100 * rel_tbl$B8441 / nrow(expr), 2)
rel_tbl$MMC1_pct  <- round(100 * rel_tbl$MMC1  / nrow(expr), 2)
print(rel_tbl)
write.csv(rel_tbl, file.path(OUT_DIR, "relationship_summary_STRANDED.csv"), row.names = FALSE)
cross_tbl <- as.data.frame(table(expr$cross_strain_relationship))
print(cross_tbl)
write.csv(cross_tbl, file.path(OUT_DIR, "cross_strain_summary_STRANDED.csv"), row.names = FALSE)

## -----------------------------------------------------------------------------
## Chance baseline: shuffle PCG partners among tested pairs, 1,000 times.
## "If lncRNAs were paired with random PCGs, how often would a pair look concordant?"
## -----------------------------------------------------------------------------
baseline <- function(df, strain, n_perm = 1000) {
  t <- df[!(df[[paste0(strain, "_relationship")]] %in% c("lncRNA not tested", "PCG not tested")), ]
  l_sig <- !is.na(t[[paste0(strain, "_lncRNA_padj")]]) & t[[paste0(strain, "_lncRNA_padj")]] < SIG
  p_sig <- !is.na(t[[paste0(strain, "_PCG_padj")]])    & t[[paste0(strain, "_PCG_padj")]]    < SIG
  l_dir <- sign(t[[paste0(strain, "_lncRNA_log2FC")]])
  p_dir <- sign(t[[paste0(strain, "_PCG_log2FC")]])
  observed <- mean(t[[paste0(strain, "_relationship")]] %in% conc)
  null <- replicate(n_perm, { k <- sample(nrow(t)); mean(l_sig & p_sig[k] & l_dir == p_dir[k]) })
  data.frame(strain = strain, n_tested_pairs = nrow(t),
             observed_concordant_fraction = observed,
             chance_mean = mean(null),
             chance_95pct_low = quantile(null, 0.025), chance_95pct_high = quantile(null, 0.975),
             fold_over_chance = observed / mean(null),
             empirical_p = (sum(null >= observed) + 1) / (n_perm + 1))
}
# Run for all pairs, and separately within overlapping and non-overlapping pairs
base_tbl <- bind_rows(lapply(c("B8441", "MMC1"), function(s) bind_rows(
  cbind(subset = "all pairs",       baseline(expr, s)),
  cbind(subset = "overlapping",     baseline(expr[expr$distance_category == "Overlapping", ], s)),
  cbind(subset = "not overlapping", baseline(expr[expr$distance_category != "Overlapping", ], s))
)))
rownames(base_tbl) <- NULL
print(base_tbl)
write.csv(base_tbl, file.path(OUT_DIR, "chance_baseline_STRANDED.csv"), row.names = FALSE)

## Diagnostic: do overlapping pairs still move together much more than nearby ones?
diag_tbl <- bind_rows(lapply(c("B8441", "MMC1"), function(s) {
  expr %>% filter(!is.na(.data[[paste0(s, "_lncRNA_log2FC")]]), !is.na(.data[[paste0(s, "_PCG_log2FC")]])) %>%
    group_by(distance_category) %>%
    summarise(strain = s, n = n(),
              log2FC_correlation = cor(.data[[paste0(s, "_lncRNA_log2FC")]], .data[[paste0(s, "_PCG_log2FC")]]),
              concordant_fraction = mean(.data[[paste0(s, "_relationship")]] %in% conc),
              .groups = "drop")
}))
print(as.data.frame(diag_tbl))
write.csv(diag_tbl, file.path(OUT_DIR, "overlap_diagnostic_STRANDED.csv"), row.names = FALSE)

## -----------------------------------------------------------------------------
## Step 9: rank concordant pairs (same hierarchical scheme as before)
## -----------------------------------------------------------------------------
cat_order <- c("Overlapping", "<=1 kb", "1-5 kb", "5-10 kb", ">10 kb")
rank_strain <- function(s) {
  expr %>% filter(.data[[paste0(s, "_relationship")]] %in% conc) %>%
    mutate(worst_padj = pmax(.data[[paste0(s, "_lncRNA_padj")]], .data[[paste0(s, "_PCG_padj")]]),
           min_abs_log2FC = pmin(abs(.data[[paste0(s, "_lncRNA_log2FC")]]), abs(.data[[paste0(s, "_PCG_log2FC")]])),
           cat_rank = match(distance_category, cat_order)) %>%
    arrange(cat_rank, distance_bp, worst_padj, desc(min_abs_log2FC)) %>%
    mutate(rank = row_number()) %>% select(-cat_rank)
}
ranked_B8441 <- rank_strain("B8441")
ranked_MMC1  <- rank_strain("MMC1")
ranked_shared <- expr %>%
  filter(cross_strain_relationship %in% c("Same co-upregulated in both strains",
                                          "Same co-downregulated in both strains")) %>%
  mutate(worst_padj_across_both_strains = pmax(B8441_lncRNA_padj, B8441_PCG_padj, MMC1_lncRNA_padj, MMC1_PCG_padj),
         min_abs_log2FC_across_both_strains = pmin(abs(B8441_lncRNA_log2FC), abs(B8441_PCG_log2FC),
                                                   abs(MMC1_lncRNA_log2FC), abs(MMC1_PCG_log2FC)),
         cat_rank = match(distance_category, cat_order)) %>%
  arrange(cat_rank, distance_bp, worst_padj_across_both_strains, desc(min_abs_log2FC_across_both_strains)) %>%
  mutate(rank = row_number()) %>% select(-cat_rank)

write.csv(ranked_B8441,  file.path(OUT_DIR, "B8441_ranked_coregulated_pairs_STRANDED.csv"), row.names = FALSE)
write.csv(ranked_MMC1,   file.path(OUT_DIR, "MMC1_ranked_coregulated_pairs_STRANDED.csv"),  row.names = FALSE)
write.csv(ranked_shared, file.path(OUT_DIR, "Shared_strain_ranked_coregulated_pairs_STRANDED.csv"), row.names = FALSE)
cat("\nTop 10 shared candidates (stranded):\n")
print(head(ranked_shared[, c("rank", "lncRNA", "PCG", "overlap_bp", "cross_strain_relationship")], 10))

dist_tbl <- bind_rows(
  data.frame(set = "B8441",  table(factor(ranked_B8441$distance_category,  cat_order))),
  data.frame(set = "MMC1",   table(factor(ranked_MMC1$distance_category,   cat_order))),
  data.frame(set = "Shared", table(factor(ranked_shared$distance_category, cat_order))))
write.csv(dist_tbl, file.path(OUT_DIR, "ranked_distance_summary_STRANDED.csv"), row.names = FALSE)

## -----------------------------------------------------------------------------
## Old vs new: what happened to each pair, and to the old top 10 / histone pairs
## -----------------------------------------------------------------------------
old <- read.csv(OLD_CLASSIFIED, stringsAsFactors = FALSE)
cmp <- old[, c("lncRNA", "PCG", "B8441_relationship", "MMC1_relationship", "cross_strain_relationship")] %>%
  rename(PCG_old = PCG, B8441_old = B8441_relationship, MMC1_old = MMC1_relationship,
         cross_old = cross_strain_relationship) %>%
  left_join(expr[, c("lncRNA", "PCG", "overlap_bp", "B8441_relationship", "MMC1_relationship",
                     "cross_strain_relationship")], by = "lncRNA") %>%
  mutate(partner_changed = PCG_old != PCG)
write.csv(cmp, file.path(OUT_DIR, "old_vs_stranded_per_pair.csv"), row.names = FALSE)
cat("\nPartner PCG changed for", sum(cmp$partner_changed), "lncRNAs\n")
print(table(old = cmp$cross_old, stranded = cmp$cross_strain_relationship))

old_top <- c("LNCRNA.1367", "LNCRNA.1431", "LNCRNA.5455", "LNCRNA.2073", "LNCRNA.2826",
             "LNCRNA.5133", "LNCRNA.5135", "LNCRNA.278", "LNCRNA.339", "LNCRNA.5136")
old_top_now <- cmp[match(old_top, cmp$lncRNA), ] %>%
  left_join(ranked_shared[, c("lncRNA", "rank")], by = "lncRNA") %>%
  rename(new_shared_rank = rank)
cat("\nOld top-10 shared candidates under stranded counts (histones: B9J08_002073, _005144, _005145):\n")
print(old_top_now)
write.csv(old_top_now, file.path(OUT_DIR, "old_top10_under_stranded.csv"), row.names = FALSE)

cat("\nDone. Outputs in", OUT_DIR, "\n")
