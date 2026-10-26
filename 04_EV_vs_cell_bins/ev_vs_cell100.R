## ev_vs_cell100.R -- Option-1 100-nt branch, step 5/6 (2026-10-09): direct EV-vs-Cell sRNA comparison of each bin's share of
## its transcript, then combination with the within-compartment hotspot calls. Decisions (Justine, 2026-10-09):
##   test = option C using test B: edgeR replicate-aware differential-usage test, restricted to transcripts with enough reads
##   read rule = transcript >= 10 PHYSICAL reads in >= 2 of 3 EV AND >= 2 of 3 Cell sRNA libraries of the same group (D6 value)
##   layout = the layout of the hotspot results it is combined with (main: merged/N; sensitivity: unmerged/U)
##   comparisons = EV vs Cell sRNA within each strain x treatment group, never pooled
##   no FDR level fixed: q-values reported; counts at 0.01/0.05/0.10 for reference only (Justine chooses after review)
## Implementation:
##   counts = bin_count (incl. Caro's neighbouring-tie +1s) from the tie-rule tables; excluded bins (D5: competitor-only rDNA
##   lncRNA bins, lncRNA bins overlapping another same-strand transcript) removed before testing; bins with 0 reads in all 6
##   libraries not tested; transcripts with < 2 remaining bins not tested (a usage test needs >= 2 features).
##   edgeR: DGEList (lib.size = column sums over all non-excluded bins), design ~ group (Cell reference), glmQLFit (defaults:
##   legacy = FALSE, robust = FALSE), diffSplice.DGEGLM(coef = "groupEV") = the edgeR-authors' recommended replacement for
##   diffSpliceDGE (installed edgeR 4.8.2 docs; method of Baldoni et al. 2025). Legacy diffSpliceDGE run alongside as an
##   agreement check only.
##   q_dataset = BH over all tested bins of the group comparison (PCG + lncRNA together); q_within_transcript = BH within
##   transcript (ranking). Transcript-level: Simes p (edgeR) + BH.
##   Observed shares: per library bin_count / transcript physical reads (as D2); libraries with 0 transcript reads excluded
##   from the mean (G1); fold enrichment = mean EV share / mean Cell share; DIRECTION FROM OBSERVED SHARES.
## Usage: Rscript ev_vs_cell100.R <EV_tables_dir> <Cell_tables_dir> <annotation.tsv> <EV_manifest> <Cell_manifest>
##        <EV_hotspot_dir> <Cell_hotspot_dir> <outdir> <real|null_sim> <seed>
## r2 (2026-10-09): shares equal to within 1e-12 (relative) are labelled "equal" with fold enrichment 1 (r1 let float
##   rounding of mathematically equal shares decide the direction; found by the independent check, 7 untested bins).
## null_sim: per transcript, counts redrawn per library from the POOLED share vector (multinomial, same per-library totals),
##   so EV and Cell have no true difference; the same pipeline is run (calibration check of the pipeline, informational).
suppressMessages({ library(limma); library(edgeR) })
a <- commandArgs(TRUE)
EVD <- a[1]; CD <- a[2]; ANN <- a[3]; MEV <- a[4]; MC <- a[5]; HEV <- a[6]; HC <- a[7]; OUT <- a[8]; MODE <- a[9]; seed <- as.integer(a[10])
stopifnot(MODE %in% c("real", "null_sim"), !dir.exists(OUT)); dir.create(OUT, recursive = TRUE)
MINREADS <- 10; MINREPS <- 2; LEVELS <- c(0.01, 0.05, 0.10)
ann <- read.delim(ANN, colClasses = "character")
excl_rdna <- ann$bin_id[ann$role == "competitor_only"]
excl_ovl <- ann$bin_id[ann$same_strand_binned_model_overlap != "" & ann$role != "competitor_only"]
mev <- read.delim(MEV, colClasses = "character"); mc <- read.delim(MC, colClasses = "character")
gates <- data.frame(gate = character(), result = character(), detail = character())
gate <- function(n, ok, d = "") gates[nrow(gates) + 1, ] <<- list(n, ifelse(isTRUE(ok), "PASS", "FAIL"), as.character(d))
rd <- function(D, L) read.delim(gzfile(file.path(D, paste0(L, ".bin100_final_table.tsv.gz"))), colClasses = "character")
set.seed(seed)
SUM <- list(); CAT <- list()
for (g in unique(mev$group)) {
  le <- mev$library[mev$group == g]; lc <- mc$library[mc$group == g]
  stopifnot(length(le) == 3, length(lc) == 3)
  libs <- c(lc, le)                                                   # Cell first (reference), then EV
  Tb <- c(lapply(lc, function(L) rd(CD, L)), lapply(le, function(L) rd(EVD, L)))
  base <- Tb[[1]][, c("transcript_id", "gene_id", "rna_class", "role", "bin_id", "bin_number", "n_bins", "bin_tx_start0", "bin_tx_end", "bin_length")]
  gate(paste0("E1_same_bins_all_6_libraries_", g), all(vapply(Tb, function(t) identical(t$bin_id, base$bin_id), TRUE)) &&
       setequal(base$bin_id, ann$bin_id), nrow(base))
  cnt <- sapply(Tb, function(t) as.numeric(t$bin_count)); colnames(cnt) <- libs
  txr <- sapply(Tb, function(t) as.numeric(t$transcript_physical_reads)); colnames(txr) <- libs
  if (MODE == "null_sim") {                                           # redraw from pooled shares: no true EV-Cell difference
    for (tx in unique(base$transcript_id)) {
      r <- which(base$transcript_id == tx); tot <- colSums(cnt[r, , drop = FALSE]); p <- rowSums(cnt[r, , drop = FALSE])
      if (sum(p) == 0) next
      for (j in seq_along(libs)) if (tot[j] > 0) cnt[r, j] <- rmultinom(1, tot[j], p / sum(p))[, 1]
    }
  }
  # ---- status (read rule on physical transcript reads; same values for all bins of a transcript)
  ev_ok <- rowSums(txr[, le, drop = FALSE] >= MINREADS) >= MINREPS
  ce_ok <- rowSums(txr[, lc, drop = FALSE] >= MINREADS) >= MINREPS
  status <- ifelse(base$bin_id %in% excl_rdna, "excluded_rDNA_competitor_only",
            ifelse(base$bin_id %in% excl_ovl, "excluded_same_strand_overlap",
            ifelse(!(ev_ok & ce_ok), "not_tested_read_rule",
            ifelse(rowSums(cnt) == 0, "not_tested_zero_in_all_6", "candidate"))))
  ncand <- ave(status == "candidate", base$transcript_id, FUN = sum)
  status[status == "candidate" & ncand < 2] <- "not_tested_lt2_bins"
  status[status == "candidate"] <- "tested"
  te <- status == "tested"
  # ---- observed shares and fold enrichment (direction from these)
  share <- ifelse(txr > 0, cnt / txr, NA)
  mE <- rowMeans(share[, le, drop = FALSE], na.rm = TRUE); mC <- rowMeans(share[, lc, drop = FALSE], na.rm = TRUE)
  mE[is.nan(mE)] <- NA; mC[is.nan(mC)] <- NA
  fe <- ifelse(!is.na(mE) & !is.na(mC) & mC > 0, mE / mC, ifelse(!is.na(mE) & !is.na(mC) & mE > 0 & mC == 0, Inf, NA))
  eqtol <- !is.na(mE) & !is.na(mC) & abs(mE - mC) <= 1e-12 * pmax(mE, mC)   # r2: float-rounding ties are equal
  dirn <- ifelse(is.na(mE) | is.na(mC), NA, ifelse(eqtol, "equal", ifelse(mE > mC, "higher_in_EV", "higher_in_Cell")))
  fe[eqtol & mC > 0] <- 1
  # ---- edgeR fit and differential-usage test
  keep_lib <- !(base$bin_id %in% c(excl_rdna, excl_ovl))
  grp <- factor(c(rep("Cell", 3), rep("EV", 3)), levels = c("Cell", "EV")); design <- model.matrix(~ grp)
  colnames(design) <- c("(Intercept)", "groupEV")
  y <- DGEList(counts = cnt[te, , drop = FALSE], lib.size = colSums(cnt[keep_lib, , drop = FALSE]),
               genes = data.frame(transcript_id = base$transcript_id[te], bin_id = base$bin_id[te]))
  fit <- glmQLFit(y, design)
  ds <- diffSplice(fit, coef = "groupEV", geneid = "transcript_id", exonid = "bin_id", verbose = FALSE)
  m <- match(base$bin_id, ds$genes$bin_id)
  col1 <- function(x) as.numeric(as.matrix(x)[, 1])
  coefU <- col1(ds$coefficients)[m]; tU <- col1(ds$t)[m]; p <- col1(ds$p.value)[m]
  p[!te] <- NA; coefU[!te] <- NA; tU[!te] <- NA
  q <- rep(NA_real_, length(p)); q[te] <- p.adjust(p[te], method = "BH")
  qw <- rep(NA_real_, length(p)); qw[te] <- ave(p[te], base$transcript_id[te], FUN = function(x) p.adjust(x, method = "BH"))
  # legacy function, agreement check only
  dl <- diffSpliceDGE(fit, coef = 2, geneid = "transcript_id", exonid = "bin_id", verbose = FALSE)
  pl <- dl$exon.p.value[match(base$bin_id, dl$genes$bin_id)]; pl[!te] <- NA
  # ---- gates
  gate(paste0("E2_tested_bins_once_in_output_", g), sum(te) == nrow(ds$genes) && !anyNA(m[te]) && !anyDuplicated(ds$genes$bin_id), sum(te))
  gate(paste0("E3_tested_transcripts_ge2_bins_nonzero_", g), all(rowSums(cnt[te, , drop = FALSE]) > 0) &&
       all(table(base$transcript_id[te]) >= 2), length(unique(base$transcript_id[te])))
  gate(paste0("E4_excluded_and_untested_have_no_p_", g), all(is.na(p[!te])) && !any(te & base$bin_id %in% c(excl_rdna, excl_ovl)), sum(!te))
  gate(paste0("E5_p_in_0_1_", g), all(p[te] >= 0 & p[te] <= 1) && !anyNA(p[te]))
  o <- order(p[te]); pp <- p[te][o]; mm <- length(pp); qq <- rev(cummin(rev(pmin(1, pp * mm / seq_len(mm))))); q2 <- numeric(mm); q2[o] <- qq
  gate(paste0("E6_BH_recomputed_", g), isTRUE(all.equal(q2, q[te], tolerance = 1e-12)), mm)
  sg <- sign(coefU[te]); so <- ifelse(dirn[te] == "higher_in_EV", 1, ifelse(dirn[te] == "higher_in_Cell", -1, 0))
  # ---- combine with hotspot calls (step 6)
  out <- data.frame(group = g, base, status = status)
  for (j in seq_along(libs)) out[[paste0("count_", libs[j])]] <- cnt[, j]
  for (j in seq_along(libs)) out[[paste0("transcript_reads_", libs[j])]] <- txr[, j]
  out$mean_share_EV <- mE; out$mean_share_Cell <- mC; out$fold_enrichment_EV_over_Cell <- fe; out$direction_observed <- dirn
  out$usage_coef_EV_vs_Cell <- coefU; out$usage_t <- tU; out$p <- p; out$q_dataset <- q; out$q_within_transcript <- qw
  out$p_legacy_diffSpliceDGE <- pl
  rdh <- function(D) { f <- file.path(D, paste0(g, ".group_calls.tsv.gz")); if (file.exists(f)) read.delim(gzfile(f), colClasses = "character") else NULL }
  he <- rdh(HEV); hc <- rdh(HC)
  if (!is.null(he) && !is.null(hc)) {
    gate(paste0("E7_hotspot_bins_match_layout_", g), setequal(he$bin_id, base$bin_id) && setequal(hc$bin_id, base$bin_id))
    for (lv in LEVELS) {
      evh <- he[[paste0("group_call_2of3_q", lv)]][match(base$bin_id, he$bin_id)] == "TRUE"
      ceh <- hc[[paste0("group_call_2of3_q", lv)]][match(base$bin_id, hc$bin_id)] == "TRUE"
      sig <- te & !is.na(q) & q < lv
      cat_ <- ifelse(evh,
                ifelse(!te, "EV_hotspot__comparison_not_tested",
                ifelse(sig & dirn == "higher_in_EV", ifelse(ceh, "EV_hotspot__also_Cell_hotspot__EV_share_higher", "EV_hotspot__EV_share_higher"),
                ifelse(sig & dirn == "higher_in_Cell", "EV_hotspot__Cell_share_higher",
                ifelse(ceh, "EV_hotspot__also_Cell_hotspot__no_significant_difference", "EV_hotspot__no_Cell_hotspot__no_significant_difference")))),
              ifelse(ceh & sig & dirn == "higher_in_Cell", "Cell_hotspot__Cell_share_higher",
              ifelse(ceh, "Cell_hotspot__other",
              ifelse(sig, paste0("no_hotspot__share_", dirn), "no_hotspot__no_difference"))))
      out[[paste0("EV_hotspot_q", lv)]] <- evh; out[[paste0("Cell_hotspot_q", lv)]] <- ceh; out[[paste0("category_q", lv)]] <- cat_
      tb <- table(cat_); CAT[[paste(g, lv)]] <- data.frame(group = g, level = lv, category = names(tb), bins = as.integer(tb))
    }
  }
  gz <- gzfile(file.path(OUT, paste0(g, ".ev_vs_cell_bins.tsv.gz")), "w"); write.table(out, gz, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"); close(gz)
  nb <- table(base$transcript_id[te])
  tt <- data.frame(group = g, transcript_id = ds$gene.genes$transcript_id,
                   n_bins_tested = as.integer(nb[ds$gene.genes$transcript_id]),
                   gene_F_p = col1(ds$gene.F.p.value), simes_p = col1(ds$gene.simes.p.value))
  tt$simes_q <- p.adjust(tt$simes_p, method = "BH")
  gz <- gzfile(file.path(OUT, paste0(g, ".ev_vs_cell_transcripts.tsv.gz")), "w"); write.table(tt, gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz)
  st <- table(factor(status, levels = c("tested", "not_tested_read_rule", "not_tested_zero_in_all_6", "not_tested_lt2_bins",
                                        "excluded_rDNA_competitor_only", "excluded_same_strand_overlap")))
  ph <- hist(p[te], breaks = seq(0, 1, 0.05), plot = FALSE)$counts
  s <- data.frame(group = g, mode = MODE, transcripts_tested = length(unique(base$transcript_id[te])), bins_tested = st[["tested"]],
                  bins_not_tested_read_rule = st[["not_tested_read_rule"]], bins_not_tested_zero = st[["not_tested_zero_in_all_6"]],
                  bins_not_tested_lt2_bins = st[["not_tested_lt2_bins"]], bins_excluded_rDNA = st[["excluded_rDNA_competitor_only"]],
                  bins_excluded_overlap = st[["excluded_same_strand_overlap"]],
                  p_lt_0.05_fraction = round(mean(p[te] < 0.05), 4), p_histogram_20bins = paste(ph, collapse = ","),
                  sign_coef_agrees_observed_direction = round(mean(sg[so != 0] == so[so != 0]), 4),
                  spearman_p_new_vs_legacy = round(cor(p[te], pl[te], method = "spearman"), 4))
  for (lv in LEVELS) {
    s[[paste0("calls_q", lv, "_higher_in_EV")]] <- sum(te & q < lv & dirn == "higher_in_EV", na.rm = TRUE)
    s[[paste0("calls_q", lv, "_higher_in_Cell")]] <- sum(te & q < lv & dirn == "higher_in_Cell", na.rm = TRUE)
  }
  SUM[[g]] <- s
}
S <- do.call(rbind, SUM); write.table(S, file.path(OUT, "summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (length(CAT)) write.table(do.call(rbind, CAT), file.path(OUT, "category_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gates, file.path(OUT, "gates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
print(gates); print(S)
ok <- all(gates$result == "PASS"); cat(if (ok) "E_ALL_GATES_PASS\n" else "E_GATE_FAILURE\n"); quit(status = if (ok) 0 else 1)
