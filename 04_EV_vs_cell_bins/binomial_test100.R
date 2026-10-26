## binomial_test100.R -- Option-1 100-nt branch, one-sided exact binomial test (2026-10-08). Specification = decision
## register entry 2026-10-08 (second) / decision sheet D1-D8 (Justine):
##   per bin per library: p = P(X >= k), X ~ Binomial(n, bin_length / transcript_length)            [read counts]
##   k = bin count incl. tie +1 (D2); n = transcript PHYSICAL reads in that library (D2); merged final bin (D1, annotation)
##   tested only if n >= 10 (D6, transcript total only); excluded: competitor-only rDNA lncRNAs and lncRNA bins overlapping
##   another same-strand transcript (D5) -- excluded from testing and from both BH families
##   BH dataset-wide per library over all tested bins, PCG + lncRNA together (calling) + BH within transcript x library
##   (ranking) (D3); no FDR level fixed (D8): q-values reported, call counts at 0.01/0.05/0.10 for reference only;
##   group call = called in >= 2 of 3 replicates, untested replicate = not called (D4) -- reference levels only.
## Usage: Rscript binomial_test100.R <final_tables_dir> <merged_annotation.tsv> <manifest.tsv> <outdir> <seed>
## r2 (2026-10-09): works with either bin layout (merged_final_bin column optional); adds final_bin = last bin of its
##   transcript (bin_number == n_bins) to results and group calls; group_summary adds final_bins_among_calls_q0.05.
## Validation inside: binom.test() == pbinom() on a seeded random sample of tested bins per library; p.adjust BH recomputed.
a <- commandArgs(TRUE); TD <- a[1]; ANN <- a[2]; MAN <- a[3]; OUT <- a[4]; seed <- as.integer(a[5])
stopifnot(!dir.exists(OUT)); dir.create(OUT, recursive = TRUE)
LEVELS <- c(0.01, 0.05, 0.10); MINREADS <- 10
ann <- read.delim(ANN, colClasses = "character")
excl_overlap <- ann$bin_id[ann$same_strand_binned_model_overlap != ""]
man <- read.delim(MAN, colClasses = "character")
gates <- data.frame(gate = character(), result = character(), detail = character())
gate <- function(n, ok, d = "") gates[nrow(gates) + 1, ] <<- list(n, ifelse(isTRUE(ok), "PASS", "FAIL"), as.character(d))
set.seed(seed)
summ <- list(); per <- list()
for (i in seq_len(nrow(man))) {
  L <- man$library[i]
  t <- read.delim(gzfile(file.path(TD, paste0(L, ".bin100_final_table.tsv.gz"))), colClasses = "character")
  num <- c("bin_number", "n_bins", "bin_tx_start0", "bin_tx_end", "bin_length", "bin_count", "transcript_physical_reads")
  for (c in num) t[[c]] <- as.numeric(t[[c]])
  tl <- ann$transcript_length[match(t$bin_id, ann$bin_id)]
  stopifnot(!anyNA(tl)); t$transcript_length <- as.numeric(tl)
  t$status <- ifelse(t$role == "competitor_only", "excluded_rDNA_competitor_only",
              ifelse(t$bin_id %in% excl_overlap, "excluded_same_strand_overlap",
              ifelse(t$transcript_physical_reads == 0, "not_tested_no_reads",
              ifelse(t$transcript_physical_reads < MINREADS, "not_tested_lt10_reads", "tested"))))
  te <- t$status == "tested"
  t$k <- t$bin_count; t$n <- t$transcript_physical_reads; t$p0 <- t$bin_length / t$transcript_length
  t$expected <- t$n * t$p0; t$obs_over_exp <- ifelse(t$expected > 0, t$k / t$expected, NA)
  t$p <- NA_real_; t$q_dataset <- NA_real_; t$q_within_transcript <- NA_real_
  t$p[te] <- pbinom(t$k[te] - 1, t$n[te], t$p0[te], lower.tail = FALSE)
  t$q_dataset[te] <- p.adjust(t$p[te], method = "BH")
  t$q_within_transcript[te] <- ave(t$p[te], t$transcript_id[te], FUN = function(x) p.adjust(x, method = "BH"))
  # gates
  gate(paste0("S1_k_le_n_and_p_in_0_1_", L), all(t$k[te] <= t$n[te]) && all(t$p[te] >= 0 & t$p[te] <= 1), sum(te))
  smp <- which(te); smp <- smp[sample.int(length(smp), min(2000, length(smp)))]
  bt <- vapply(smp, function(j) binom.test(t$k[j], t$n[j], t$p0[j], alternative = "greater")$p.value, 0)
  gate(paste0("S2_binom.test_eq_pbinom_", L), isTRUE(all.equal(bt, t$p[smp], tolerance = 1e-10)), paste("n_checked", length(smp), "max_abs_diff", signif(max(abs(bt - t$p[smp])), 3)))
  o <- order(t$p[te]); pp <- t$p[te][o]; m <- length(pp)
  qq <- rev(cummin(rev(pmin(1, pp * m / seq_len(m))))); q2 <- numeric(m); q2[o] <- qq
  gate(paste0("S3_BH_dataset_recomputed_", L), isTRUE(all.equal(q2, t$q_dataset[te], tolerance = 1e-12)), m)
  gate(paste0("S4_excluded_bins_untested_", L), all(is.na(t$p[!te])), sum(!te))
  keep <- c("library", "strain", "treatment", "replicate", "transcript_id", "gene_id", "rna_class", "bin_id", "bin_number", "n_bins",
            "bin_tx_start0", "bin_tx_end", "bin_length", "transcript_length", "status", "k", "n", "p0", "expected", "obs_over_exp",
            "p", "q_dataset", "q_within_transcript")
  t$merged_final_bin <- if ("merged_final_bin" %in% names(ann)) ann$merged_final_bin[match(t$bin_id, ann$bin_id)] else ""
  t$final_bin <- ifelse(t$bin_number == t$n_bins, "yes", "")
  out <- t[, c(keep, "merged_final_bin", "final_bin")]
  gz <- gzfile(file.path(OUT, paste0(L, ".binomial_results.tsv.gz")), "w"); write.table(out, gz, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"); close(gz)
  per[[L]] <- out
  st <- table(factor(t$status, levels = c("tested", "not_tested_lt10_reads", "not_tested_no_reads", "excluded_rDNA_competitor_only", "excluded_same_strand_overlap")))
  tx <- tapply(t$status == "tested", t$transcript_id, any)
  ph <- hist(t$p[te], breaks = seq(0, 1, 0.05), plot = FALSE)$counts
  summ[[L]] <- data.frame(library = L, group = man$group[i], replicate = man$replicate[i], transcripts_tested = sum(tx),
                          bins_tested = st[["tested"]], bins_not_tested_lt10 = st[["not_tested_lt10_reads"]],
                          bins_not_tested_no_reads = st[["not_tested_no_reads"]], bins_excluded_rDNA = st[["excluded_rDNA_competitor_only"]],
                          bins_excluded_overlap = st[["excluded_same_strand_overlap"]],
                          calls_q_0.01 = sum(t$q_dataset < 0.01, na.rm = TRUE), calls_q_0.05 = sum(t$q_dataset < 0.05, na.rm = TRUE),
                          calls_q_0.10 = sum(t$q_dataset < 0.10, na.rm = TRUE),
                          p_lt_0.05_fraction = round(mean(t$p[te] < 0.05), 4), p_histogram_20bins = paste(ph, collapse = ","))
}
S <- do.call(rbind, summ); write.table(S, file.path(OUT, "library_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# group calls (D4): called in >= 2 of 3 replicates at each reference level; untested replicate = not called
G <- list()
for (g in unique(man$group)) {
  libs <- man$library[man$group == g]; stopifnot(length(libs) == 3)
  base <- per[[libs[1]]][, c("transcript_id", "rna_class", "bin_id", "bin_number", "n_bins", "bin_length", "merged_final_bin", "final_bin")]
  for (lv in LEVELS) {
    called <- sapply(libs, function(L) { q <- per[[L]]$q_dataset[match(base$bin_id, per[[L]]$bin_id)]; !is.na(q) & q < lv })
    base[[paste0("n_replicates_called_q", lv)]] <- rowSums(called)
    base[[paste0("group_call_2of3_q", lv)]] <- rowSums(called) >= 2
  }
  base$n_replicates_tested <- rowSums(sapply(libs, function(L) per[[L]]$status[match(base$bin_id, per[[L]]$bin_id)] == "tested"))
  gz <- gzfile(file.path(OUT, paste0(g, ".group_calls.tsv.gz")), "w"); write.table(base, gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz)
  G[[g]] <- data.frame(group = g, bins_tested_in_ge2_reps = sum(base$n_replicates_tested >= 2),
                       group_calls_2of3_q0.01 = sum(base$group_call_2of3_q0.01), group_calls_2of3_q0.05 = sum(base$group_call_2of3_q0.05),
                       group_calls_2of3_q0.1 = sum(base$group_call_2of3_q0.1),
                       transcripts_with_group_call_q0.05 = length(unique(base$transcript_id[base$group_call_2of3_q0.05])),
                       merged_final_bins_among_calls_q0.05 = sum(base$group_call_2of3_q0.05 & base$merged_final_bin == "yes"),
                       final_bins_among_calls_q0.05 = sum(base$group_call_2of3_q0.05 & base$final_bin == "yes"),
                       final_bins_among_bins_tested_ge2 = sum(base$n_replicates_tested >= 2 & base$final_bin == "yes"))
}
write.table(do.call(rbind, G), file.path(OUT, "group_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gates, file.path(OUT, "gates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
print(gates); print(S[, 1:13]); print(do.call(rbind, G))
cat(if (all(gates$result == "PASS")) "M5_ALL_GATES_PASS\n" else "M5_GATE_FAILURE\n"); quit(status = if (all(gates$result == "PASS")) 0 else 1)
