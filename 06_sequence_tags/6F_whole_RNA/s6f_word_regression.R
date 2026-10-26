## s6f_word_regression.R -- Step 6F (v2, 2026-10-13), stage F5: per-word test, EV-enriched vs all other tested genes, with
## covariates (Justine "V2", 2026-10-13 15:19 UTC = plan Q2 b; Q5 a = GC as 5 classes; Q8 a = >= 3 of 4 groups).
## Per set (filter x group x region), per 6-mer:
##   rate = count * 1000 / region length (one value per gene; each gene counts once);
##   model (rank-based regression, the covariate analogue of the plan's Wilcoxon rank test):
##     rank(rate) ~ enriched + log10(Cell polyA TPM + 1) + log10(transcript length) + GC class (5 quantile classes of the
##     region GC over the set's genes; class = ceil(rank * 5 / n), ties by model_id), ordinary least squares;
##   p = two-sided t test of the "enriched" coefficient; BH over the 4,096 words of the set; call = q < 0.05;
##   direction = sign of the enriched coefficient (adjusted); observed mean rates are reported alongside.
##   Words with the same rate in every gene: p = 1 (kept in the BH, as 6E).
## Random-label null (plan): N random "enriched" sets of the same size drawn from the comparison genes (the real enriched genes
##   removed), the rest of the comparison genes as comparison; calls counted per draw. Run for the FX sets.
## Candidates: per filter x region, q < 0.05 with the same direction in >= 3 of the 4 groups (untestable groups count as no call).
## Usage: Rscript s6f_word_regression.R <F4 dir> <outdir> <real|sim> <seed> <n_null>
a <- commandArgs(TRUE); F4 <- a[1]; OUT <- a[2]; MODE <- a[3]; SEED <- as.integer(a[4]); NNULL <- as.integer(a[5])
stopifnot(MODE %in% c("real", "sim"), !dir.exists(OUT)); dir.create(OUT, recursive = TRUE)
gates <- data.frame(gate = character(), result = character(), detail = character())
gate <- function(n, ok, d = "") gates[nrow(gates) + 1, ] <<- list(n, ifelse(isTRUE(ok), "PASS", "FAIL"), as.character(d))
B <- c("A", "C", "G", "T"); WORDS <- sort(apply(expand.grid(B, B, B, B, B, B, stringsAsFactors = FALSE), 1, paste, collapse = ""))
ss <- read.delim(file.path(F4, "sets_summary.tsv"), colClasses = "character")
olsfit <- function(Y, X, ienr) {               # OLS of every column of Y on X; returns coef, t, p of column ienr
  q <- qr(X); if (q$rank < ncol(X)) return(NULL)
  cf <- qr.coef(q, Y); res <- qr.resid(q, Y); df <- nrow(X) - ncol(X); if (df < 1) return(NULL)
  s2 <- colSums(res^2) / df; XtXi <- chol2inv(qr.R(q)); se <- sqrt(s2 * XtXi[ienr, ienr])
  tt <- cf[ienr, ] / se; p <- 2 * pt(-abs(tt), df)
  const <- apply(Y, 2, function(v) all(v == v[1])); p[const] <- 1; tt[const] <- 0
  list(coef = cf[ienr, ], t = tt, p = p, df = df)
}
SUM <- list(); WT <- list(); set.seed(SEED)
for (i in seq_len(nrow(ss))) {
  s <- ss[i, ]; d <- file.path(F4, "sets", s$set)
  gn <- read.delim(file.path(d, "genes.tsv"), colClasses = c("character", "character", "numeric", "numeric", "numeric", "numeric"))
  kc <- read.delim(gzfile(file.path(d, "kmer_counts.tsv.gz")), colClasses = c("character", "character", "numeric"))
  n <- nrow(gn); nE <- sum(gn$set == "enriched"); nC <- n - nE
  row <- data.frame(set = s$set, filter = s$filter, group = s$group, rna_class = s$rna_class, region = s$region, variant = s$variant,
                    enriched = nE, comparison = nC, df = NA, calls_higher_in_enriched = NA, calls_higher_in_comparison = NA,
                    min_q = NA, null_draws = 0, null_calls_median = NA, null_calls_mean = NA, null_calls_max = NA, null_draws_with_any_call = NA,
                    status = "tested")
  if (nE < 1 || nC < 1) { row$status <- "not_testable_empty_side"; SUM[[i]] <- row; next }
  cnt <- matrix(0, n, length(WORDS), dimnames = list(gn$model_id, WORDS))
  cnt[cbind(match(kc$model_id, gn$model_id), match(kc$word, WORDS))] <- kc$count
  gate(paste0("F5_counts_placed_", s$set), !anyNA(match(kc$word, WORDS)) && !anyNA(match(kc$model_id, gn$model_id)) &&
       sum(cnt) == sum(kc$count) && all(rowSums(cnt) <= pmax(gn$region_length - 5, 0)), sum(kc$count))
  rate <- cnt * 1000 / gn$region_length
  o <- order(gn$region_GC, gn$model_id); gcc <- integer(n); gcc[o] <- ceiling(seq_len(n) * 5 / n)
  enr <- as.numeric(gn$set == "enriched"); lt <- log10(gn$cell_polyA_TPM + 1); ll <- log10(gn$len_whole)
  X <- model.matrix(~ enr + lt + ll + factor(gcc))
  Y <- apply(rate, 2, rank)
  f <- olsfit(Y, X, 2)
  if (is.null(f)) { row$status <- "not_testable_design"; SUM[[i]] <- row; next }
  q <- p.adjust(f$p, "BH")
  # implementation check against lm() on up to 10 non-constant words
  nc <- which(f$p < 1); chkw <- if (length(nc)) nc[unique(round(seq(1, length(nc), length.out = min(10, length(nc)))))] else integer(0)
  agree <- all(vapply(chkw, function(j) { m <- summary(lm(Y[, j] ~ enr + lt + ll + factor(gcc)))$coefficients
    abs(m["enr", 4] - f$p[j]) < 1e-8 * max(1, m["enr", 4]) && abs(m["enr", 1] - f$coef[j]) < 1e-8 * max(1, abs(m["enr", 1])) }, TRUE))
  gate(paste0("F5_matches_lm_", s$set), agree, length(chkw))
  oo <- order(f$p); qq <- rev(cummin(rev(pmin(1, f$p[oo] * 4096 / seq_len(4096))))); q2 <- numeric(4096); q2[oo] <- qq
  gate(paste0("F5_BH_recomputed_p_in_0_1_", s$set), isTRUE(all.equal(q2, unname(q), tolerance = 1e-12)) && all(f$p >= 0 & f$p <= 1))  # r5: unname (r3 failed on names only)
  dirn <- ifelse(f$coef > 0, "higher_in_enriched", ifelse(f$coef < 0, "higher_in_comparison", "none"))
  mE <- colMeans(rate[enr == 1, , drop = FALSE]); mC <- colMeans(rate[enr == 0, , drop = FALSE])
  wt <- data.frame(set = s$set, filter = s$filter, group = s$group, rna_class = s$rna_class, region = s$region, variant = s$variant,
                   word = WORDS, mean_rate_per_kb_enriched = mE, mean_rate_per_kb_comparison = mC,
                   genes_with_word_enriched = colSums(cnt[enr == 1, , drop = FALSE] > 0), genes_with_word_comparison = colSums(cnt[enr == 0, , drop = FALSE] > 0),
                   coef_enriched_rank = f$coef, t = f$t, p = f$p, q = q, direction = dirn,
                   observed_mean_agrees = sign(mE - mC) == sign(f$coef))
  gz <- gzfile(file.path(OUT, paste0(s$set, ".words.tsv.gz")), "w"); write.table(wt, gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz)
  WT[[s$set]] <- wt[wt$q < 0.05, ]
  row$df <- f$df; row$calls_higher_in_enriched <- sum(q < 0.05 & f$coef > 0); row$calls_higher_in_comparison <- sum(q < 0.05 & f$coef < 0)
  row$min_q <- signif(min(q), 4)
  if (MODE == "real" && s$filter == "FX" && nC > nE) {
    ic <- which(enr == 0); nl <- integer(NNULL); Y2 <- apply(rate[ic, , drop = FALSE], 2, rank)
    for (b in seq_len(NNULL)) {
      pick <- sample(ic, nE); e2 <- as.numeric(ic %in% pick)
      X2 <- model.matrix(~ e2 + lt[ic] + ll[ic] + factor(gcc[ic]))
      f2 <- olsfit(Y2, X2, 2); nl[b] <- if (is.null(f2)) NA else sum(p.adjust(f2$p, "BH") < 0.05)
    }
    row$null_draws <- NNULL; row$null_calls_median <- median(nl, na.rm = TRUE); row$null_calls_mean <- mean(nl, na.rm = TRUE)
    row$null_calls_max <- max(nl, na.rm = TRUE); row$null_draws_with_any_call <- sum(nl > 0, na.rm = TRUE)
    gate(paste0("F5_random_label_null_median_0_", s$set), isTRUE(median(nl, na.rm = TRUE) == 0),
         paste0("draws with any call ", sum(nl > 0, na.rm = TRUE), "/", NNULL, "; max ", max(nl, na.rm = TRUE)))
  } else if (MODE == "real" && s$filter == "FX") row$status <- "tested_null_not_run_comparison_le_enriched"
  SUM[[i]] <- row
}
S <- do.call(rbind, SUM); write.table(S, file.path(OUT, "summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
SIG <- do.call(rbind, WT); write.table(SIG, file.path(OUT, "significant_words.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (MODE == "sim") {
  sd <- list(planted_word = readLines(file.path(F4, "sim_planted_word.txt"))[1])
  pl <- S$set[S$variant == "planted"]; gs <- S$set[S$variant == "gcshift"]
  wp <- read.delim(gzfile(file.path(OUT, paste0(pl, ".words.tsv.gz"))))
  r <- wp[wp$word == sd$planted_word, ]
  gate("F5_SIM_planted_word_called_higher_in_enriched", r$q < 0.05 && r$direction == "higher_in_enriched",
       paste0(sd$planted_word, " q=", signif(r$q, 3), " rank ", rank(wp$p)[wp$word == sd$planted_word], "; other calls ", sum(wp$q < 0.05) - (r$q < 0.05)))
  gate("F5_SIM_gcshift_alone_0_calls", S$calls_higher_in_enriched[S$set == gs] + S$calls_higher_in_comparison[S$set == gs] == 0,
       paste0("calls ", S$calls_higher_in_enriched[S$set == gs], " up / ", S$calls_higher_in_comparison[S$set == gs], " down"))
} else {
  CAND <- list()
  for (fl in unique(S$filter)) for (rg in unique(paste(S$rna_class, S$region))) {
    if (is.null(SIG)) next
    sub <- SIG[SIG$filter == fl & paste(SIG$rna_class, SIG$region) == rg, ]
    if (!nrow(sub)) next
    tb <- aggregate(group ~ word + direction, sub, function(x) length(unique(x)))
    tb <- tb[tb$group >= 3, ]
    if (nrow(tb)) CAND[[paste(fl, rg)]] <- data.frame(filter = fl, rna_class_region = rg, word = tb$word, direction = tb$direction, groups_called = tb$group)
  }
  C <- if (length(CAND)) do.call(rbind, CAND) else data.frame(filter = character(), rna_class_region = character(), word = character(), direction = character(), groups_called = integer())
  write.table(C, file.path(OUT, "tag_candidates_ge3of4.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  print(C)
}
write.table(gates, file.path(OUT, "gates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
print(gates[gates$result != "PASS" | grepl("SIM|null", gates$gate), ]); print(S)
ok <- all(gates$result == "PASS"); cat(if (ok) "F5_ALL_GATES_PASS\n" else "F5_GATE_FAILURE\n"); quit(status = if (ok) 0 else 1)
