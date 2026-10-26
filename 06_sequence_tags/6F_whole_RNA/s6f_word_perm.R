## s6f_word_perm.R -- Step 6F (v2), stage F5 r6 (2026-10-14): per-set word test with PERMUTATION p-values (Justine, option a,
## 2026-10-14 07:28 UTC), replacing the t-test p of r5, whose random-label null failed for sparse UTR words.
## Model (unchanged from r5, Justine "V2"): rank(rate per kb) ~ enriched + log10(Cell polyA TPM + 1) + log10(transcript length)
##   + GC class (5 quantile classes of region GC; class = ceil(rank * 5 / n), ties by model_id); statistic = t of "enriched".
## Permutation p (Freedman & Lane 1983, the standard permutation scheme with covariates; Anderson & Legendre 1999; Winkler et al.
##   2014): reduced model rank ~ covariates -> fitted F0 + residuals R0; each permutation shuffles the rows of R0 (one shuffle shared
##   by all words), Y* = F0 + R0[perm, ], t* from the full model; exceedance = |t*| >= |t_obs|.
## Sequential stopping (Besag & Clifford 1991): a word stops after H exceedances, p = H / L (L = permutations used); words that
##   never reach H by BMAX permutations get p = (g + 1) / (BMAX + 1). Computing settings (not calling thresholds): H = 10,
##   BMAX = 100,000 (smallest attainable p ~1e-5, needed for BH over 4,096 words at q < 0.05). BH over the 4,096 words; call =
##   q < 0.05; direction = sign of the adjusted coefficient. The r5 t-test p is kept as a column (p_ttest) for comparison.
## Random-label null (plan): NNULL random "enriched" sets of the enriched-set size from the comparison genes (real enriched genes
##   removed), same permutation test; calls counted per draw; gate = median calls 0 (FX real sets only).
## r7 (2026-10-14, Justine 07:32 via the project assistant, after Caro's feedback): COV = GC (abundance + length + GC class) or
##   noGC (abundance + length only); both reported side by side. GC correction is an additional control (EV RNA may be GC-rich).
## Usage: Rscript s6f_word_perm.R <F4 dir> <set> <outdir> <seed> <n_null> <BMAX> <H> <GC|noGC>
a <- commandArgs(TRUE); F4 <- a[1]; SET <- a[2]; OUT <- a[3]; SEED <- as.integer(a[4]); NNULL <- as.integer(a[5])
BMAX <- as.integer(a[6]); H <- as.integer(a[7]); COV <- a[8]; stopifnot(COV %in% c("GC", "noGC"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE); stopifnot(!file.exists(file.path(OUT, paste0(SET, ".summary.tsv"))))
gates <- data.frame(gate = character(), result = character(), detail = character())
gate <- function(n, ok, d = "") gates[nrow(gates) + 1, ] <<- list(n, ifelse(isTRUE(ok), "PASS", "FAIL"), as.character(d))
B <- c("A", "C", "G", "T"); WORDS <- sort(apply(expand.grid(B, B, B, B, B, B, stringsAsFactors = FALSE), 1, paste, collapse = ""))
ss <- read.delim(file.path(F4, "sets_summary.tsv"), colClasses = "character"); s <- ss[ss$set == SET, ]; stopifnot(nrow(s) == 1)
flp <- function(Y, X, X0, seed) {      # Freedman-Lane sequential permutation test of column 2 of X, every column of Y
  q <- qr(X); q0 <- qr(X0); df <- nrow(X) - ncol(X)
  if (q$rank < ncol(X) || df < 1) return(NULL)
  vi <- chol2inv(qr.R(q))[2, 2]
  tstat <- function(Yc) { cf <- qr.coef(q, Yc); r <- qr.resid(q, Yc); cf[2, ] / sqrt(colSums(r^2) / df * vi) }
  const <- apply(Y, 2, function(v) all(v == v[1]))
  t0 <- rep(0, ncol(Y)); t0[!const] <- tstat(Y[, !const, drop = FALSE]); const <- const | is.nan(t0); t0[is.nan(t0)] <- 0
  cf0 <- rep(0, ncol(Y)); cf0[!const] <- qr.coef(q, Y[, !const, drop = FALSE])[2, ]
  pt_ <- ifelse(const, 1, 2 * pt(-abs(t0), df))
  F0 <- qr.fitted(q0, Y); R0 <- Y - F0
  ex <- integer(ncol(Y)); L <- integer(ncol(Y)); act <- which(!const); b <- 0; set.seed(seed)
  while (length(act) && b < BMAX) {
    b <- b + 1; pm <- sample.int(nrow(Y))
    ts <- tstat(F0[, act, drop = FALSE] + R0[pm, act, drop = FALSE])
    hit <- abs(ts) >= abs(t0[act]) * (1 - 1e-9); hit[is.na(hit)] <- FALSE
    ex[act] <- ex[act] + hit; L[act] <- L[act] + 1L; act <- act[ex[act] < H]
  }
  p <- ifelse(const, 1, ifelse(ex >= H, H / L, (ex + 1) / (L + 1)))
  list(t = t0, coef = cf0, p = p, p_ttest = pt_, ex = ex, L = L, df = df, F0 = F0, R0 = R0, const = const)
}
d <- file.path(F4, "sets", SET)
gn <- read.delim(file.path(d, "genes.tsv"), colClasses = c("character", "character", "numeric", "numeric", "numeric", "numeric"))
kc <- read.delim(gzfile(file.path(d, "kmer_counts.tsv.gz")), colClasses = c("character", "character", "numeric"))
n <- nrow(gn); nE <- sum(gn$set == "enriched"); nC <- n - nE
row <- data.frame(set = SET, covariates = COV, filter = s$filter, group = s$group, rna_class = s$rna_class, region = s$region, variant = s$variant,
                  enriched = nE, comparison = nC, df = NA, calls_higher_in_enriched = NA, calls_higher_in_comparison = NA, min_q = NA,
                  words_reaching_BMAX = NA, perms_max_used = NA, ttest_calls_r5_style = NA,
                  null_draws = 0, null_calls_median = NA, null_calls_mean = NA, null_calls_max = NA, null_draws_with_any_call = NA,
                  null_calls_all = "", status = "tested")
fin <- function() {
  write.table(row, file.path(OUT, paste0(SET, ".summary.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(gates, file.path(OUT, paste0(SET, ".gates.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  print(gates); print(t(row)); ok <- all(gates$result == "PASS"); cat(if (ok) "F5SET_PASS\n" else "F5SET_GATE_FAILURE\n")
  quit(status = 0)                      # gate results are judged by the aggregator; the array task itself succeeds
}
if (nE < 1 || nC < 1) { row$status <- "not_testable_empty_side"; fin() }
cnt <- matrix(0, n, length(WORDS), dimnames = list(gn$model_id, WORDS))
cnt[cbind(match(kc$model_id, gn$model_id), match(kc$word, WORDS))] <- kc$count
gate(paste0("F5_counts_placed_", SET), !anyNA(match(kc$word, WORDS)) && !anyNA(match(kc$model_id, gn$model_id)) &&
     sum(cnt) == sum(kc$count) && all(rowSums(cnt) <= pmax(gn$region_length - 5, 0)), sum(kc$count))
rate <- cnt * 1000 / gn$region_length
o <- order(gn$region_GC, gn$model_id); gcc <- integer(n); gcc[o] <- ceiling(seq_len(n) * 5 / n)
enr <- as.numeric(gn$set == "enriched"); lt <- log10(gn$cell_polyA_TPM + 1); ll <- log10(gn$len_whole)
X <- if (COV == "GC") model.matrix(~ enr + lt + ll + factor(gcc)) else model.matrix(~ enr + lt + ll); X0 <- X[, -2, drop = FALSE]
Y <- apply(rate, 2, rank)
f <- flp(Y, X, X0, SEED)
if (is.null(f)) { row$status <- "not_testable_design"; fin() }
q <- p.adjust(f$p, "BH")
# check 1: observed t / coefficient == lm() for up to 10 non-constant words
nc <- which(!f$const); chk <- if (length(nc)) nc[unique(round(seq(1, length(nc), length.out = min(10, length(nc)))))] else integer(0)
fm <- function(y) if (COV == "GC") lm(y ~ enr + lt + ll + factor(gcc)) else lm(y ~ enr + lt + ll)
agree <- all(vapply(chk, function(j) { m <- summary(fm(Y[, j]))$coefficients
  abs(m["enr", 3] - f$t[j]) < 1e-8 * max(1, abs(m["enr", 3])) && abs(m["enr", 1] - f$coef[j]) < 1e-8 * max(1, abs(m["enr", 1])) }, TRUE))
gate(paste0("F5_observed_t_matches_lm_", SET), agree, length(chk))
# check 2: brute-force Freedman-Lane with lm() and the same permutation sequence reproduces exceedance counts / stop points
bf <- which(!f$const & f$L <= 300); bf <- bf[unique(round(seq(1, length(bf), length.out = min(3, length(bf)))))]
okbf <- TRUE
if (length(bf)) {
  set.seed(SEED); Lmax <- max(f$L[bf]); exb <- integer(length(bf)); Lb <- integer(length(bf)); live <- rep(TRUE, length(bf))
  for (b in seq_len(Lmax)) {
    pm <- sample.int(n)
    for (k in which(live)) { j <- bf[k]; ys <- f$F0[, j] + f$R0[pm, j]
      tb <- summary(fm(ys))$coefficients["enr", 3]
      exb[k] <- exb[k] + (abs(tb) >= abs(f$t[j]) * (1 - 1e-9)); Lb[k] <- Lb[k] + 1L; if (exb[k] >= H) live[k] <- FALSE }
  }
  okbf <- all(exb == f$ex[bf] & Lb == f$L[bf])
}
gate(paste0("F5_sequential_perm_matches_bruteforce_lm_", SET), okbf, paste0(length(bf), " words; L = ", paste(f$L[bf], collapse = ",")))
oo <- order(f$p); qq <- rev(cummin(rev(pmin(1, f$p[oo] * 4096 / seq_len(4096))))); q2 <- numeric(4096); q2[oo] <- qq
gate(paste0("F5_BH_recomputed_p_in_0_1_", SET), isTRUE(all.equal(q2, unname(q), tolerance = 1e-12)) && all(f$p > 0 & f$p <= 1))
dirn <- ifelse(f$coef > 0, "higher_in_enriched", ifelse(f$coef < 0, "higher_in_comparison", "none"))
mE <- colMeans(rate[enr == 1, , drop = FALSE]); mC <- colMeans(rate[enr == 0, , drop = FALSE])
wt <- data.frame(set = SET, covariates = COV, filter = s$filter, group = s$group, rna_class = s$rna_class, region = s$region, variant = s$variant,
                 word = WORDS, mean_rate_per_kb_enriched = mE, mean_rate_per_kb_comparison = mC,
                 genes_with_word_enriched = colSums(cnt[enr == 1, , drop = FALSE] > 0), genes_with_word_comparison = colSums(cnt[enr == 0, , drop = FALSE] > 0),
                 coef_enriched_rank = f$coef, t = f$t, perm_exceedances = f$ex, perms_used = f$L, p = f$p, q = q, direction = dirn,
                 p_ttest = f$p_ttest, q_ttest = p.adjust(f$p_ttest, "BH"), observed_mean_agrees = sign(mE - mC) == sign(f$coef))
gz <- gzfile(file.path(OUT, paste0(SET, ".words.tsv.gz")), "w"); write.table(wt, gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz)
row$df <- f$df; row$calls_higher_in_enriched <- sum(q < 0.05 & f$coef > 0); row$calls_higher_in_comparison <- sum(q < 0.05 & f$coef < 0)
row$min_q <- signif(min(q), 4); row$words_reaching_BMAX <- sum(f$L >= BMAX); row$perms_max_used <- max(f$L)
row$ttest_calls_r5_style <- sum(wt$q_ttest < 0.05)
if (s$variant == "real" && s$filter == "FX") {
  if (nC > nE) {
    ic <- which(enr == 0); nl <- integer(NNULL); Y2 <- apply(rate[ic, , drop = FALSE], 2, rank); set.seed(SEED + 1e6)
    seeds <- sample.int(1e8, NNULL); picks <- lapply(seq_len(NNULL), function(b) sample(length(ic), nE))
    for (b in seq_len(NNULL)) {
      e2 <- as.numeric(seq_along(ic) %in% picks[[b]])
      X2 <- if (COV == "GC") model.matrix(~ e2 + lt[ic] + ll[ic] + factor(gcc[ic])) else model.matrix(~ e2 + lt[ic] + ll[ic]); f2 <- flp(Y2, X2, X2[, -2, drop = FALSE], seeds[b])
      nl[b] <- if (is.null(f2)) NA else sum(p.adjust(f2$p, "BH") < 0.05)
    }
    row$null_draws <- NNULL; row$null_calls_median <- median(nl, na.rm = TRUE); row$null_calls_mean <- mean(nl, na.rm = TRUE)
    row$null_calls_max <- max(nl, na.rm = TRUE); row$null_draws_with_any_call <- sum(nl > 0, na.rm = TRUE); row$null_calls_all <- paste(nl, collapse = ",")
    gate(paste0("F5_random_label_null_median_0_", SET), isTRUE(median(nl, na.rm = TRUE) == 0),
         paste0("draws with any call ", sum(nl > 0, na.rm = TRUE), "/", NNULL, "; max ", max(nl, na.rm = TRUE)))
  } else row$status <- "tested_null_not_run_comparison_le_enriched"
}
fin()
