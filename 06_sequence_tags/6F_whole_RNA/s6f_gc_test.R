## s6f_gc_test.R -- Step 6F (v2), stage F7 (2026-10-14; Justine 07:32 via the project assistant, after Caro's feedback that EV-enriched
## RNA may itself be GC-rich): is the region GC fraction different in EV-enriched vs comparison genes?
## Per set (filter x group x region; same gene sets as the word test): rank(region GC) ~ enriched + log10(Cell polyA TPM + 1)
##   + log10(transcript length); t of "enriched"; Freedman-Lane permutation p with Besag-Clifford sequential stopping (H, BMAX as F5).
##   BH over the 20 group x region tests of a filter; call = q < 0.05; direction = sign of the adjusted coefficient (median GC
##   reported alongside). Rule (plan Q8 a): same direction, q < 0.05, in >= 3 of 4 groups, per filter x region.
## Synthetic gate: the GC-shift set (F4_sim, enriched 3' UTRs A/T -> G/C at 10%) must be called higher in enriched.
## Random-label null (informational): 100 random enriched-size draws from the comparison genes per FX set; fraction with p < 0.05.
## Usage: Rscript s6f_gc_test.R <F4_sets dir> <F4_sim dir> <outdir> <seed> <n_null> <BMAX> <H>
a <- commandArgs(TRUE); FR <- a[1]; FS <- a[2]; OUT <- a[3]; SEED <- as.integer(a[4]); NNULL <- as.integer(a[5]); BMAX <- as.integer(a[6]); H <- as.integer(a[7])
stopifnot(!dir.exists(OUT)); dir.create(OUT)
gates <- data.frame(gate = character(), result = character(), detail = character())
gate <- function(n, ok, d = "") gates[nrow(gates) + 1, ] <<- list(n, ifelse(isTRUE(ok), "PASS", "FAIL"), as.character(d))
fl1 <- function(y, X, seed) {           # one response: observed t, coefficient, Freedman-Lane sequential permutation p
  q <- qr(X); q0 <- qr(X[, -2, drop = FALSE]); df <- nrow(X) - ncol(X); if (q$rank < ncol(X) || df < 1 || all(y == y[1])) return(NULL)
  vi <- chol2inv(qr.R(q))[2, 2]; ts <- function(v) qr.coef(q, v)[2] / sqrt(sum(qr.resid(q, v)^2) / df * vi)
  t0 <- ts(y); F0 <- qr.fitted(q0, y); R0 <- y - F0; ex <- 0L; L <- 0L; set.seed(seed)
  while (ex < H && L < BMAX) { L <- L + 1L; if (abs(ts(F0 + R0[sample.int(length(y))])) >= abs(t0) * (1 - 1e-9)) ex <- ex + 1L }
  list(t = unname(t0), coef = unname(qr.coef(q, y)[2]), p = if (ex >= H) H / L else (ex + 1) / (L + 1), L = L, df = df)
}
run <- function(dir, set, seed, nnull) {
  gn <- read.delim(file.path(dir, "sets", set, "genes.tsv"), colClasses = c("character", "character", "numeric", "numeric", "numeric", "numeric"))
  enr <- as.numeric(gn$set == "enriched"); lt <- log10(gn$cell_polyA_TPM + 1); ll <- log10(gn$len_whole); y <- rank(gn$region_GC)
  out <- data.frame(set = set, enriched = sum(enr), comparison = sum(1 - enr), median_GC_enriched = median(gn$region_GC[enr == 1]),
                    median_GC_comparison = median(gn$region_GC[enr == 0]), coef = NA, t = NA, p = NA, perms_used = NA, p_ttest = NA,
                    null_draws = 0, null_fraction_p_lt_0.05 = NA, status = "tested")
  if (sum(enr) < 1 || sum(1 - enr) < 1) { out$status <- "not_testable_empty_side"; return(out) }
  X <- model.matrix(~ enr + lt + ll); f <- fl1(y, X, seed)
  if (is.null(f)) { out$status <- "not_testable_design"; return(out) }
  m <- summary(lm(y ~ enr + lt + ll))$coefficients
  gate(paste0("F7_t_matches_lm_", set), abs(m["enr", 3] - f$t) < 1e-8 * max(1, abs(f$t)))
  out$coef <- f$coef; out$t <- f$t; out$p <- f$p; out$perms_used <- f$L; out$p_ttest <- m["enr", 4]
  if (nnull > 0 && sum(1 - enr) > sum(enr)) {
    ic <- which(enr == 0); set.seed(seed + 1e6); pn <- numeric(nnull)
    for (b in seq_len(nnull)) { e2 <- as.numeric(seq_along(ic) %in% sample(length(ic), sum(enr)))
      f2 <- fl1(rank(gn$region_GC[ic]), model.matrix(~ e2 + lt[ic] + ll[ic]), seed + b); pn[b] <- if (is.null(f2)) NA else f2$p }
    out$null_draws <- nnull; out$null_fraction_p_lt_0.05 <- mean(pn < 0.05, na.rm = TRUE)
  }
  out
}
ss <- read.delim(file.path(FR, "sets_summary.tsv"), colClasses = "character")
R <- do.call(rbind, lapply(seq_len(nrow(ss)), function(i) cbind(ss[i, c("filter", "group", "rna_class", "region")],
       run(FR, ss$set[i], SEED + i, if (ss$filter[i] == "FX") NNULL else 0))))
R$q <- NA; for (fl in unique(R$filter)) { k <- R$filter == fl & !is.na(R$p); R$q[k] <- p.adjust(R$p[k], "BH") }
R$direction <- ifelse(is.na(R$coef), NA, ifelse(R$coef > 0, "GC_higher_in_enriched", "GC_higher_in_comparison"))
R$call <- !is.na(R$q) & R$q < 0.05
write.table(R, file.path(OUT, "gc_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
C <- aggregate(call ~ filter + rna_class + region + direction, R[!is.na(R$direction), ], sum)
C$ge3of4 <- C$call >= 3; write.table(C, file.path(OUT, "gc_rule_ge3of4.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
sim <- read.delim(file.path(FS, "sets_summary.tsv"), colClasses = "character"); gs <- sim$set[sim$variant == "gcshift"]
rs <- run(FS, gs, SEED, 0); rr <- run(FS, sim$set[sim$variant == "planted"], SEED, 0)
gate("F7_SIM_gcshift_called_GC_higher_in_enriched", !is.na(rs$p) && rs$p < 0.05 && rs$coef > 0,
     paste0("p=", signif(rs$p, 3), " median GC enriched ", round(rs$median_GC_enriched, 3), " vs ", round(rs$median_GC_comparison, 3),
            "; unshifted (planted-word) set p=", signif(rr$p, 3)))
write.table(gates, file.path(OUT, "gates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
print(gates[gates$result != "PASS" | grepl("SIM", gates$gate), ]); print(R[, c("filter", "group", "region", "rna_class", "enriched", "comparison",
  "median_GC_enriched", "median_GC_comparison", "p", "q", "direction", "null_fraction_p_lt_0.05")]); print(C[C$ge3of4, ])
ok <- all(gates$result == "PASS"); cat(if (ok) "F7_ALL_GATES_PASS\n" else "F7_GATE_FAILURE\n"); quit(status = if (ok) 0 else 1)
