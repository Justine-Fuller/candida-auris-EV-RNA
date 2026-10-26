## kmer_edger6e.R -- Step 6E (v2, 2026-10-12), stage K3: EV vs Cell sRNA test of each 6-mer's SHARE, per group and RNA type.
## Approved plan (00_PLAN_DRAFT_6D_6E_v2_2026-10-12.md): edgeR quasi-likelihood test, 3 EV vs 3 Cell sRNA (unpaired),
## BH over the 4,096 words per type per group, q < 0.05, direction from observed shares.
##  primary   (Q4 approach 1): Cell sRNA = GC-matched draw 1 (K2 "__m"), EV = reads in matchable GC bins ("__m");
##            "library size" = all 6-mers of that type in that library -> tests the word's share.
##  stability: the same with Cell draws 2 and 3 (information only).
##  crosscheck (Q4 approach 2): all reads, no matching; each word's offset = log(total 6-mers of the word's GC class
##            (G+C count 0-6) in that library) -> tests the word's share WITHIN its GC class (GC class as covariate).
##  naive     : all reads, no GC control (information only: shows what GC control removes).
##  null      : primary data, the 9 label swaps that exchange one EV and one Cell library (2 EV + 1 Cell vs 1 EV + 2 Cell).
## Words with 0 counts in all tested libraries are not fitted (p = 1, kept in the BH over 4,096). Libraries with 0 6-mers of
## the type are dropped; a test needs >= 2 libraries per side. edgeR defaults otherwise (glmQLFit legacy = FALSE, robust = FALSE;
## no TMM: the library size IS the denominator of the share).
## Usage: Rscript kmer_edger6e.R <manifest.tsv> <K2 dir> <outdir> [types=PCG_CDS,lncRNA,tRNA,5S_rRNA]
suppressMessages({ library(limma); library(edgeR) })
a <- commandArgs(TRUE); MAN <- a[1]; K2 <- a[2]; OUT <- a[3]
TYPES <- if (length(a) > 3) strsplit(a[4], ",")[[1]] else c("PCG_CDS", "lncRNA", "tRNA", "5S_rRNA")
stopifnot(!dir.exists(OUT)); dir.create(OUT, recursive = TRUE)
man <- read.delim(MAN, colClasses = "character")
gates <- data.frame(gate = character(), result = character(), detail = character())
gate <- function(n, ok, d = "") gates[nrow(gates) + 1, ] <<- list(n, ifelse(isTRUE(ok), "PASS", "FAIL"), as.character(d))
qltest <- function(cnt, isEV, libsize = NULL, offset = NULL) {
  # cnt: 4096 x libs. Returns p, logFC (edgeR coef, log2) and observed mean shares, direction.
  tot <- if (is.null(offset)) libsize else NULL
  okl <- if (is.null(offset)) libsize > 0 else colSums(cnt) > 0
  n <- nrow(cnt); p <- rep(1, n); lfc <- rep(NA_real_, n); tested <- rep(FALSE, n)
  den <- if (is.null(offset)) matrix(libsize, n, ncol(cnt), byrow = TRUE) else exp(offset)
  sh <- ifelse(den > 0, cnt / den, NA)
  mE <- rowMeans(sh[, isEV & okl, drop = FALSE], na.rm = TRUE); mC <- rowMeans(sh[, !isEV & okl, drop = FALSE], na.rm = TRUE)
  if (sum(isEV & okl) >= 2 && sum(!isEV & okl) >= 2) {
    c2 <- cnt[, okl, drop = FALSE]; ev2 <- isEV[okl]; keep <- rowSums(c2) > 0
    grp <- factor(ifelse(ev2, "EV", "Cell"), levels = c("Cell", "EV")); design <- model.matrix(~ grp)
    y <- DGEList(counts = c2[keep, , drop = FALSE], lib.size = if (is.null(offset)) libsize[okl] else colSums(c2))
    if (!is.null(offset)) { o <- offset[keep, okl, drop = FALSE]; o[!is.finite(o)] <- 0; y$offset <- o }
    fit <- glmQLFit(y, design); tt <- glmQLFTest(fit, coef = 2)
    p[keep] <- tt$table$PValue; lfc[keep] <- tt$table$logFC; tested[keep] <- TRUE
  } else p[] <- NA
  q <- p.adjust(p, "BH")
  eq <- !is.na(mE) & !is.na(mC) & abs(mE - mC) <= 1e-12 * pmax(mE, mC)
  dirn <- ifelse(is.na(mE) | is.na(mC), NA, ifelse(eq, "equal", ifelse(mE > mC, "higher_in_EV", "higher_in_Cell")))
  list(p = p, q = q, lfc = lfc, mE = mE, mC = mC, dir = dirn, tested = tested, nEV = sum(isEV & okl), nCell = sum(!isEV & okl))
}
calls <- function(r, d) sum(!is.na(r$q) & r$q < 0.05 & r$dir == d, na.rm = TRUE)
SUM <- list(); NUL <- list()
for (G in unique(man$group)) for (T in TYPES) {
  f <- file.path(K2, paste0(G, ".", T, ".counts.tsv.gz")); X <- read.delim(gzfile(f), colClasses = c(word = "character"))
  ev <- man$library[man$group == G & man$compartment == "EV"]; ce <- man$library[man$group == G & man$compartment == "Cell_sRNA"]
  libs <- c(ce, ev); isEV <- libs %in% ev
  gate(paste0("K3_4096_words_", G, "_", T), nrow(X) == 4096 && !anyDuplicated(X$word), nrow(X))
  M <- function(suf) { cn <- paste0(libs, suf); cn <- ifelse(isEV & suf %in% c("__m2", "__m3"), paste0(libs, "__m"), cn)
                       if (!all(cn %in% colnames(X))) return(NULL); m <- as.matrix(X[, cn]); colnames(m) <- libs; m }
  out <- data.frame(group = G, type = T, word = X$word, gc_class = X$gc_class)
  res <- list()
  Pm <- M("__m")
  if (!is.null(Pm)) {
    res$primary <- qltest(Pm, isEV, libsize = colSums(Pm))
    for (d in c("__m2", "__m3")) { Md <- M(d); if (!is.null(Md)) res[[paste0("draw", substr(d, 4, 4))]] <- qltest(Md, isEV, libsize = colSums(Md)) }
  }
  Pa <- M("__all")
  ctot <- apply(Pa, 2, function(v) tapply(v, X$gc_class, sum)); off <- log(ctot[as.character(X$gc_class), , drop = FALSE])
  res$gcclass <- qltest(Pa, isEV, offset = off)
  res$naive <- qltest(Pa, isEV, libsize = colSums(Pa))
  for (k in names(res)) { r <- res[[k]]
    out[[paste0("p_", k)]] <- r$p; out[[paste0("q_", k)]] <- r$q; out[[paste0("log2FC_", k)]] <- r$lfc
    out[[paste0("share_EV_", k)]] <- r$mE; out[[paste0("share_Cell_", k)]] <- r$mC; out[[paste0("direction_", k)]] <- r$dir
    tq <- r$q[!is.na(r$p)]; if (length(tq)) { pp <- r$p[!is.na(r$p)]
      gate(paste0("K3_BH_over_4096_", k, "_", G, "_", T), length(pp) == 4096 && isTRUE(all.equal(p.adjust(pp, "BH"), tq)) && all(pp >= 0 & pp <= 1), sum(r$tested)) }
  }
  for (j in seq_along(libs)) out[[paste0("count_primary_", libs[j])]] <- if (!is.null(Pm)) Pm[, j] else NA
  gz <- gzfile(file.path(OUT, paste0(G, ".", T, ".words.tsv.gz")), "w"); write.table(out, gz, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"); close(gz)
  s <- data.frame(group = G, type = T, kmers_primary_EV = if (!is.null(Pm)) sum(Pm[, isEV]) else NA,
                  kmers_primary_Cell = if (!is.null(Pm)) sum(Pm[, !isEV]) else NA)
  for (k in names(res)) { r <- res[[k]]
    s[[paste0(k, "_libs_EV_Cell")]] <- paste0(r$nEV, "v", r$nCell); s[[paste0(k, "_words_fitted")]] <- sum(r$tested)
    s[[paste0(k, "_q05_higher_in_EV")]] <- calls(r, "higher_in_EV"); s[[paste0(k, "_q05_higher_in_Cell")]] <- calls(r, "higher_in_Cell") }
  SUM[[length(SUM) + 1]] <- s
  # ---- label-swap null on the primary data
  if (!is.null(Pm)) {
    for (i in seq_along(ev)) for (j in seq_along(ce)) {
      isEVs <- isEV; isEVs[libs == ev[i]] <- FALSE; isEVs[libs == ce[j]] <- TRUE
      r <- qltest(Pm, isEVs, libsize = colSums(Pm))
      NUL[[length(NUL) + 1]] <- data.frame(group = G, type = T, swapped = paste(ev[i], ce[j], sep = "<->"),
        calls_q05 = sum(!is.na(r$q) & r$q < 0.05), calls_primary_q05 = sum(!is.na(res$primary$q) & res$primary$q < 0.05),
        min_q = suppressWarnings(min(r$q, na.rm = TRUE)))
    }
  }
  cat(G, T, "done\n")
}
S <- do.call(rbind, SUM); write.table(S, file.path(OUT, "K3_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (length(NUL)) write.table(do.call(rbind, NUL), file.path(OUT, "K3_label_swap_null.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gates, file.path(OUT, "gates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
print(S); ok <- all(gates$result == "PASS"); cat(if (ok) "K3_ALL_GATES_PASS\n" else "K3_GATE_FAILURE\n"); quit(status = if (ok) 0 else 1)
