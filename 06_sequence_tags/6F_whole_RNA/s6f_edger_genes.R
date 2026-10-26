## s6f_edger_genes.R -- Step 6F (v2, 2026-10-13), stage F3: gene-level EV vs Cell sRNA test per group (plan Q1 a), and the
## EV-enriched / comparison-pool gene counts per region subset, reported BEFORE any word test (plan Q8 note).
## Plan rules: EV-enriched = q < 0.05 and log2(EV/Cell share) > 1; comparison pool = q > 0.5 and |log2(EV/Cell share)| < 0.5
## (Q2 a; matching to the enriched set is a later stage). Shares = observed (read count / library total over all
## non-competitor gene models), mean of the 3 replicates per side; DIRECTION AND FOLD FROM OBSERVED SHARES (project rule).
## Not fixed by the plan, so BOTH are run and reported (Justine chooses):
##   filter FE = step-E read rule (>= 10 reads in >= 2/3 EV AND >= 2/3 Cell sRNA of the group);
##   filter FX = edgeR filterByExpr(y, group) defaults.
## Implementation: edgeR 4.8.2; DGEList(lib.size = column sums over all non-competitor models), TMM (normLibSizes),
##   design ~ group (Cell reference), glmQLFit defaults, glmQLFTest(coef = "groupEV"); BH within group over all tested genes
##   (PCG + lncRNA together). Label-swap null (informational): the 9 swaps exchanging one EV and one Cell library.
## Usage: Rscript s6f_edger_genes.R <manifest_6E.tsv> <counts_dir> <gene_sequences.tsv.gz> <tpm_with_rRNA_per_library.tsv.gz>
##        <manifest_all36.tsv> <outdir>
suppressMessages({ library(limma); library(edgeR) })
a <- commandArgs(TRUE); MAN <- a[1]; CD <- a[2]; SEQ <- a[3]; TPM <- a[4]; M36 <- a[5]; OUT <- a[6]
stopifnot(!dir.exists(OUT)); dir.create(OUT, recursive = TRUE)
Q_EN <- 0.05; LFC_EN <- 1; Q_CMP <- 0.5; LFC_CMP <- 0.5
gates <- data.frame(gate = character(), result = character(), detail = character())
gate <- function(n, ok, d = "") gates[nrow(gates) + 1, ] <<- list(n, ifelse(isTRUE(ok), "PASS", "FAIL"), as.character(d))
man <- read.delim(MAN, colClasses = "character"); m36 <- read.delim(M36, colClasses = "character")
sq <- read.delim(gzfile(SEQ), colClasses = "character")[, c("model_id", "gene_id", "rna_class", "role", "len_UTR5", "len_CDS", "len_UTR3", "len_whole")]
for (k in c("len_UTR5", "len_CDS", "len_UTR3", "len_whole")) sq[[k]] <- as.integer(sq[[k]])
tpm <- read.delim(gzfile(TPM), colClasses = "character"); tpm <- tpm[tpm$compartment == "Cell_polyA", ]
rdc <- function(L) read.delim(file.path(CD, paste0(L, ".gene_counts.tsv")), colClasses = c("character", "character", "numeric"))
obs <- function(cnt, le, lc) {
  sh <- sweep(cnt, 2, colSums(cnt), "/"); mE <- rowMeans(sh[, le, drop = FALSE]); mC <- rowMeans(sh[, lc, drop = FALSE])
  eq <- abs(mE - mC) <= 1e-12 * pmax(mE, mC)
  l2 <- ifelse(eq, 0, log2(mE / mC))   # +Inf when Cell share 0, -Inf when EV share 0, NaN when both 0
  list(mE = mE, mC = mC, l2 = l2)
}
fitq <- function(cnt, keep, lib, grp) {
  design <- model.matrix(~ grp); colnames(design) <- c("(Intercept)", "groupEV")
  y <- DGEList(counts = cnt[keep, , drop = FALSE], lib.size = lib); y <- normLibSizes(y)
  fit <- glmQLFit(y, design); t <- glmQLFTest(fit, coef = "groupEV")
  p <- t$table$PValue; list(p = p, q = p.adjust(p, "BH"), logFC = t$table$logFC, nf = y$samples$norm.factors)
}
SUM <- list(); NUL <- list()
for (g in unique(man$group)) {
  lc <- man$library[man$group == g & man$compartment != "EV"]; le <- man$library[man$group == g & man$compartment == "EV"]
  stopifnot(length(le) == 3, length(lc) == 3); libs <- c(lc, le)
  Tb <- lapply(libs, rdc)
  gate(paste0("F3_same_models_all_6_", g), all(vapply(Tb, function(t) identical(t$model_id, Tb[[1]]$model_id), TRUE)) &&
       setequal(Tb[[1]]$model_id, sq$model_id), nrow(Tb[[1]]))
  ids <- Tb[[1]]$model_id; cls <- Tb[[1]]$rna_class
  cnt <- sapply(Tb, function(t) t$reads); colnames(cnt) <- libs; rownames(cnt) <- ids
  comp <- sq$role[match(ids, sq$model_id)] == "competitor_only"
  gate(paste0("F3_competitor_models_have_0_reads_", g), all(cnt[comp, ] == 0), sum(comp))
  cnt <- cnt[!comp, ]; ids <- ids[!comp]; cls <- cls[!comp]
  lib <- colSums(cnt); grp <- factor(c(rep("Cell", 3), rep("EV", 3)), levels = c("Cell", "EV"))
  o <- obs(cnt, le, lc)
  # Cell polyA TPM (mean of the group's 3 polyA libraries; for later abundance matching)
  lp <- m36$library[m36$group == g & m36$kind == "polyA"]; stopifnot(length(lp) == 3)
  tp <- tpm[tpm$library %in% lp, ]; tpw <- tapply(as.numeric(tp$TPM_with_rRNA), tp$unit, mean)
  ptpm <- as.numeric(tpw[ids])
  gate(paste0("F3_polyA_TPM_found_for_all_models_", g), !anyNA(ptpm), sum(is.na(ptpm)))
  ev_ok <- rowSums(cnt[, le] >= 10) >= 2; ce_ok <- rowSums(cnt[, lc] >= 10) >= 2
  y0 <- DGEList(counts = cnt, lib.size = lib, group = grp)
  FILT <- list(FE = ev_ok & ce_ok, FX = filterByExpr(y0, group = grp))
  for (fn in names(FILT)) {
    keep <- FILT[[fn]]; r <- fitq(cnt, keep, lib, grp)
    p <- rep(NA_real_, length(ids)); q <- p; lfc <- p; p[keep] <- r$p; q[keep] <- r$q; lfc[keep] <- r$logFC
    gate(paste0("F3_p_in_0_1_", g, "_", fn), all(r$p >= 0 & r$p <= 1) && !anyNA(r$p), sum(keep))
    oo <- order(r$p); mm <- length(r$p); qq <- rev(cummin(rev(pmin(1, r$p[oo] * mm / seq_len(mm))))); q2 <- numeric(mm); q2[oo] <- qq
    gate(paste0("F3_BH_recomputed_", g, "_", fn), isTRUE(all.equal(q2, r$q, tolerance = 1e-12)), mm)
    dn <- sign(o$l2[keep]); ok <- is.finite(dn) & dn != 0
    agree <- mean(sign(r$logFC)[ok] == dn[ok])
    enr <- keep & !is.na(q) & q < Q_EN & !is.na(o$l2) & o$l2 > LFC_EN
    cmp <- keep & !is.na(q) & q > Q_CMP & !is.na(o$l2) & abs(o$l2) < LFC_CMP
    tab <- data.frame(group = g, filter = fn, model_id = ids, rna_class = cls, tested = keep, cnt,
                      mean_share_EV = o$mE, mean_share_Cell = o$mC, log2_share_EV_over_Cell = o$l2, edgeR_logFC = lfc,
                      p = p, q = q, ev_enriched = enr, comparison_pool = cmp, cell_polyA_TPM_mean = ptpm, check.names = FALSE)
    tab <- cbind(tab, sq[match(ids, sq$model_id), c("len_UTR5", "len_CDS", "len_UTR3", "len_whole")])
    gz <- gzfile(file.path(OUT, paste0(g, ".", fn, ".gene_tests.tsv.gz")), "w")
    write.table(tab, gz, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"); close(gz)
    P <- cls == "PCG_CDS"; L <- cls == "lncRNA"; u5 <- tab$len_UTR5 > 0; u3 <- tab$len_UTR3 > 0
    n <- function(s) c(enr = sum(enr & s), cmp = sum(cmp & s), tested = sum(keep & s))
    sub <- list(PCG_any = P, PCG_with_UTR5 = P & u5, PCG_with_UTR3 = P & u3, PCG_with_both_UTRs = P & u5 & u3, lncRNA = L)
    s <- data.frame(group = g, filter = fn, genes_tested = sum(keep), PCG_tested = sum(keep & P), lncRNA_tested = sum(keep & L),
                    q005_any_direction = sum(keep & q < 0.05, na.rm = TRUE),
                    q005_EV_higher_any_fold = sum(keep & q < 0.05 & o$l2 > 0, na.rm = TRUE),
                    edgeR_sign_agrees_observed = round(agree, 4),
                    TMM_factors = paste(round(r$nf, 3), collapse = ","))
    for (k in names(sub)) { v <- n(sub[[k]]); s[[paste0(k, "_EV_enriched")]] <- v[["enr"]]; s[[paste0(k, "_comparison_pool")]] <- v[["cmp"]] }
    s$EV_enriched_with_polyA_TPM_0 <- sum(enr & ptpm == 0)
    # label-swap null (informational): enriched calls when one EV and one Cell library swap labels
    nl <- c()
    for (i in 1:3) for (j in 1:3) {
      ce2 <- lc; ev2 <- le; ce2[i] <- le[j]; ev2[j] <- lc[i]; c2 <- cnt[, c(ce2, ev2)]
      o2 <- obs(c2, ev2, ce2); k2 <- if (fn == "FE") rowSums(c2[, ev2] >= 10) >= 2 & rowSums(c2[, ce2] >= 10) >= 2 else
        filterByExpr(DGEList(counts = c2, lib.size = colSums(c2), group = grp), group = grp)
      r2 <- fitq(c2, k2, colSums(c2), grp); q2b <- rep(NA_real_, length(ids)); q2b[k2] <- r2$q
      nl <- c(nl, sum(k2 & q2b < Q_EN & o2$l2 > LFC_EN, na.rm = TRUE))
    }
    s$label_swap_null_EV_enriched_calls_9_swaps <- paste(nl, collapse = ",")
    SUM[[paste(g, fn)]] <- s
  }
}
S <- do.call(rbind, SUM); write.table(S, file.path(OUT, "summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gates, file.path(OUT, "gates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
print(gates); print(t(S))
ok <- all(gates$result == "PASS"); cat(if (ok) "F3_ALL_GATES_PASS\n" else "F3_GATE_FAILURE\n"); quit(status = if (ok) 0 else 1)
