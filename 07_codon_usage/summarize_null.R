#!/usr/bin/env Rscript
# summarize_null.R  (v2 Step 6C (3) null, 2026-10-12)
# Empirical fraction per tracked term = share of the 1,151 alternative codon sets in which the real
# p is strictly smaller than the alternative-set p (p = 1 when the term has no foreground gene).
# Specific = fraction >= 0.95 (approved spec). Also: completeness gate and planted-term gate.
# Args: <real_run_dir> <null_dir> <null_sets.tsv> <outdir>
args <- commandArgs(trailingOnly = TRUE)
rdir <- args[1]; ndir <- args[2]; setsf <- args[3]; odir <- args[4]
n_expected <- if (length(args) >= 5) as.integer(args[5]) else 1151L   # r3: number of alternative sets
dir.create(odir, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) read.delim(f, quote = "", comment.char = "", stringsAsFactors = FALSE)

sets <- rd(setsf)
summ <- do.call(rbind, lapply(list.files(ndir, "^null_summary_task\\d+\\.tsv$", full.names = TRUE), rd))
pv <- do.call(rbind, lapply(list.files(ndir, "^null_p_task\\d+\\.tsv$", full.names = TRUE), rd))
tracked <- rd(file.path(rdir, "tracked_terms.tsv"))

n_sets <- length(unique(summ$set_id)); n_err <- sum(summ$status != "ok")
exp_rows <- n_sets * (nrow(tracked) + 1)
complete <- n_sets == n_expected && setequal(summ$set_id, sets$set_id) && n_err == 0 && !any(duplicated(summ$set_id)) &&
  nrow(pv) == exp_rows && !any(is.na(pv$pvalue))

real <- do.call(rbind, lapply(c("BP", "MF", "CC", "KEGG"), function(ont) {
  tag <- if (ont == "KEGG") "KEGG" else paste0("GO_", ont)
  d <- rd(file.path(rdir, paste0(tag, "_all.tsv")))
  d$ontology <- ont
  d[, c("ontology", "ID", "Description", "GeneRatio", "BgRatio", "FoldEnrichment", "pvalue", "p.adjust", "Count", "geneID")]
}))
tr <- merge(tracked, real, by = c("ontology", "ID"), all.x = TRUE)
tr$pvalue[is.na(tr$pvalue)] <- 1   # term had no foreground gene in the real set

key <- paste(pv$ontology, pv$ID)
byk <- split(pv$pvalue, key)   # r3: one split instead of a full scan per term (same values; 30x more sets)
tr$n_alt_sets <- NA; tr$n_alt_p_greater <- NA; tr$n_alt_p_le_real <- NA
for (i in seq_len(nrow(tr))) {
  alt <- byk[[paste(tr$ontology[i], tr$ID[i])]]; if (is.null(alt)) alt <- numeric(0)
  tr$n_alt_sets[i] <- length(alt)
  tr$n_alt_p_greater[i] <- sum(tr$pvalue[i] < alt)
  tr$n_alt_p_le_real[i] <- sum(alt <= tr$pvalue[i])
}
tr$empirical_fraction <- tr$n_alt_p_greater / tr$n_alt_sets
tr$specific <- ifelse(tr$empirical_fraction >= 0.95, "yes", "no")
tr <- tr[order(tr$ontology, tr$pvalue), ]
write.table(tr, file.path(odir, "specificity_tracked_terms.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
sig <- tr[tr$real_sig, ]
write.table(sig, file.path(odir, "specificity_significant_terms.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# planted
plr <- readLines(file.path(rdir, "planted_term_real.txt"))
p_real_pl <- as.numeric(sub(".*real p = ", "", plr[2]))
alt_pl <- pv$pvalue[pv$ontology == "PLANTED"]
frac_pl <- mean(p_real_pl < alt_pl)
planted_pass <- frac_pl >= 0.95

# number of significant terms per ontology: real vs alternative sets
nsig_real <- sapply(c("BP", "MF", "CC", "KEGG"), function(o) sum(real$ontology == o & real$p.adjust < 0.05))
nsig <- data.frame(ontology = names(nsig_real), real = as.integer(nsig_real),
  alt_median = sapply(names(nsig_real), function(o) median(summ[[paste0("nsig_", o)]])),
  alt_min = sapply(names(nsig_real), function(o) min(summ[[paste0("nsig_", o)]])),
  alt_max = sapply(names(nsig_real), function(o) max(summ[[paste0("nsig_", o)]])),
  alt_sets_with_ge_real = sapply(names(nsig_real), function(o) sum(summ[[paste0("nsig_", o)]] >= nsig_real[o])))
write.table(nsig, file.path(odir, "n_significant_real_vs_null.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

con <- file(file.path(odir, "null_gates_and_summary.txt"), "w")
w <- function(...) writeLines(sprintf(...), con)
w("Null completeness: sets with output %d / %d; errors %d; p rows %d (expected %d) -> %s",
  n_sets, n_expected, n_err, nrow(pv), exp_rows, if (complete) "PASS" else "FAIL")
w("Planted term (top 50 real-foreground genes, BP universe): real p = %.4g; empirical fraction = %.4f -> %s",
  p_real_pl, frac_pl, if (planted_pass) "PASS (specific)" else "FAIL")
for (o in c("BP", "MF", "CC", "KEGG")) {
  s <- sig[sig$ontology == o, ]
  w("%s: %d significant (BH < 0.05); %d specific (fraction >= 0.95)", o, nrow(s), sum(s$specific == "yes"))
}
w("Total: %d significant; %d specific", nrow(sig), sum(sig$specific == "yes"))
close(con)
cat(readLines(file.path(odir, "null_gates_and_summary.txt")), sep = "\n")
print(nsig)
