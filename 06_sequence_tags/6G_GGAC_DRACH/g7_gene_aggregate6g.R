## g7_gene_aggregate6g.R -- Step 6G (v2, 2026-10-14), gene level aggregation (after g6_gene_perm6g.R, one COV at a time).
## Gates: every set has output and per-set gates; synthetic planted GGACT -> GGAC AND DRACH called higher in enriched; GC shift alone
## -> 0 calls (GC version; INFO for noGC); random-label null median 0 is a per-set gate (from g6).
## Results: GGAC and DRACH rows per set (q < 0.05 over 257 units; direction = adjusted coefficient sign), GGAC rank among the 256
## 4-mers by p; >= 3 of 4 groups rule (plan Q8 a / 6G fixed method) per region for GGAC and DRACH.
## Usage: Rscript g7_gene_aggregate6g.R <F4_sim dir> <F4_sets dir> <per-set dir> <outdir> <GC|noGC>
a <- commandArgs(TRUE); FS <- a[1]; FR <- a[2]; D <- a[3]; OUT <- a[4]; COV <- a[5]
stopifnot(!dir.exists(OUT)); dir.create(OUT)
sets <- c(read.delim(file.path(FS, "sets_summary.tsv"))$set, read.delim(file.path(FR, "sets_summary.tsv"))$set)
miss <- sets[!file.exists(file.path(D, paste0(sets, ".summary.tsv")))]
G <- do.call(rbind, lapply(sets[!sets %in% miss], function(x) read.delim(file.path(D, paste0(x, ".gates.tsv")), colClasses = "character")))
G <- rbind(data.frame(gate = "G6_all_sets_have_output", result = ifelse(length(miss) == 0, "PASS", "FAIL"), detail = paste(miss, collapse = ";")), G)
S <- do.call(rbind, lapply(sets[!sets %in% miss], function(x) read.delim(file.path(D, paste0(x, ".summary.tsv")), colClasses = "character")))
write.table(S, file.path(OUT, "summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
rw <- function(x) read.delim(gzfile(file.path(D, paste0(x, ".words.tsv.gz"))))
pl <- S$set[S$variant == "planted"]; gs <- S$set[S$variant == "gcshift"]
wp <- rw(pl); r <- wp[wp$word %in% c("GGAC", "DRACH"), ]
G <- rbind(G, data.frame(gate = "G6_SIM_planted_GGACT_GGAC_and_DRACH_called_higher_in_enriched",
  result = ifelse(nrow(r) == 2 && all(r$q < 0.05 & r$direction == "higher_in_enriched"), "PASS", "FAIL"),
  detail = paste(paste0(r$word, " q=", signif(r$q, 3), " ", r$direction), collapse = "; ")))
ng <- as.integer(S$calls_higher_in_enriched[S$set == gs]) + as.integer(S$calls_higher_in_comparison[S$set == gs])
G <- rbind(G, data.frame(gate = if (COV == "GC") "G6_SIM_gcshift_alone_0_calls" else "INFO_noGC_SIM_gcshift_calls",
  result = if (COV == "GC") ifelse(ng == 0, "PASS", "FAIL") else "PASS", detail = paste0("calls ", ng)))
real <- S[S$variant == "real" & grepl("^tested", S$status), ]
M <- do.call(rbind, lapply(real$set, function(x) { w <- rw(x); k4 <- w[w$word != "DRACH", ]
  rk <- rank(k4$p, ties.method = "min"); g <- w[w$word %in% c("GGAC", "DRACH"), ]
  g$GGAC_rank_among_256_by_p <- ifelse(g$word == "GGAC", rk[k4$word == "GGAC"], NA)
  g$calls_in_set <- sum(w$q < 0.05); g }))
keep <- c("set", "covariates", "group", "rna_class", "region", "word", "mean_rate_per_kb_enriched", "mean_rate_per_kb_comparison",
          "genes_with_word_enriched", "genes_with_word_comparison", "coef_enriched_rank", "p", "q", "direction", "p_ttest",
          "GGAC_rank_among_256_by_p", "calls_in_set")
M <- M[, keep]; write.table(M, file.path(OUT, "GGAC_DRACH_gene_level.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
R3 <- aggregate(cbind(EV_enriched_higher = q < 0.05 & direction == "higher_in_enriched",
                      comparison_higher = q < 0.05 & direction == "higher_in_comparison") ~ word + rna_class + region, M, sum)
R3$ge3of4 <- ifelse(R3$EV_enriched_higher >= 3, "higher_in_EV_enriched", ifelse(R3$comparison_higher >= 3, "higher_in_comparison", "none"))
write.table(R3, file.path(OUT, "GGAC_DRACH_rule_ge3of4.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(G, file.path(OUT, "gates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
print(G[G$result != "PASS" | grepl("SIM|null|all_sets", G$gate), ]); print(R3); print(M[, c("set", "word", "q", "direction", "GGAC_rank_among_256_by_p")])
ok <- all(G$result == "PASS"); cat(if (ok) "G6_ALL_GATES_PASS\n" else "G6_GATE_FAILURE\n"); quit(status = if (ok) 0 else 1)
