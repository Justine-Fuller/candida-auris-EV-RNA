## s6f_word_aggregate.R -- Step 6F (v2), stage F5 r6 aggregation: all per-set gates, synthetic gates (planted word called;
## GC shift alone 0 calls), significant words, tag candidates (plan Q8 a: q < 0.05, same direction, >= 3 of 4 groups, per filter
## x region; untestable groups count as no call).
## r7: COV = GC or noGC. For noGC the GC-shift gate is INFO only (no GC control by design).
## Usage: Rscript s6f_word_aggregate.R <F4_sim dir> <F4_sets dir> <F5 per-set dir> <outdir> <GC|noGC>
a <- commandArgs(TRUE); FS <- a[1]; FR <- a[2]; D <- a[3]; OUT <- a[4]; COV <- a[5]
stopifnot(!dir.exists(OUT)); dir.create(OUT)
sets <- c(read.delim(file.path(FS, "sets_summary.tsv"))$set, read.delim(file.path(FR, "sets_summary.tsv"))$set)
miss <- sets[!file.exists(file.path(D, paste0(sets, ".summary.tsv")))]
G <- do.call(rbind, lapply(sets[!sets %in% miss], function(x) read.delim(file.path(D, paste0(x, ".gates.tsv")), colClasses = "character")))
G <- rbind(data.frame(gate = "F5_all_sets_have_output", result = ifelse(length(miss) == 0, "PASS", "FAIL"), detail = paste(miss, collapse = ";")), G)
S <- do.call(rbind, lapply(sets[!sets %in% miss], function(x) read.delim(file.path(D, paste0(x, ".summary.tsv")), colClasses = "character")))
write.table(S, file.path(OUT, "summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
pw <- readLines(file.path(FS, "sim_planted_word.txt"))[1]
pl <- S$set[S$variant == "planted"]; gs <- S$set[S$variant == "gcshift"]
wp <- read.delim(gzfile(file.path(D, paste0(pl, ".words.tsv.gz"))))
r <- wp[wp$word == pw, ]
G <- rbind(G, data.frame(gate = "F5_SIM_planted_word_called_higher_in_enriched",
  result = ifelse(r$q < 0.05 && r$direction == "higher_in_enriched", "PASS", "FAIL"),
  detail = paste0(pw, " p=", signif(r$p, 3), " q=", signif(r$q, 3), "; other calls ", sum(wp$q < 0.05) - (r$q < 0.05))))
ng <- as.integer(S$calls_higher_in_enriched[S$set == gs]) + as.integer(S$calls_higher_in_comparison[S$set == gs])
G <- rbind(G, data.frame(gate = if (COV == "GC") "F5_SIM_gcshift_alone_0_calls" else "INFO_noGC_SIM_gcshift_calls",
  result = if (COV == "GC") ifelse(ng == 0, "PASS", "FAIL") else "PASS", detail = paste0("calls ", ng)))
real <- S[S$variant == "real" & S$status != "not_testable_empty_side" & S$status != "not_testable_design", ]
SIG <- do.call(rbind, lapply(real$set, function(x) { w <- read.delim(gzfile(file.path(D, paste0(x, ".words.tsv.gz")))); w[w$q < 0.05, ] }))
write.table(SIG, file.path(OUT, "significant_words.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
CAND <- list()
if (!is.null(SIG) && nrow(SIG)) for (fl in unique(SIG$filter)) for (rg in unique(paste(SIG$rna_class, SIG$region))) {
  sub <- SIG[SIG$filter == fl & paste(SIG$rna_class, SIG$region) == rg, ]; if (!nrow(sub)) next
  tb <- aggregate(group ~ word + direction, sub, function(x) length(unique(x))); tb <- tb[tb$group >= 3, ]
  if (nrow(tb)) CAND[[paste(fl, rg)]] <- data.frame(filter = fl, rna_class_region = rg, word = tb$word, direction = tb$direction, groups_called = tb$group)
}
C <- if (length(CAND)) do.call(rbind, CAND) else data.frame(filter = character(), rna_class_region = character(), word = character(), direction = character(), groups_called = integer())
write.table(C, file.path(OUT, "tag_candidates_ge3of4.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(G, file.path(OUT, "gates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
print(G[G$result != "PASS" | grepl("SIM|null|all_sets", G$gate), ]); print(C)
print(S[, c("set", "enriched", "comparison", "calls_higher_in_enriched", "calls_higher_in_comparison", "min_q", "ttest_calls_r5_style",
            "null_calls_median", "null_calls_max", "null_draws_with_any_call", "words_reaching_BMAX", "status")])
ok <- all(G$result == "PASS"); cat(if (ok) "F5_ALL_GATES_PASS\n" else "F5_GATE_FAILURE\n"); quit(status = if (ok) 0 else 1)
