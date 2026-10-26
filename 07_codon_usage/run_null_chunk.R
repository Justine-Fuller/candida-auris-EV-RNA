#!/usr/bin/env Rscript
# run_null_chunk.R  (v2 Step 6C (3) null, 2026-10-12)
# One SLURM array task: for its share of the 1,151 alternative codon sets, run the same ORA
# (GO BP/MF/CC + KEGG, frozen gene sets and universes from the real run) and record
#   - p for every tracked term (p = 1 when the term has no foreground gene, as clusterProfiler drops it)
#   - p for the planted term
#   - number of terms with BH p.adjust < 0.05 per ontology
# Args: <null_foregrounds.tsv> <real_run_dir> <outdir> <task_id> <n_tasks> <libdir>
args <- commandArgs(trailingOnly = TRUE)
fgf <- args[1]; rdir <- args[2]; odir <- args[3]; task <- as.integer(args[4]); ntask <- as.integer(args[5]); libdir <- args[6]
suppressMessages({ library(clusterProfiler); library(org.Hs.eg.db) })
source(file.path(libdir, "ora_lib.R"))
dir.create(odir, showWarnings = FALSE, recursive = TRUE)

gs_all <- readRDS(file.path(rdir, "gene_sets_frozen.rds"))
planted <- readRDS(file.path(rdir, "planted_term.rds"))
tracked <- read.delim(file.path(rdir, "tracked_terms.tsv"), quote = "", comment.char = "")
lines <- readLines(fgf)
mine <- lines[(seq_along(lines) - 1) %% ntask == (task - 1)]

out_p <- file.path(odir, sprintf("null_p_task%03d.tsv", task))
out_s <- file.path(odir, sprintf("null_summary_task%03d.tsv", task))
if (file.exists(out_p) || file.exists(out_s)) stop("output exists; refusing to overwrite")
pl <- list(); sl <- list()
for (ln in mine) {
  f <- strsplit(ln, "\t")[[1]]; sid <- f[1]; ens <- strsplit(f[2], ",")[[1]]
  st <- tryCatch({
    g_go <- fg_for_go(ens); g_kegg <- fg_for_kegg(ens)
    nsig <- c()
    for (ont in c("BP", "MF", "CC", "KEGG")) {
      res <- ora_terms(if (ont == "KEGG") g_kegg else g_go, gs_all[[ont]])
      nsig[ont] <- sum(res$p.adjust < 0.05)
      tr <- tracked$ID[tracked$ontology == ont]
      p <- res$pvalue[match(tr, res$ID)]; p[is.na(p)] <- 1
      pl[[length(pl) + 1]] <- data.frame(set_id = sid, ontology = ont, ID = tr, pvalue = p)
    }
    pr <- ora_terms(g_go, planted)
    pp <- if (nrow(pr) == 0) 1 else pr$pvalue
    pl[[length(pl) + 1]] <- data.frame(set_id = sid, ontology = "PLANTED", ID = "PLANTED_top50", pvalue = pp)
    data.frame(set_id = sid, status = "ok", n_fg = length(ens), n_fg_go = length(unique(g_go)),
               n_fg_kegg = length(unique(g_kegg)), nsig_BP = nsig["BP"], nsig_MF = nsig["MF"],
               nsig_CC = nsig["CC"], nsig_KEGG = nsig["KEGG"])
  }, error = function(e) data.frame(set_id = sid, status = paste("ERROR:", conditionMessage(e)), n_fg = length(ens),
                                    n_fg_go = NA, n_fg_kegg = NA, nsig_BP = NA, nsig_MF = NA, nsig_CC = NA, nsig_KEGG = NA))
  sl[[length(sl) + 1]] <- st
}
write.table(do.call(rbind, pl), out_p, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, sl), out_s, sep = "\t", quote = FALSE, row.names = FALSE)
cat("task", task, "sets", length(mine), "errors", sum(do.call(rbind, sl)$status != "ok"), "\n")
