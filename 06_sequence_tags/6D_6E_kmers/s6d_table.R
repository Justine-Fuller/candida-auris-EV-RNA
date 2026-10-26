## s6d_table.R -- Step 6D (v2, 2026-10-12): annotation table of the FL r3 EV full-length PCGs, grouped like v1 6D.
## GO Biological Process terms as goslim_yeast words (v1 6D method; the v1 goslim_yeast.obo file is reused, data-version
## go/releases/2026-07-26). v1 used the QuickGO API (is_a + part_of ancestors, term itself included); the HPC has no internet,
## so ancestors come from GO.db (GOBPANCESTOR) -- the agreement with the v1 mapping is measured and reported (G6D_2/3, INFO).
## Plan gate: the group rule must reproduce the v1 groups of the 15 v1 genes (G6D_1).
## Groups (rule, so it can be re-run): "Histones" = eggNOG preferred name HHT/HHF/HTA/HTB/HTZ or V3 product "histone H2A/H2B/H3/H4";
## "Ribosomal proteins" = "ribosomal protein" in the eggNOG description or V3 product, or preferred name RPS*/RPL*/RPP*; else "Other".
## Within a group genes are ranked by mean B8441 caspofungin EV TPM. Words: all slim BP terms except the root; at most 8 shown
## (alphabetical, "+N more"), as v1.
## Usage: Rscript s6d_table.R <fl_genes.tsv> <goslim_yeast.obo> <v1 go_to_slim tsv> <v1 table tsv> <v1 eggnog_15_genes.tsv> <outdir>
suppressPackageStartupMessages({ library(GO.db); library(AnnotationDbi) })
a <- commandArgs(TRUE); FG <- a[1]; OBO <- a[2]; V1M <- a[3]; V1T <- a[4]; V1E <- a[5]; OUT <- a[6]
stopifnot(!dir.exists(OUT)); dir.create(OUT)
gates <- data.frame(gate = character(), result = character(), detail = character())
gate <- function(n, ok, d = "") gates[nrow(gates) + 1, ] <<- list(n, ifelse(isTRUE(ok), "PASS", "FAIL"), as.character(d))
info <- function(n, d) gates[nrow(gates) + 1, ] <<- list(n, "INFO", as.character(d))   # agreement measured, no threshold set
obo <- readLines(OBO); st <- grep("^\\[Term\\]", obo); en <- c(st[-1] - 1, length(obo))
terms <- do.call(rbind, lapply(seq_along(st), function(i) { b <- obo[st[i]:en[i]]
  data.frame(id = sub("^id: ", "", grep("^id: ", b, value = TRUE)[1]), name = sub("^name: ", "", grep("^name: ", b, value = TRUE)[1]),
             ns = sub("^namespace: ", "", grep("^namespace: ", b, value = TRUE)[1]), obs = any(b == "is_obsolete: true")) }))
slim <- terms[terms$ns == "biological_process" & !terms$obs & terms$id != "GO:0008150", ]
anc <- as.list(GOBPANCESTOR)
toslim <- function(g) { if (!(g %in% names(anc)) && !(g %in% slim$id)) return(character(0)); intersect(c(g, anc[[g]]), slim$id) }
bp_ids <- function(gl) { gl <- unique(gl[gl != ""]); ont <- suppressMessages(AnnotationDbi::select(GO.db, keys = gl, columns = "ONTOLOGY", keytype = "GOID"))
                         ont$GOID[!is.na(ont$ONTOLOGY) & ont$ONTOLOGY == "BP"] }
words <- function(gstr) { ids <- bp_ids(strsplit(gstr, "[,;]")[[1]]); s <- unique(unlist(lapply(ids, toslim)))
  sort(slim$name[match(s, slim$id)]) }
fmt <- function(w) if (!length(w)) "" else if (length(w) > 8) paste0(paste(w[1:8], collapse = "; "), " (+", length(w) - 8, " more)") else paste(w, collapse = "; ")
grp <- function(pref, desc, v3) ifelse(grepl("^(HHT|HHF|HTA|HTB|HTZ)[0-9]*$", pref, ignore.case = TRUE) | grepl("histone H(2A|2B|3|4)\\b", v3, ignore.case = TRUE), "Histones",
                             ifelse(grepl("ribosomal protein", paste(desc, v3), ignore.case = TRUE) | grepl("^RP[SLP]", pref, ignore.case = TRUE), "Ribosomal proteins", "Other"))
# ---- gate 1: group rule on the v1 15 genes reproduces the v1 groups
v1t <- read.delim(V1T, colClasses = "character"); v1e <- read.delim(V1E, colClasses = "character")
fg <- read.delim(FG, colClasses = "character")
m <- match(v1t$gene_id, fg$gene_id)
g1 <- grp(fg$eggNOG_Preferred_name[m], fg$eggNOG_Description[m], fg$V3_product[m])
gate("G6D_1_v1_15_genes_in_v2_list", all(!is.na(m)), sum(!is.na(m)))
gate("G6D_1_group_rule_reproduces_v1_groups", identical(g1, v1t$group), paste(sum(g1 == v1t$group), "/ 15"))
# ---- gate 2: GO.db slim mapping vs the v1 QuickGO mapping (per raw GO id seen in v1; BP ids only, root removed)
v1m <- read.delim(V1M, colClasses = "character"); v1m$slim_ids[is.na(v1m$slim_ids)] <- ""
bp <- bp_ids(v1m$go_id); v1m <- v1m[v1m$go_id %in% bp, ]
same <- vapply(seq_len(nrow(v1m)), function(i) { a1 <- setdiff(strsplit(v1m$slim_ids[i], ",")[[1]], "GO:0008150"); setequal(a1, toslim(v1m$go_id[i])) }, TRUE)
info("G6D_2_slim_mapping_agreement_with_v1_per_GO_id", paste0(sum(same), " / ", length(same), " BP ids identical"))
write.table(data.frame(go_id = v1m$go_id, v1 = v1m$slim_ids, v2 = vapply(v1m$go_id, function(g) paste(toslim(g), collapse = ","), ""), same = same),
            file.path(OUT, "slim_mapping_v1_vs_GOdb.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# ---- gate 3: slim words of the v1 15 genes from their v1 eggNOG GO ids vs the v1 table (total count + shown words)
w15 <- lapply(v1e$gene_id, function(g) words(v1e$GOs[v1e$gene_id == g])); names(w15) <- v1e$gene_id
cmp <- data.frame(gene_id = v1t$gene_id, v1_n = as.integer(v1t$n_go_bp_terms_total), v2_n = lengths(w15[v1t$gene_id]),
                  v1_words = v1t$go_bp_terms, v2_words = vapply(w15[v1t$gene_id], fmt, ""))
cmp$identical <- cmp$v1_words == cmp$v2_words
write.table(cmp, file.path(OUT, "v1_15_genes_words_v1_vs_v2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
info("G6D_3_v1_15_genes_slim_words_identical", paste(sum(cmp$identical), "/ 15 identical"))
# ---- the table
fg$group_v1_rule <- grp(fg$eggNOG_Preferred_name, fg$eggNOG_Description, fg$V3_product)
wl <- lapply(fg$eggNOG_GOs, words)
fg$n_go_bp_slim_terms <- lengths(wl); fg$go_bp_slim_terms <- vapply(wl, fmt, "")
fg$in_v1_15 <- fg$gene_id %in% v1t$gene_id
fg <- fg[order(match(fg$group_v1_rule, c("Histones", "Ribosomal proteins", "Other")), -as.numeric(fg$mean_TPM_B8441_caspofungin_EV)), ]
keep <- c("group_v1_rule", "gene_id", "full_length_groups", "eggNOG_Preferred_name", "V3_product", "eggNOG_Description",
          "Candida_albicans_orthologs", "Saccharomyces_cerevisiae_orthologs", "SignalP_prediction", "DeepTMHMM_has_TM", "PredGPI_prediction",
          "mean_TPM_B8441_caspofungin_EV", grep("^reads_", names(fg), value = TRUE), "n_go_bp_slim_terms", "go_bp_slim_terms", "in_v1_15")
write.table(fg[, keep], file.path(OUT, "fl_EV_PCG_annotation_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# slim-term counts across the B8441 caspofungin genes (descriptive)
bc <- grepl("B8441_caspofungin", fg$full_length_groups); tab <- sort(table(unlist(wl[bc])), decreasing = TRUE)
write.table(data.frame(slim_term = names(tab), genes = as.integer(tab), of = sum(bc)), file.path(OUT, "slim_term_counts_B8441_caspofungin.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(table(fg$group_v1_rule, bc)); write.table(gates, file.path(OUT, "gates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE); print(gates)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat(if (all(gates$result %in% c("PASS", "INFO"))) "S6D_TABLE_GATES_PASS\n" else "S6D_TABLE_GATE_FAILURE\n")
