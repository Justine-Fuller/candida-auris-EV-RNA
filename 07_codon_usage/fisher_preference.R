#!/usr/bin/env Rscript
# fisher_preference.R  (v2 Step 6C (2), 2026-10-12)
# Two-sided Fisher exact test: do the 6 chosen families decode the preferred codon more/less often
# than the other 37 families? Families with no defined value (CUG under human code) are left out.
# Args: <fisher_counts.tsv> <out.tsv>
a <- commandArgs(trailingOnly = TRUE)
d <- read.delim(a[1])
out <- do.call(rbind, lapply(c("human", "cauris"), function(sp) {
  x <- d[d$species == sp, ]; c6 <- x[grepl("^chosen", x$group), ]; o <- x[grepl("^other", x$group), ]; al <- x[x$group == "all43", ]  # r3: any chosen/other size
  m <- matrix(c(c6$n_decode_preferred, c6$n_not, o$n_decode_preferred, o$n_not), 2, byrow = TRUE)
  ft <- fisher.test(m)
  data.frame(species = sp, chosen6_pref = c6$n_decode_preferred, chosen6_n = c6$n_decode_preferred + c6$n_not,
             other37_pref = o$n_decode_preferred, other37_n = o$n_decode_preferred + o$n_not,
             all43_pref = al$n_decode_preferred, all43_n = al$n_decode_preferred + al$n_not, all43_NA = al$n_NA,
             odds_ratio = unname(ft$estimate), p_two_sided = ft$p.value)
}))
write.table(out, a[2], sep = "\t", quote = FALSE, row.names = FALSE)
print(out)
