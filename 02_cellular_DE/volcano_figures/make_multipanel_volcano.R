## ============================================================================
## Multi-panel volcano summary figure: Cell sRNA / Cell mRNA / EV RNA (B8441, MMC1)
## ============================================================================
## Reads only EXISTING, already-saved DESeq2 result files. Does NOT rerun any
## DESeq2 analysis. Panel structure is A/B/C/D (4 panels) because the EV
## differential-expression analysis was performed separately per strain
## (B8441, MMC1) -- no combined-strain EV result was ever saved, so it is not
## legitimate to show EV as a single panel.
##
## Source files for each panel (verified against Results.docx headline numbers
## before use -- see accompanying audit):
##   A. Cell sRNA   : Differential_expression_analysis/03_Cellular_DE_analysis/
##                    FINAL_Cell_sRNA_transcript_only_harmonized/02_Results/
##                    FINAL_Cell_sRNA_DESeq2_Results.csv
##                    (THE FINAL, STAR-harmonized, name-safe, transcript-only
##                    analysis -- NOT the historical or provisional/diagnostic
##                    Cell sRNA results)
##   B. Cell mRNA   : Differential_expression_analysis/03_Cellular_DE_analysis/
##                    mRNA_DE_analysis/02_Results/03_mRNA_Cells_DESeq2_Results.csv
##   C. B8441 EV RNA: Differential_expression_analysis/04_EV_DE_analysis/
##                    01_mRNA_EV_analysis/01_Results/
##                    EV_DESeq2_B8441_Treated_vs_Untreated.csv
##   D. MMC1 EV RNA : Differential_expression_analysis/04_EV_DE_analysis/
##                    01_mRNA_EV_analysis/01_Results/
##                    EV_DESeq2_MMC1_Treated_vs_Untreated.csv
##
## All four analyses use design ~Treatment (Cell mRNA: ~Treatment*Strain, main
## Treatment effect extracted), Untreated as the releveled reference, and
## contrast = c("Treatment","Treated","Untreated") -- so positive log2FC means
## higher abundance in Treated in every panel. Verified empirically for Cell
## sRNA (normalized counts checked directly); verified from the actual
## saved contrast/releveling code for Cell mRNA and both EV analyses (no
## DESeqDataSet survives for those three to re-check empirically without
## rerunning DESeq2, which was out of scope for this figure).
##
## All four panels use the SAME significance criterion: padj < 0.05, with no
## fold-change threshold baked into "significant" -- confirmed identical
## across all four source analyses before this script was written, so no
## values are silently harmonised here.
##
## Output (both written to the permanent project tree, NOT scratchpad):
##   Differential_expression_analysis/00_Summary_figures/
##     DESeq2_volcano_summary_ABCD.pdf   (vector, for the thesis)
##     DESeq2_volcano_summary_ABCD.png   (high-res raster, for quick preview)
## ============================================================================

CADIR <- "/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV"
DE <- file.path(CADIR, "Differential_expression_analysis")
OUT <- file.path(DE, "00_Summary_figures")
dir.create(OUT, showWarnings = FALSE)

files <- list(
  A = list(path = file.path(DE, "03_Cellular_DE_analysis/FINAL_Cell_sRNA_transcript_only_harmonized/02_Results/FINAL_Cell_sRNA_DESeq2_Results.csv"),
           title = "A. Cell sRNA"),
  B = list(path = file.path(DE, "03_Cellular_DE_analysis/mRNA_DE_analysis/02_Results/03_mRNA_Cells_DESeq2_Results.csv"),
           title = "B. Cell mRNA"),
  C = list(path = file.path(DE, "04_EV_DE_analysis/01_mRNA_EV_analysis/01_Results/EV_DESeq2_B8441_Treated_vs_Untreated.csv"),
           title = "C. B8441 EV RNA"),
  D = list(path = file.path(DE, "04_EV_DE_analysis/01_mRNA_EV_analysis/01_Results/EV_DESeq2_MMC1_Treated_vs_Untreated.csv"),
           title = "D. MMC1 EV RNA")
)

SIG_COL   <- rgb(0.80, 0.10, 0.10, 0.55)   # red, significant
NONSIG_COL <- rgb(0.55, 0.55, 0.55, 0.30)  # grey, not significant
SIG_COL_SOLID <- rgb(0.80, 0.10, 0.10, 0.9)
NONSIG_COL_SOLID <- rgb(0.55, 0.55, 0.55, 0.9)
FONT_FAMILY <- "sans"

draw_panel <- function(d, panel_title) {
  is_sig <- !is.na(d$padj) & d$padj < 0.05
  col <- ifelse(is_sig, SIG_COL, NONSIG_COL)
  xr <- range(d$log2FoldChange, na.rm = TRUE)
  yr <- range(-log10(d$padj), na.rm = TRUE, finite = TRUE)
  plot(d$log2FoldChange, -log10(d$padj), pch = 16, cex = 0.55, col = col,
       xlab = expression(log[2]~"fold change"), ylab = expression(-log[10]~"(adjusted p)"),
       family = FONT_FAMILY, cex.lab = 1.0, cex.axis = 0.85, font.main = 1,
       xlim = xr, ylim = yr)
  abline(h = -log10(0.05), col = "red", lty = 2, lwd = 0.8)
  title(main = panel_title, adj = 0, cex.main = 1.15, family = FONT_FAMILY, font.main = 2)
  n_tested <- nrow(d)
  n_sig <- sum(is_sig)
  mtext(sprintf("n = %d tested, %d significant", n_tested, n_sig),
        side = 3, line = 0.2, adj = 0, cex = 0.65, family = FONT_FAMILY, col = "grey30")
}

render_figure <- function() {
  # reserve a dedicated top outer-margin strip for the shared legend, so it
  # sits ABOVE all four panels (still the upper-right of the overall figure)
  # and can never overlap any panel's data, regardless of that panel's own
  # point distribution.
  par(mfrow = c(2, 2), mar = c(4.2, 4.2, 3.2, 1.5), oma = c(0.5, 0.5, 2.6, 0.5),
      family = FONT_FAMILY, bg = "white")
  for (p in names(files)) {
    d <- read.csv(files[[p]]$path, stringsAsFactors = FALSE)
    draw_panel(d, files[[p]]$title)
  }
  # single shared legend for the whole figure, placed in the reserved top
  # margin strip, right-aligned -- outside every panel's plot region
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  legend(x = 0.99, y = 1.0, xjust = 1, yjust = 1,
         legend = c("Significant (padj < 0.05)", "Not significant"),
         col = c(SIG_COL_SOLID, NONSIG_COL_SOLID), pch = 16, pt.cex = 1.2,
         cex = 0.85, bty = "o", bg = "white", box.col = "grey50",
         xpd = NA, horiz = FALSE, x.intersp = 0.8, seg.len = 1)
}

## ---- PDF (vector, thesis) ----
pdf(file.path(OUT, "DESeq2_volcano_summary_ABCD.pdf"), width = 9, height = 8)
render_figure()
dev.off()

## ---- PNG (high-res raster, preview) ----
png(file.path(OUT, "DESeq2_volcano_summary_ABCD.png"), width = 9, height = 8, units = "in", res = 300)
render_figure()
dev.off()

cat("Saved:\n  ", file.path(OUT, "DESeq2_volcano_summary_ABCD.pdf"), "\n  ",
    file.path(OUT, "DESeq2_volcano_summary_ABCD.png"), "\n")

## ---- print headline numbers actually used, for the record ----
for (p in names(files)) {
  d <- read.csv(files[[p]]$path, stringsAsFactors = FALSE)
  cat(sprintf("%s: tested=%d valid_padj=%d sig(padj<0.05)=%d  [%s]\n",
      files[[p]]$title, nrow(d), sum(!is.na(d$padj)), sum(d$padj < 0.05, na.rm = TRUE), files[[p]]$path))
}
