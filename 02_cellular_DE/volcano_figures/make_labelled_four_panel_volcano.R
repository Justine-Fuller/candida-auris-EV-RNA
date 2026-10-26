## ============================================================================
## FINAL four-panel, gene-labelled DESeq2 volcano summary figure
## A. Cell sRNA | B. Cell mRNA | C. B8441 EV RNA | D. MMC1 EV RNA
## ============================================================================
## Permanent, reproducible script (project-relative paths via PROJECT_ROOT
## below -- edit only that line if the project folder is ever moved).
##
## Reads ONLY the existing, already-finalised DESeq2 result CSVs -- does NOT
## rerun DESeq2, does NOT alter any DESeq2 result. Labels exactly the 24
## genes (6 per panel) approved in the prior label-selection audit; log2FC
## and padj for the label table below are pulled live from the canonical
## result files at run time (never hand-typed), so the labelled-gene table
## and the plotted points can never drift out of sync with the source data.
##
## Significance convention (identical to every other DE figure in this
## project): padj < 0.05. NA padj is never classified as significant.
## No fold-change threshold is used for the significance colouring.
##
## Positive log2FoldChange = higher abundance in Treated relative to
## Untreated (contrast = c("Treatment","Treated","Untreated"), Untreated
## releveled as reference in all four underlying DESeq2 analyses).
##
## Axis limits: PANEL-SPECIFIC, not shared across panels. Ranges inspected
## before plotting:
##   A (Cell sRNA):  log2FC [-5.3,  3.9], max -log10(padj)  56.6
##   B (Cell mRNA):  log2FC [-8.3,  7.0], max -log10(padj) 245.6
##   C (B8441 EV):   log2FC [-7.0, 14.9], max -log10(padj) 155.4
##   D (MMC1 EV):    log2FC [-9.5, 16.0], max -log10(padj) 195.5
## Forcing identical limits across all four panels would compress panel A
## (both axes) and panel C (x-axis) into a small fraction of the plotting
## area relative to panels B/D; panel-specific limits were used instead so
## every panel's own distribution is legible. Stated explicitly in the
## figure caption.
##
## Outputs (new files; nothing pre-existing is overwritten):
##   ../FINAL_DESeq2_four_panel_labelled_volcano.pdf
##   ../FINAL_DESeq2_four_panel_labelled_volcano.png
##   ./labelled_genes_table.csv   (Panel, Gene_ID, plot_label, log2FoldChange,
##                                 padj, Direction, annotation_source)
## ============================================================================

## The plot labels use a few Unicode characters (superscript "+", beta,
## prime). R on this machine defaults to the "C" (ASCII-only) locale, which
## corrupts multibyte UTF-8 text at parse time -- set UTF-8 explicitly before
## anything else runs so the labels are read and rendered correctly
## regardless of the caller's shell environment.
invisible(tryCatch(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), error = function(e) NULL))

suppressMessages({
  library(tidyverse)
  library(ggrepel)
  library(patchwork)
})

PROJECT_ROOT <- "/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV"
DE  <- file.path(PROJECT_ROOT, "Differential_expression_analysis")
OUT <- file.path(DE, "00_Summary_figures")
SCRIPTS_DIR <- file.path(OUT, "scripts")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(SCRIPTS_DIR, showWarnings = FALSE, recursive = TRUE)

files <- list(
  A = list(path = file.path(DE, "03_Cellular_DE_analysis/FINAL_Cell_sRNA_transcript_only_harmonized/02_Results/FINAL_Cell_sRNA_DESeq2_Results.csv"),
           title = "A. Cell sRNA", annotation_source = "full_proteome_annotation.emapper.annotations"),
  B = list(path = file.path(DE, "03_Cellular_DE_analysis/mRNA_DE_analysis/02_Results/03_mRNA_Cells_DESeq2_Results.csv"),
           title = "B. Cell mRNA", annotation_source = "full_proteome_annotation.emapper.annotations"),
  C = list(path = file.path(DE, "04_EV_DE_analysis/01_mRNA_EV_analysis/01_Results/EV_DESeq2_B8441_Treated_vs_Untreated.csv"),
           title = "C. B8441 EV RNA", annotation_source = "ev_all_annotation.emapper.csv"),
  D = list(path = file.path(DE, "04_EV_DE_analysis/01_mRNA_EV_analysis/01_Results/EV_DESeq2_MMC1_Treated_vs_Untreated.csv"),
           title = "D. MMC1 EV RNA", annotation_source = "ev_all_annotation.emapper.csv")
)
for (p in names(files)) stopifnot(file.exists(files[[p]]$path))

## ---- the 24 approved gene IDs / labels (exact text, approved in the prior audit) ----
labels <- tribble(
  ~Panel, ~Gene_ID, ~plot_label,
  "A","B9J08_003046","PLP-dependent enzyme (CHA1)",
  "A","B9J08_005075","Energy-metabolism regulator (RGI1)",
  "A","B9J08_003088","Aldehyde dehydrogenase (ALD5)",
  "A","B9J08_002279","Carbon-utilization TF (GSM1)",
  "A","B9J08_000483","Hsp70 chaperone (SSA3)",
  "A","B9J08_002073","Histone H4 (HHF1)",
  "B","B9J08_002202","Na⁺/phosphate symporter (PHO89)",
  "B","B9J08_001464","Phosphate transporter (PHO84)",
  "B","B9J08_002660","MFS sugar transporter",
  "B","B9J08_001366","Cell wall protein (CCW12)",
  "B","B9J08_002073","Histone H4 (HHF1)",
  "B","B9J08_001487","Siderophore-iron transporter (ARN family)",
  "C","B9J08_002066_transcript","Git3/GPR_Gpa2_C-domain protein (GPR1)",
  "C","B9J08_000918_transcript","β-1,3-glucan-elongating enzyme (PHR1)",
  "C","B9J08_004489_transcript","COPI coatomer subunit (SEC21)",
  "C","B9J08_001105_transcript","MCM-family protein (MCM6)",
  "C","B9J08_004496_transcript","Aldo-keto reductase (plr1)",
  "C","B9J08_000872_transcript","COG-complex transport protein (COG4)",
  "D","B9J08_002066_transcript","Git3/GPR_Gpa2_C-domain protein (GPR1)",
  "D","B9J08_000918_transcript","β-1,3-glucan-elongating enzyme (PHR1)",
  "D","B9J08_004496_transcript","Aldo-keto reductase (plr1)",
  "D","B9J08_001105_transcript","MCM-family protein (MCM6)",
  "D","B9J08_001677_transcript","V-ATPase subunit (STV1)",
  "D","B9J08_001879_transcript","tRNA 2′-O-ribose methyltransferase (TRM7)"
)
stopifnot(nrow(labels) == 24)
stopifnot(all(table(labels$Panel) == 6))

SIG_COL <- "#C41E1E"      # red, significant
NONSIG_COL <- "grey70"    # grey, not significant
## "Helvetica Neue" is used (rather than the generic "sans" alias) because it
## is confirmed, by direct test render, to include the superscript-plus
## glyph (U+207A) used in one label ("Na⁺..."); the generic sans
## fallback on this machine does not.
FONT_FAMILY <- "Helvetica Neue"

theme_volcano <- theme_bw(base_size = 11, base_family = FONT_FAMILY) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "grey93", linewidth = 0.3),
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    axis.title = element_text(size = 10.5),
    axis.text = element_text(size = 9, colour = "black"),
    legend.position = "none",
    plot.margin = margin(6, 10, 4, 4)
  )

## builds the panel plot AND returns the resolved label rows (with live
## log2FoldChange/padj pulled from this same data frame) so the permanent
## labelled_genes_table.csv can never drift from what is actually plotted
panel_result <- list()

## Panel D (MMC1 EV RNA): all 6 labelled points sit low on the y-axis
## (padj not as extreme as in panel C) while the panel's y-range extends to
## ~200, leaving a lot of empty space above them unused by default ggrepel
## placement -- which, combined with default repel strength, left labels
## overlapping and one label clipped on first render. A small upward nudge
## plus a much stronger repel force fixes this (confirmed by isolated
## test render); default (no nudge, standard force) is kept for A-C, which
## rendered cleanly without it.
NUDGE_Y <- c(A = 0, B = 0, C = 0, D = 40)

build_panel <- function(panel_id) {
  info <- files[[panel_id]]
  d <- read.csv(info$path, stringsAsFactors = FALSE)
  colnames(d)[1] <- "Feature"
  d$is_sig <- !is.na(d$padj) & d$padj < 0.05
  d$Sig <- factor(ifelse(d$is_sig, "Significant (padj < 0.05)", "Not significant"),
                   levels = c("Significant (padj < 0.05)", "Not significant"))
  d$neglog10padj <- -log10(d$padj)

  panel_labels <- labels %>% filter(Panel == panel_id)
  idx <- match(panel_labels$Gene_ID, d$Feature)
  stopifnot(!any(is.na(idx)))   # every approved Gene_ID must resolve -- stop hard if not

  lab_d <- d[idx, ]
  lab_d$plot_label <- panel_labels$plot_label
  lab_d$Panel <- panel_id
  lab_d$annotation_source <- info$annotation_source
  lab_d$Direction <- ifelse(lab_d$log2FoldChange > 0, "Up in Treated", "Down in Treated")

  panel_result[[panel_id]] <<- lab_d %>%
    select(Panel, Gene_ID = Feature, plot_label, log2FoldChange, padj, Direction, annotation_source)

  ggplot(d, aes(x = log2FoldChange, y = neglog10padj, colour = Sig)) +
    geom_point(size = 1.1, alpha = 0.5, na.rm = TRUE) +
    scale_colour_manual(values = c("Significant (padj < 0.05)" = SIG_COL, "Not significant" = NONSIG_COL),
                         name = NULL, drop = FALSE) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), colour = "grey40", linewidth = 0.4, linetype = "dashed") +
    geom_point(data = lab_d, aes(x = log2FoldChange, y = neglog10padj), inherit.aes = FALSE,
               colour = "black", size = 1.6, shape = 21, fill = SIG_COL, stroke = 0.5) +
    geom_label_repel(
      data = lab_d, aes(x = log2FoldChange, y = neglog10padj, label = plot_label), inherit.aes = FALSE,
      size = 2.7, family = FONT_FAMILY, colour = "black",
      fill = alpha("white", 0.92), label.size = 0.15,
      segment.size = 0.3, segment.colour = "grey30", min.segment.length = 0,
      box.padding = 0.6, point.padding = 0.3, force = 60, force_pull = 0.2,
      max.overlaps = Inf, seed = 42, max.iter = 50000, max.time = 8,
      direction = "both", nudge_y = NUDGE_Y[[panel_id]]
    ) +
    labs(title = info$title, x = expression(log[2]~"fold change (Treated vs Untreated)"),
         y = expression(-log[10]~"(adjusted p-value)")) +
    theme_volcano
}

pA <- build_panel("A")
pB <- build_panel("B")
pC <- build_panel("C")
pD <- build_panel("D")

## ---- permanent labelled-gene table (values pulled live from source above) ----
labelled_genes_table <- bind_rows(panel_result[["A"]], panel_result[["B"]], panel_result[["C"]], panel_result[["D"]])
stopifnot(nrow(labelled_genes_table) == 24)
table_path <- file.path(SCRIPTS_DIR, "labelled_genes_table.csv")
write.csv(labelled_genes_table, table_path, row.names = FALSE)
cat("Wrote labelled-gene table:", table_path, "\n")

## ---- assemble 2x2 layout with ONE shared legend, kept OUTSIDE all panels ----
## (avoids the earlier unlabelled-figure issue where an in-panel legend
## briefly overlapped real data points)
combined <- (pA + pB) / (pC + pD) +
  plot_layout(guides = "collect") &
  theme(legend.position = "top", legend.justification = "right",
        legend.text = element_text(size = 9), legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm"), legend.margin = margin(0, 0, 2, 0))

W <- 11; H <- 9.5

pdf_path <- file.path(OUT, "FINAL_DESeq2_four_panel_labelled_volcano.pdf")
png_path <- file.path(OUT, "FINAL_DESeq2_four_panel_labelled_volcano.png")

## cairo_pdf requires X11/cairo system libraries that are not installed on
## this machine (XQuartz absent), and ragg (used for the PNG below) has no
## vector PDF device. Use macOS's native Quartz PDF device instead --
## dependency-free on this platform and, combined with the locale fix and
## font choice above, confirmed by direct test render to draw all label
## glyphs (including the superscript "+") correctly.
quartz_pdf_device <- function(filename, width, height, ...) {
  grDevices::quartz(type = "pdf", file = filename, width = width, height = height)
}
ggsave(pdf_path, combined, width = W, height = H, device = quartz_pdf_device)
ggsave(png_path, combined, width = W, height = H, dpi = 400, bg = "white", device = ragg::agg_png)

cat("Saved:\n ", pdf_path, "\n ", png_path, "\n")
cat("Figure dimensions:", W, "x", H, "in (PNG at 400 dpi)\n")
