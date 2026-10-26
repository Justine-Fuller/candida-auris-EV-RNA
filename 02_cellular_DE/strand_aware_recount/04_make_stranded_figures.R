## ============================================================================
## Figures for the STRANDED cellular mRNA results (2026-10-14)
##   1. Four-panel volcano (same as FINAL_DESeq2_four_panel_labelled_volcano),
##      with panel B switched to the stranded B8441 cell mRNA result.
##      A, C, D are read unchanged from their existing result files.
##   2. New cell mRNA figure: A. B8441 volcano  B. MMC1 volcano  C. PCA
## Same styling, significance rule (padj < 0.05, NA never significant) and
## panel-B gene labels as make_labelled_four_panel_volcano.R. log2FC/padj for
## the labels are pulled live from the results files.
## Outputs go to 00_Summary_figures/ with new names; nothing is overwritten.
## ============================================================================
invisible(tryCatch(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), warning = function(w) Sys.setlocale("LC_CTYPE", "C.UTF-8"),
                   error = function(e) NULL))
suppressMessages({ library(tidyverse); library(ggrepel); library(patchwork); library(DESeq2) })

PROJECT_ROOT <- Sys.getenv("CANDIDAEV_ROOT",
  "/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV")
FONT_FAMILY <- Sys.getenv("FIG_FONT", "Helvetica Neue")
DE  <- file.path(PROJECT_ROOT, "Differential_expression_analysis")
OUT <- file.path(DE, "00_Summary_figures")
STR <- file.path(DE, "03_Cellular_DE_analysis/mRNA_DE_analysis/STRANDED_recount_2026-10-14/results")

SIG_COL <- "#C41E1E"; NONSIG_COL <- "grey70"
STRAIN_COL <- c(B8441 = "#2a78d6", MMC1 = "#eb6834")   # validated categorical pair (CVD dE 24.7)

theme_volcano <- theme_bw(base_size = 11, base_family = FONT_FAMILY) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey93", linewidth = 0.3),
        plot.title = element_text(face = "bold", size = 12, hjust = 0),
        axis.title = element_text(size = 10.5),
        axis.text = element_text(size = 9, colour = "black"),
        legend.position = "none", plot.margin = margin(6, 10, 4, 4))

cell_mrna_labels <- tribble(
  ~Gene_ID, ~plot_label,
  "B9J08_002202", "Na⁺/phosphate symporter (PHO89)",
  "B9J08_001464", "Phosphate transporter (PHO84)",
  "B9J08_002660", "MFS sugar transporter",
  "B9J08_001366", "Cell wall protein (CCW12)",
  "B9J08_002073", "Histone H4 (HHF1)",
  "B9J08_001487", "Siderophore-iron transporter (ARN family)")
old_labels <- tribble(
  ~Panel, ~Gene_ID, ~plot_label,
  "A","B9J08_003046","PLP-dependent enzyme (CHA1)",
  "A","B9J08_005075","Energy-metabolism regulator (RGI1)",
  "A","B9J08_003088","Aldehyde dehydrogenase (ALD5)",
  "A","B9J08_002279","Carbon-utilization TF (GSM1)",
  "A","B9J08_000483","Hsp70 chaperone (SSA3)",
  "A","B9J08_002073","Histone H4 (HHF1)",
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
  "D","B9J08_001879_transcript","tRNA 2′-O-ribose methyltransferase (TRM7)")

volcano <- function(path, title, labs_tbl, nudge_y = 0) {
  d <- read.csv(path, stringsAsFactors = FALSE); colnames(d)[1] <- "Feature"
  d$Sig <- factor(ifelse(!is.na(d$padj) & d$padj < 0.05, "Significant (padj < 0.05)", "Not significant"),
                  levels = c("Significant (padj < 0.05)", "Not significant"))
  d$y <- -log10(d$padj)
  idx <- match(labs_tbl$Gene_ID, d$Feature); stopifnot(!any(is.na(idx)))
  lab <- d[idx, ]; lab$plot_label <- labs_tbl$plot_label
  ggplot(d, aes(log2FoldChange, y, colour = Sig)) +
    geom_point(size = 1.1, alpha = 0.5, na.rm = TRUE) +
    scale_colour_manual(values = c("Significant (padj < 0.05)" = SIG_COL, "Not significant" = NONSIG_COL),
                        name = NULL, drop = FALSE) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), colour = "grey40", linewidth = 0.4, linetype = "dashed") +
    geom_point(data = lab, aes(log2FoldChange, y), inherit.aes = FALSE,
               colour = "black", size = 1.6, shape = 21, fill = SIG_COL, stroke = 0.5) +
    geom_label_repel(data = lab, aes(log2FoldChange, y, label = plot_label), inherit.aes = FALSE,
      size = 2.7, family = FONT_FAMILY, colour = "black", fill = alpha("white", 0.92), label.size = 0.15,
      segment.size = 0.3, segment.colour = "grey30", min.segment.length = 0,
      box.padding = 0.6, point.padding = 0.3, force = 60, force_pull = 0.2,
      max.overlaps = Inf, seed = 42, max.iter = 50000, max.time = 8, direction = "both", nudge_y = nudge_y) +
    labs(title = title, x = expression(log[2]~"fold change (Treated vs Untreated)"),
         y = expression(-log[10]~"(adjusted p-value)")) +
    theme_volcano
}
collect <- function(p) p + plot_layout(guides = "collect") &
  theme(legend.position = "top", legend.justification = "right",
        legend.text = element_text(size = 9), legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm"), legend.margin = margin(0, 0, 2, 0))

## ---- 1. four-panel figure, panel B = stranded B8441 ----
pA <- volcano(file.path(DE, "03_Cellular_DE_analysis/FINAL_Cell_sRNA_transcript_only_harmonized/02_Results/FINAL_Cell_sRNA_DESeq2_Results.csv"),
              "A. Cell sRNA", filter(old_labels, Panel == "A"))
pB <- volcano(file.path(STR, "new_stranded_B8441_results.csv"), "B. Cell mRNA (B8441)", cell_mrna_labels)
pC <- volcano(file.path(DE, "04_EV_DE_analysis/01_mRNA_EV_analysis/01_Results/EV_DESeq2_B8441_Treated_vs_Untreated.csv"),
              "C. B8441 EV RNA", filter(old_labels, Panel == "C"))
pD <- volcano(file.path(DE, "04_EV_DE_analysis/01_mRNA_EV_analysis/01_Results/EV_DESeq2_MMC1_Treated_vs_Untreated.csv"),
              "D. MMC1 EV RNA", filter(old_labels, Panel == "D"), nudge_y = 40)
four <- collect((pA + pB) / (pC + pD))
ggsave(file.path(OUT, "FINAL_DESeq2_four_panel_labelled_volcano_STRANDED_2026-10-14.png"),
       four, width = 11, height = 9.5, dpi = 400, bg = "white", device = ragg::agg_png)

## ---- 2. cell mRNA figure: both strains + PCA ----
vB <- volcano(file.path(STR, "new_stranded_B8441_results.csv"), "A. Cell mRNA, B8441", cell_mrna_labels)
## MMC1 labels sit low on a tall y-axis; nudge them up (as for EV panel D) so none overlap or clip
vM <- volcano(file.path(STR, "new_stranded_MMC1_results.csv"),  "B. Cell mRNA, MMC1",  cell_mrna_labels, nudge_y = 45)

dds <- readRDS(file.path(STR, "new_stranded_dds.rds"))
vsd <- vst(dds, blind = TRUE)
pca <- plotPCA(vsd, intgroup = c("Strain", "Treatment"), ntop = 500, returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"), 1)
pca$Treatment <- factor(pca$Treatment, levels = c("Untreated", "Treated"))
pP <- ggplot(pca, aes(PC1, PC2, colour = Strain, shape = Treatment)) +
  geom_point(size = 3, stroke = 0.9) +
  scale_colour_manual(values = STRAIN_COL) +
  scale_shape_manual(values = c(Untreated = 1, Treated = 16)) +
  labs(title = "C. PCA (top 500 variable genes)",
       x = paste0("PC1 (", pv[1], "% variance)"), y = paste0("PC2 (", pv[2], "% variance)")) +
  theme_volcano + theme(legend.position = "right", legend.title = element_text(size = 9),
                        legend.text = element_text(size = 9), aspect.ratio = 1)

cell <- vB + vM + pP + plot_layout(ncol = 3, widths = c(1, 1, 0.9), guides = "collect") &
  theme(legend.position = "bottom", legend.title = element_blank())
ggsave(file.path(OUT, "Cell_mRNA_STRANDED_volcano_PCA_2026-10-14.png"),
       cell, width = 14, height = 5.2, dpi = 400, bg = "white", device = ragg::agg_png)
cat("Saved both figures to", OUT, "\n")
