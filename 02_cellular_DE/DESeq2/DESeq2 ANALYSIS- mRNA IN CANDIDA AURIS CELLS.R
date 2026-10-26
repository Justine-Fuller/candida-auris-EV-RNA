# ============================================================================
# DESeq2 ANALYSIS: mRNA IN CANDIDA AURIS CELLS
# Caspofungin Treatment Response
# ============================================================================
# This script identifies messenger RNAs (mRNA) that are differentially expressed
# between untreated and caspofungin-treated C. auris cells
# mRNAs are the protein-coding transcripts

# ============================================================================
# STEP 1: LOAD REQUIRED LIBRARIES AND DATA
# ============================================================================

library(DESeq2)      # Main package for differential expression analysis
library(tidyverse)   # Data manipulation and visualization
library(pheatmap)    # Heatmap visualization

# Load the transcript-level count matrix and metadata from the setup script
# These RDS files were created in the data preparation step
# 
# CHOOSE ONE OF THE FOLLOWING:
# Option 1: For TRANSCRIPT-level analysis (uncomment below):
 count_matrix <- readRDS("/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline/count_matrix_transcripts.rds")
 cat("=== mRNA CELLS ANALYSIS - TRANSCRIPT LEVEL ===\n")

# Option 2: For NON-CODING FEATURE analysis (exons, introns, UTRs) (uncomment below):
# count_matrix <- readRDS("/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline/count_matrix_ncfeatures.rds")
# cat("=== mRNA CELLS ANALYSIS - NON-CODING FEATURES ===\n")

# Load sample metadata (same for both analyses)
sample_metadata <- readRDS("/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline/sample_metadata.rds")

cat("Loaded count matrix with", nrow(count_matrix), "features\n")
cat("Loaded metadata with", nrow(sample_metadata), "samples\n")

# ============================================================================
# STEP 2: FILTER SAMPLES FOR THIS ANALYSIS
# ============================================================================
# We want to subset the data to include ONLY:
#   - RNA_type == "mRNA" (messenger RNA / protein-coding)
#   - Compartment == "Cell" (cells, not extracellular vesicles)

# Create a logical vector indicating which samples match our criteria
samples_to_keep <- (sample_metadata$RNA_type == "mRNA") & 
  (sample_metadata$Compartment == "Cell")

# Display which samples will be included
cat("\n=== SAMPLE SELECTION ===\n")
cat("Criteria: mRNA in Cells\n")
cat("Samples selected:", sum(samples_to_keep), "\n")
cat("Sample IDs:\n")
print(rownames(sample_metadata)[samples_to_keep])

# Subset count matrix to keep only selected samples
# The double bracket syntax [rows, columns] keeps all features (rows) 
# and only the selected samples (columns)
count_matrix_subset <- count_matrix[, samples_to_keep]

# Subset metadata to match the count matrix
metadata_subset <- sample_metadata[samples_to_keep, ]

cat("\nSubset count matrix dimensions:", dim(count_matrix_subset), "\n")
cat("Subset metadata dimensions:", dim(metadata_subset), "\n")

# Display the metadata for this analysis
cat("\nSample metadata for this analysis:\n")
print(metadata_subset)

# ============================================================================
# STEP 3: FILTER OUT FEATURES WITH NO COUNTS IN THIS SUBSET
# ============================================================================
# After subsetting samples, some features may have zero counts across all
# samples in this subset. We need to remove these to avoid errors.

# Count total reads per feature across all selected samples
# Convert to matrix to ensure all values are numeric
count_matrix_numeric <- as.matrix(count_matrix_subset)
feature_totals <- rowSums(count_matrix_numeric, na.rm = TRUE)

# Keep only features with at least 1 read total
features_to_keep <- feature_totals > 0

cat("\n=== FEATURE FILTERING FOR SUBSET ===\n")
cat("Features with at least 1 read in this subset:", sum(features_to_keep, na.rm = TRUE), "\n")
cat("Features with zero reads (removed):", sum(!features_to_keep, na.rm = TRUE), "\n")

# Apply filter
count_matrix_subset <- count_matrix_subset[features_to_keep, ]

cat("Final feature count for analysis:", nrow(count_matrix_subset), "\n")

# ============================================================================
# STEP 4: PREPARE DATA FOR DESeq2
# ============================================================================
# DESeq2 requires a specific data structure called a DESeqDataSet (dds)
# This object combines:
#   - Count matrix (raw integer counts)
#   - Sample metadata (phenotype information)
#   - Design formula (statistical model for the analysis)

# INTERACTION MODEL (Treatment * Strain)
# Tests for Treatment effect, Strain effect, AND their interaction
# Allows Treatment effect to differ between strains
 design_formula <- ~Treatment * Strain
 cat("\n=== DESIGN FORMULA: INTERACTION MODEL ===\n")
 cat("Formula: ~Treatment * Strain\n")
 cat("This tests if Treatment effect differs by Strain\n")

# Create the DESeqDataSet object
# colData = sample metadata (information about each column/sample)
# design = the statistical model to fit
dds <- DESeqDataSetFromMatrix(
  # Count matrix - must have integer values, not decimals
  countData = count_matrix_subset,
  # Sample metadata - describes each sample
  colData = metadata_subset,
  # Design formula - specifies the statistical model
  design = design_formula
)

# Set the reference levels for variables
# This determines which condition is the "control" or "baseline"
dds$Treatment <- factor(dds$Treatment, levels = c("Untreated", "Treated"))
dds$Strain <- factor(dds$Strain, levels = c("B8441", "MMC1"))

cat("\n=== DESeqDataSet CREATED ===\n")
cat("Design formula:", as.character(design(dds)), "\n")
cat("Treatment levels:", levels(dds$Treatment), "\n")
cat("Strain levels:", levels(dds$Strain), "\n")
cat("DESeqDataSet dimensions:", dim(dds), "\n")

# ============================================================================
# STEP 5: RUN DESeq2 ANALYSIS
# ============================================================================
# The DESeq function performs three key steps:
#   1. Estimation of size factors (library size normalization)
#   2. Estimation of dispersions (variance-mean relationship)
#   3. Negative binomial generalized linear model fitting and Wald testing

cat("\n=== RUNNING DESeq2 ===\n")

# Run the DESeq2 analysis pipeline
# This function automatically runs: estimateSizeFactors, estimateDispersions, 
# and nbinomWaldTest
dds <- DESeq(dds)

cat("✓ DESeq2 analysis complete\n")

# ============================================================================
# STEP 6: EXTRACT AND ORGANIZE RESULTS
# ============================================================================

# Extract results from the DESeq2 object
# contrast specifies which comparison to extract
# c("Treatment", "Treated", "Untreated") means: Treated vs Untreated
results <- results(
  dds,
  # Specify the contrast: compare Treated to Untreated
  contrast = c("Treatment", "Treated", "Untreated"),
  # Calculate p-values using Wald test
  test = "Wald"
)

# Convert results to a dataframe and add feature names
results_df <- as.data.frame(results)
results_df$Feature <- rownames(results_df)

# Reorder columns for clarity
results_df <- results_df %>%
  mutate(Feature = str_remove(Feature, "_transcript$")) %>%
  # Select columns in this order: Feature, baseMean, log2FoldChange, etc.
  select(Feature, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
  # Sort by adjusted p-value (most significant first)
  arrange(padj)

cat("\n=== RESULTS SUMMARY ===\n")
cat("Total features tested:", nrow(results_df), "\n")

# Count significant features at different thresholds
sig_005 <- sum(results_df$padj < 0.05, na.rm = TRUE)
sig_001 <- sum(results_df$padj < 0.01, na.rm = TRUE)
sig_upregulated <- sum(results_df$padj < 0.05 & results_df$log2FoldChange > 0, na.rm = TRUE)
sig_downregulated <- sum(results_df$padj < 0.05 & results_df$log2FoldChange < 0, na.rm = TRUE)

cat("Features significant at padj < 0.05:", sig_005, "\n")
cat("Features significant at padj < 0.01:", sig_001, "\n")
cat("  - Upregulated (in Treated):", sig_upregulated, "\n")
cat("  - Downregulated (in Treated):", sig_downregulated, "\n")

# Display top 10 differentially expressed features
cat("\n=== TOP 10 DIFFERENTIALLY EXPRESSED mRNAs ===\n")
print(head(results_df, 10))

# Display top 10 differentially expressed lncRNAs
cat("\n=== TOP 10 DIFFERENTIALLY EXPRESSED lncRNAs ===\n")

# Filter results for lncRNA features (rownames starting with LNCRNA)
lncrna_results <- results_df %>%
  filter(grepl("^LNCRNA", Feature)) %>%
  # Sort by adjusted p-value (most significant first)
  arrange(padj)

if (nrow(lncrna_results) > 0) {
  print(head(lncrna_results, 10))
} else {
  cat("No lncRNAs found in results\n")
}

# Display top 10 differentially expressed coding transcripts (excludes lncRNAs)
cat("\n=== TOP 10 DIFFERENTIALLY EXPRESSED CODING mRNAs ===\n")
# Filter results to exclude lncRNA features (rownames starting with LNCRNA)
coding_results <- results_df %>%
  filter(!grepl("^LNCRNA", Feature)) %>%
  # Sort by adjusted p-value (most significant first)
  arrange(padj)

if (nrow(coding_results) > 0) {
  print(head(coding_results, 10))
} else {
  cat("No coding transcripts found in results\n")
}

# ============================================================================
# STEP 7: QUALITY CONTROL VISUALIZATIONS
# ============================================================================

# Create output directory if it doesn't exist
dir.create("/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline/mRNA_CELL", showWarnings = FALSE)

# PLOT 1: MA-plot (M vs A plot)
# This shows log fold change (M-axis) vs mean expression (A-axis)
# Red points are significant features (padj < 0.05)
pdf("/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline/03_mRNA_Cells_MA_plot.pdf", width = 8, height = 6)
plotMA(
  results,
  # Color points with adjusted p-value < 0.05 in red
  alpha = 0.05,
  # Main title for the plot
  main = "mRNA in C. auris Cells: MA Plot\n(Treated vs Untreated)"
)
dev.off()
cat("\n✓ Saved MA-plot to: 03_mRNA_Cells_MA_plot.pdf\n")

# PLOT 2: Volcano plot
# This shows log fold change (X-axis) vs -log10(padj) (Y-axis)
# Points farther right or left have larger fold changes
# Points higher up have lower (more significant) p-values
pdf("/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline/mRNA_CELL/nc_features_Volcano_plot_Treatment+Strain.pdf", width = 8, height = 6)

# Create the volcano plot manually for more control
plot(
  # X-axis: log2 fold change
  x = results_df$log2FoldChange,
  # Y-axis: -log10 of adjusted p-value (log scale)
  y = -log10(results_df$padj),
  # Plot as points (scatter plot)
  type = "p",
  # Use semi-transparent points to see overlapping points
  col = rgb(0, 0, 0, 0.3),
  # Point size
  pch = 16,
  # X-axis label
  xlab = "log2(Fold Change)",
  # Y-axis label
  ylab = "-log10(adjusted p-value)",
  # Plot title
  main = "mRNA in C. auris Cells: Volcano Plot\n(Treated vs Untreated)"
)

# Add threshold lines
# Vertical line at log2FC = ±1 (2-fold change)
abline(v = c(-1, 1), col = "blue", lty = 2)
# Horizontal line at padj = 0.05
abline(h = -log10(0.05), col = "red", lty = 2)

# Add legend
legend("topleft", 
       legend = c("log2FC = ±1", "padj = 0.05"),
       col = c("blue", "red"),
       lty = 2)

dev.off()
cat("✓ Saved Volcano plot to: 03_mRNA_Cells_Volcano_plot.pdf\n")

# PLOT 3: PCA plot of samples
# PCA (Principal Component Analysis) shows how samples cluster
# Samples from same condition should cluster together
pdf("/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline/03_mRNA_Cells_PCA_plot.pdf", width = 8, height = 6)

# Transform data for visualization
# log2 transformation helps stabilize variance
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

# Extract PCA results
# ntop = use top 500 most variable genes for PCA
pca <- prcomp(t(assay(vsd)), scale. = TRUE)

# Calculate proportion of variance explained by each PC
var_explained <- (pca$sdev^2) / sum(pca$sdev^2) * 100

# Create plot
plot(
  # PC1 scores for each sample
  x = pca$x[, 1],
  # PC2 scores for each sample
  y = pca$x[, 2],
  # Color by treatment
  col = ifelse(metadata_subset$Treatment == "Treated", "red", "blue"),
  # Point style
  pch = 16,
  # Point size
  cex = 2,
  # X-axis label with variance explained
  xlab = paste0("PC1 (", round(var_explained[1], 1), "%)"),
  # Y-axis label with variance explained
  ylab = paste0("PC2 (", round(var_explained[2], 1), "%)"),
  # Plot title
  main = "mRNA in C. auris Cells: PCA Plot"
)

# Add sample labels to points
text(pca$x[, 1], pca$x[, 2], labels = rownames(metadata_subset), 
     pos = 3, cex = 0.8)

# Add legend
legend("topright",
       legend = c("Untreated", "Treated"),
       col = c("blue", "red"),
       pch = 16)

dev.off()
cat("✓ Saved PCA plot to: 03_mRNA_Cells_PCA_plot.pdf\n")

# ============================================================================
# STEP 8: HEATMAP OF TOP DIFFERENTIALLY EXPRESSED mRNAs
# ============================================================================
# Heatmap visualization shows expression patterns across samples
# This helps identify sample clustering and expression trends

pdf("/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline/03_mRNA_Cells_Heatmap_Top20.pdf", width = 10, height = 8)

# Select top 20 significantly changed mRNAs
top_n <- 20
top_features <- results_df %>%
  # Remove NA values
  filter(!is.na(padj)) %>%
  # Get top 20 by adjusted p-value
  slice_head(n = top_n) %>%
  # Extract feature names
  pull(Feature)

# Extract normalized counts for top features
# Use variance-stabilizing transformation for better visualization
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
vsd_counts <- assay(vsd)[top_features, ]

# Create heatmap
pheatmap(
  # Count matrix for top features
  vsd_counts,
  # Cluster rows (features) based on expression patterns
  cluster_rows = TRUE,
  # Cluster columns (samples) based on expression patterns
  cluster_cols = TRUE,
  # Display feature names on left side
  show_rownames = TRUE,
  # Display sample names on bottom
  show_colnames = TRUE,
  # Color scale from blue (low) to red (high)
  color = colorRampPalette(c("blue", "white", "red"))(50),
  # Title for the heatmap
  main = paste("Top", top_n, "Differentially Expressed mRNAs\nin C. auris Cells"),
  # Font size for row names
  fontsize_row = 8,
  # Font size for column names
  fontsize_col = 10,
  # Scale expression values by row (feature) for better visualization
  scale = "row"
)

dev.off()
cat("✓ Saved heatmap to: 03_mRNA_Cells_Heatmap_Top20.pdf\n")

# ============================================================================
# STEP 9: SAVE COMPLETE RESULTS
# ============================================================================

# Save all results as CSV file (for viewing in Excel)
write.csv(results_df, 
          file = "/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline/03_mRNA_Cells_DESeq2_Results.csv",
          # Don't add row numbers as first column
          row.names = FALSE)
cat("✓ Saved all results to: 03_mRNA_Cells_DESeq2_Results.csv\n")

# Save just the significant features (padj < 0.05)
sig_results <- results_df %>%
  filter(padj < 0.05) %>%
  arrange(padj)

write.csv(sig_results,
          file = "/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline/03_mRNA_Cells_Significant_Features.csv",
          row.names = FALSE)
cat("✓ Saved significant features to: 03_mRNA_Cells_Significant_Features.csv\n")

# Save the DESeq2 object for later use
saveRDS(dds, file = "/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline/03_mRNA_Cells_DESeqDataSet.rds")
cat("✓ Saved DESeqDataSet object to: 03_mRNA_Cells_DESeqDataSet.rds\n")

# ============================================================================
# STEP 10: SUMMARY STATISTICS
# ============================================================================

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Analysis: mRNA in C. auris Cells\n")
cat("Comparison: Caspofungin-Treated vs Untreated\n")
cat("Total transcripts analyzed:", nrow(results_df), "\n")
cat("Significantly changed (padj < 0.05):", sig_005, "\n")
cat("  - Upregulated in treated:", sig_upregulated, "\n")
cat("  - Downregulated in treated:", sig_downregulated, "\n")
cat("\nOutput files created:\n")
cat("  1. 03_mRNA_Cells_DESeq2_Results.csv (all results)\n")
cat("  2. 03_mRNA_Cells_Significant_Features.csv (padj < 0.05)\n")
cat("  3. 03_mRNA_Cells_MA_plot.pdf\n")
cat("  4. 03_mRNA_Cells_Volcano_plot.pdf\n")
cat("  5. 03_mRNA_Cells_PCA_plot.pdf\n")
cat("  6. 03_mRNA_Cells_Heatmap_Top20.pdf\n")
cat("  7. 03_mRNA_Cells_DESeqDataSet.rds (R object)\n")

# ============================================================================
# EXTRACT SIGNIFICANTLY DE GENES FOR eggNOG MAPPER
# ============================================================================
# Get all significantly DE genes (not lncRNAs) and separate by regulation

cat("\n=== EXTRACTING DE GENES FOR eggNOG MAPPER ===\n")

# Filter for significant features (padj < 0.05)
sig_de_genes <- results_df %>%
  filter(padj < 0.05) %>%
  # Remove lncRNA features (keep only coding genes)
  filter(!grepl("^LNCRNA", Feature)) %>%
  # Sort by log2FoldChange
  arrange(desc(log2FoldChange))

cat("Total significantly DE genes (coding, padj < 0.05):", nrow(sig_de_genes), "\n")

# ============================================================================
# FILTER OUT NOVEL/UNANNOTATED GENES (MSTRG prefix)
# ============================================================================
# MSTRG genes are novel transcripts identified by StringTie but lack annotation
# eggNOG cannot annotate these, so we remove them for the analysis

cat("\n=== FILTERING OUT UNANNOTATED GENES (MSTRG) ===\n")

# Filter out MSTRG genes (keep only annotated genes)
sig_de_genes_annotated <- sig_de_genes %>%
  filter(!grepl("^MSTRG", Feature))

cat("DE genes before MSTRG removal:", nrow(sig_de_genes), "\n")
cat("DE genes after MSTRG removal:", nrow(sig_de_genes_annotated), "\n")
cat("MSTRG genes removed:", nrow(sig_de_genes) - nrow(sig_de_genes_annotated), "\n")

# Recalculate upregulated and downregulated with annotated genes only
upregulated_genes_annotated <- sig_de_genes_annotated %>%
  filter(log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange))

downregulated_genes_annotated <- sig_de_genes_annotated %>%
  filter(log2FoldChange < 0) %>%
  arrange(log2FoldChange)

# Create gene ID lists without MSTRG
all_de_gene_ids_annotated <- sig_de_genes_annotated %>% pull(Feature)
upregulated_gene_ids_annotated <- upregulated_genes_annotated %>% pull(Feature)
downregulated_gene_ids_annotated <- downregulated_genes_annotated %>% pull(Feature)

# Save updated lists for eggNOG (without MSTRG genes)
write.table(all_de_gene_ids_annotated,
            file = "/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline/Gene_IDs_All_DE_Annotated.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

write.table(upregulated_gene_ids_annotated,
            file = "/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline/Gene_IDs_Upregulated_Annotated.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

write.table(downregulated_gene_ids_annotated,
            file = "/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline/Gene_IDs_Downregulated_Annotated.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

cat("\n✓ Saved annotated gene lists (MSTRG removed) to:\n")
cat("  - Gene_IDs_All_DE_Annotated.txt\n")
cat("  - Gene_IDs_Upregulated_Annotated.txt\n")
cat("  - Gene_IDs_Downregulated_Annotated.txt\n")
# ============================================================================
# SEPARATE INTO UPREGULATED AND DOWNREGULATED
# ============================================================================

# Upregulated genes (log2FC > 0)
upregulated_genes <- sig_de_genes_annotated %>%
  filter(log2FoldChange > 0) %>%
  mutate(Feature = str_remove(Feature, "_transcript$")) %>%
  arrange(desc(log2FoldChange))
cat("Upregulated genes:", nrow(upregulated_genes), "\n")

# Downregulated genes (log2FC < 0)
downregulated_genes <- sig_de_genes %>%
  filter(log2FoldChange < 0) %>%
  mutate(Feature = str_remove(Feature, "_transcript$")) %>%
  # Sort by log2FoldChange (lowest first)
  arrange(log2FoldChange)
cat("Downregulated genes:", nrow(downregulated_genes), "\n")

# ============================================================================
# DISPLAY TOP GENES
# ============================================================================

cat("\n=== TOP 10 UPREGULATED GENES ===\n")
print(head(upregulated_genes, 10))

cat("\n=== TOP 10 DOWNREGULATED GENES ===\n")
print(head(downregulated_genes, 10))

# ============================================================================
# CREATE GENE ID LISTS FOR eggNOG MAPPER
# ============================================================================
# Extract just the gene IDs (feature names) for eggNOG input

# All significant genes
all_de_gene_ids <- sig_de_genes %>%
  pull(Feature)

# Upregulated genes only
upregulated_gene_ids <- upregulated_genes %>%
  pull(Feature)

# Downregulated genes only
downregulated_gene_ids <- downregulated_genes %>%
  pull(Feature)

cat("\n=== GENE LISTS CREATED ===\n")
cat("All DE genes:", length(all_de_gene_ids), "\n")
cat("Upregulated:", length(upregulated_gene_ids), "\n")
cat("Downregulated:", length(downregulated_gene_ids), "\n")

# ============================================================================
# SAVE AS CSV FILES FOR eggNOG MAPPER
# ============================================================================

# Save all DE genes with full statistics
write.csv(sig_de_genes,
          file = "/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline/DE_Genes_All_ForEggNOG.csv",
          row.names = FALSE)
cat("\n✓ Saved all DE genes to: DE_Genes_All_ForEggNOG.csv\n")

# Save upregulated genes
write.csv(upregulated_genes,
          file = "/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline/DE_Genes_Upregulated_ForEggNOG.csv",
          row.names = FALSE)
cat("✓ Saved upregulated genes to: DE_Genes_Upregulated_ForEggNOG.csv\n")

# Save downregulated genes
write.csv(downregulated_genes,
          file = "/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline/DE_Genes_Downregulated_ForEggNOG.csv",
          row.names = FALSE)
cat("✓ Saved downregulated genes to: DE_Genes_Downregulated_ForEggNOG.csv\n")

# ============================================================================
# SAVE GENE ID LISTS (MINIMAL FORMAT FOR eggNOG INPUT)
# ============================================================================
# These are just the gene IDs, one per line - ideal for eggNOG mapper input

write.table(all_de_gene_ids,
            file = "/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline/Gene_IDs_All_DE.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
cat("✓ Saved gene ID list (All) to: Gene_IDs_All_DE.txt\n")

write.table(upregulated_gene_ids,
            file = "/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline/Gene_IDs_Upregulated.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
cat("✓ Saved gene ID list (Upregulated) to: Gene_IDs_Upregulated.txt\n")

write.table(downregulated_gene_ids,
            file = "/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline/Gene_IDs_Downregulated.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
cat("✓ Saved gene ID list (Downregulated) to: Gene_IDs_Downregulated.txt\n")

cat("\n=== READY FOR eggNOG MAPPER ===\n")
cat("Use these files for eggNOG annotation:\n")
cat("  - Gene_IDs_All_DE.txt (all significant genes)\n")
cat("  - Gene_IDs_Upregulated.txt (upregulated only)\n")
cat("  - Gene_IDs_Downregulated.txt (downregulated only)\n")

