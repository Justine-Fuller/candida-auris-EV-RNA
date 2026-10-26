# ============================================================================
# DESeq2 ANALYSIS SETUP AND DATA PREPARATION
# Candida auris Caspofungin Treatment Response
# ============================================================================
# This script prepares the count matrix and metadata for downstream DESeq2 analyses
# It filters transcripts at the feature level and removes lowly expressed features

# ============================================================================
# STEP 1: INSTALL AND LOAD REQUIRED PACKAGES
# ============================================================================

# Install BiocManager (only needed once)
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# Install DESeq2 from Bioconductor (only needed once)
# BiocManager::install("DESeq2")

# Load required libraries
library(DESeq2)      # Main package for differential expression analysis
library(tidyverse)   # Data manipulation (includes dplyr, ggplot2, etc.)
library(pheatmap)    # For creating heatmaps of expression patterns

# Display R and package versions for reproducibility
cat("R version:", R.version$version.string, "\n")
cat("DESeq2 version:", packageVersion("DESeq2"), "\n")

# ============================================================================
# STEP 2: LOAD COUNT MATRIX AND METADATA
# ============================================================================

# Read the BEDTools feature count matrix
# This is a tab-separated file where:
#   - Row 1: Feature names (transcript/exon/UTR/lncRNA annotations)
#   - Columns 2-N: Sample IDs with raw read counts
# Read WITHOUT setting row.names (to avoid duplicate error)
count_matrix_raw <- read.delim(
  "/Users/justinefuller/Downloads/DESeq2_pipeline/BEDTools_feature_matrix.tsv",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

# Extract feature names (first column)
feature_names <- count_matrix_raw[, 1]

# Get count data (all other columns)
count_data <- count_matrix_raw[, -1]

# Create unique rownames by appending a counter to duplicates
unique_rownames <- make.unique(as.character(feature_names), sep = "_")

# Create dataframe with proper rownames built in
count_matrix <- as.data.frame(count_data)
rownames(count_matrix) <- unique_rownames

cat("\n=== COUNT MATRIX LOADED ===\n")
cat("Total features:", nrow(count_matrix), "\n")
cat("Samples:", ncol(count_matrix), "\n")
cat("Sample of features:\n")
print(head(rownames(count_matrix), 10))

# Read the sample metadata
# This contains information about each sample:
#   - Sample: Sample ID
#   - Compartment: Cell or EV (extracellular vesicle)
#   - Treatment: Treated or Untreated (caspofungin treatment)
#   - Strain: B8441 or MMC1 (fungal strain)
#   - RNA_type: sRNA or mRNA
sample_metadata <- read.delim(
  "/Users/justinefuller/Downloads/DESeq2_pipeline/sample_metadata.tsv",
  # Set first column as rownames
  row.names = 1,
  # Specify that file is tab-separated
  sep = "\t",
  # Treat first row as header
  header = TRUE,
  # Avoid converting strings to factors
  stringsAsFactors = FALSE
)

# Display basic information about loaded data
cat("\n=== SAMPLE METADATA INFORMATION ===\n")
# Show dimensions: number of samples (rows) and attributes (columns)
cat("Dimensions (Samples x Attributes):", dim(sample_metadata), "\n")
# Show first few rows of metadata to verify structure
print(head(sample_metadata))

# ============================================================================
# STEP 3: EXTRACT DIFFERENT FEATURE TYPES FOR SEPARATE ANALYSES
# ============================================================================
# Extract different feature types for separate DESeq2 analyses:
#   1. Transcript-level features (_transcript)
#   2. Non-coding features (exons, introns, 5'UTR, 3'UTR)

# Extract feature names from rownames
feature_names <- rownames(count_matrix)

cat("\n=== FEATURE EXTRACTION FOR SEPARATE ANALYSES ===\n")
cat("Total features:", nrow(count_matrix), "\n")

# ============================================================================
# EXTRACT 1: TRANSCRIPT-LEVEL FEATURES ONLY
# ============================================================================
# Keep only features that end with "_transcript"
# These represent complete transcripts (including lncRNAs)

is_transcript <- grepl("_(transcript)(_[0-9]+)?$", rownames(count_matrix))

# Create count matrix with transcript features only
count_matrix_transcripts <- count_matrix[is_transcript, ]

cat("\n=== TRANSCRIPT-LEVEL FEATURES ===\n")
cat("Transcript features:", nrow(count_matrix_transcripts), "\n")
cat("Samples:", ncol(count_matrix_transcripts), "\n")
cat("Example transcripts:\n")
print(head(rownames(count_matrix_transcripts), 5))

# ============================================================================
# EXTRACT 2: NON-CODING FEATURES (EXONS, INTRONS, UTRs)
# ============================================================================
# Keep features that are exons, introns, 5'UTR, or 3'UTR
# Do NOT include _transcript features (those are in the separate analysis)

is_ncfeature <- grepl("_(exon|intron|5_prime_UTR|3_prime_UTR)(_[0-9]+)?$", rownames(count_matrix))

# Create count matrix with non-coding features
count_matrix_ncfeatures <- count_matrix[is_ncfeature, ]

cat("\n=== NON-CODING FEATURES (exons, introns, UTRs) ===\n")
cat("Non-coding features:", nrow(count_matrix_ncfeatures), "\n")
cat("Samples:", ncol(count_matrix_ncfeatures), "\n")
cat("Breakdown by feature type:\n")

# Count each feature type
exons <- sum(grepl("_exon$", rownames(count_matrix_ncfeatures)))
introns <- sum(grepl("_intron$", rownames(count_matrix_ncfeatures)))
utr5 <- sum(grepl("_5_prime_UTR$", rownames(count_matrix_ncfeatures)))
utr3 <- sum(grepl("_3_prime_UTR$", rownames(count_matrix_ncfeatures)))

cat("  Exons:", exons, "\n")
cat("  Introns:", introns, "\n")
cat("  5' UTRs:", utr5, "\n")
cat("  3' UTRs:", utr3, "\n")

cat("Example non-coding features:\n")
print(head(rownames(count_matrix_ncfeatures), 5))

# ============================================================================
# VERIFICATION: CHECK THAT ALL FEATURES ARE ACCOUNTED FOR
# ============================================================================

total_extracted <- nrow(count_matrix_transcripts) + nrow(count_matrix_ncfeatures)

cat("\n=== FEATURE EXTRACTION SUMMARY ===\n")
cat("Total features:", nrow(count_matrix), "\n")
cat("Transcript features extracted:", nrow(count_matrix_transcripts), "\n")
cat("Non-coding features extracted:", nrow(count_matrix_ncfeatures), "\n")
cat("Total features extracted:", total_extracted, "\n")

if (total_extracted == nrow(count_matrix)) {
  cat("✓ All features successfully extracted and categorized\n")
} else {
  cat("WARNING: Mismatch in feature counts\n")
}

# ============================================================================
# STEP 4: VERIFY SAMPLE MATCHING AND DATA INTEGRITY
# ============================================================================

# Check that all samples in count matrices are in metadata
# Both matrices should have the same samples
samples_in_transcripts <- colnames(count_matrix_transcripts)
samples_in_ncfeatures <- colnames(count_matrix_ncfeatures)
samples_in_metadata <- rownames(sample_metadata)

# Verify both matrices have the same samples
if (identical(samples_in_transcripts, samples_in_ncfeatures)) {
  cat("\n✓ Both count matrices have identical sample sets\n")
  samples_in_counts <- samples_in_transcripts
} else {
  cat("\nWARNING: Sample mismatch between transcript and non-coding feature matrices\n")
  samples_in_counts <- samples_in_transcripts  # Use transcripts as primary
}

# Find any mismatches with metadata
missing_from_metadata <- setdiff(samples_in_counts, samples_in_metadata)
extra_in_metadata <- setdiff(samples_in_metadata, samples_in_counts)

if (length(missing_from_metadata) > 0) {
  cat("\nWARNING: Samples in counts but not in metadata:\n")
  print(missing_from_metadata)
}

if (length(extra_in_metadata) > 0) {
  cat("\nWARNING: Samples in metadata but not in counts:\n")
  print(extra_in_metadata)
}

if (length(missing_from_metadata) == 0 & length(extra_in_metadata) == 0) {
  cat("\n✓ All samples match between count matrices and metadata\n")
}

# Ensure metadata rows are in same order as count matrix columns
sample_metadata <- sample_metadata[samples_in_counts, ]

cat("\n=== FINAL DATA SUMMARY ===\n")
cat("Transcript features:", nrow(count_matrix_transcripts), "\n")
cat("Non-coding features:", nrow(count_matrix_ncfeatures), "\n")
cat("Samples:", ncol(count_matrix_transcripts), "\n")
cat("Sample breakdown by RNA type:\n")
print(table(sample_metadata$RNA_type))
cat("\nSample breakdown by Compartment:\n")
print(table(sample_metadata$Compartment))
cat("\nSample breakdown by Treatment:\n")
print(table(sample_metadata$Treatment))

# ============================================================================
# STEP 5: REMOVE LOWLY EXPRESSED FEATURES
# ============================================================================
# This step is important for downstream analysis because:
#   - Features with very low counts have high variance
#   - Removing them reduces multiple testing burden
#   - It improves statistical power for detecting true signals

# Define minimum count threshold
# We'll keep features that have at least 10 reads in at least 3 samples
min_count <- 10
min_samples <- 3

cat("\n=== LOW EXPRESSION FILTERING ===\n")
cat("Threshold: At least", min_count, "reads in at least", min_samples, "samples\n")

# ============================================================================
# FILTER TRANSCRIPT FEATURES
# ============================================================================

# Calculate the number of samples where each feature has sufficient counts
feature_present_transcripts <- rowSums(count_matrix_transcripts >= min_count) >= min_samples

cat("\nTranscript Features:\n")
cat("  Features before filtering:", nrow(count_matrix_transcripts), "\n")
cat("  Features passing filter:", sum(feature_present_transcripts), "\n")
cat("  Features removed:", sum(!feature_present_transcripts), "\n")

# Apply the filter to transcript count matrix
count_matrix_transcripts <- count_matrix_transcripts[feature_present_transcripts, ]

cat("  Final dimensions:", dim(count_matrix_transcripts), "\n")

# ============================================================================
# FILTER NON-CODING FEATURES
# ============================================================================

# Calculate the number of samples where each feature has sufficient counts
feature_present_ncfeatures <- rowSums(count_matrix_ncfeatures >= min_count) >= min_samples

cat("\nNon-coding Features:\n")
cat("  Features before filtering:", nrow(count_matrix_ncfeatures), "\n")
cat("  Features passing filter:", sum(feature_present_ncfeatures), "\n")
cat("  Features removed:", sum(!feature_present_ncfeatures), "\n")

# Apply the filter to non-coding count matrix
count_matrix_ncfeatures <- count_matrix_ncfeatures[feature_present_ncfeatures, ]

cat("  Final dimensions:", dim(count_matrix_ncfeatures), "\n")

# ============================================================================
# STEP 6: EXPORT PROCESSED DATA FOR DOWNSTREAM ANALYSES
# ============================================================================

# Create output directory if it doesn't exist
dir.create("/Users/justinefuller/Downloads/DESeq2_pipeline", showWarnings = FALSE)

cat("\n=== SAVING PROCESSED DATA ===\n")

# ============================================================================
# SAVE TRANSCRIPT-LEVEL COUNTS
# ============================================================================

saveRDS(count_matrix_transcripts, file = "/Users/justinefuller/Downloads/DESeq2_pipeline/count_matrix_transcripts.rds")
cat("✓ Saved transcript count matrix to: count_matrix_transcripts.rds\n")

write.csv(count_matrix_transcripts, 
          file = "/Users/justinefuller/Downloads/DESeq2_pipeline/count_matrix_transcripts.csv")
cat("✓ Saved transcript count matrix (CSV) to: count_matrix_transcripts.csv\n")

# ============================================================================
# SAVE NON-CODING FEATURE COUNTS
# ============================================================================

saveRDS(count_matrix_ncfeatures, file = "/Users/justinefuller/Downloads/DESeq2_pipeline/count_matrix_ncfeatures.rds")
cat("✓ Saved non-coding feature count matrix to: count_matrix_ncfeatures.rds\n")

write.csv(count_matrix_ncfeatures, 
          file = "/Users/justinefuller/Downloads/DESeq2_pipeline/count_matrix_ncfeatures.csv")
cat("✓ Saved non-coding feature count matrix (CSV) to: count_matrix_ncfeatures.csv\n")

# ============================================================================
# SAVE METADATA
# ============================================================================

saveRDS(sample_metadata, file = "/Users/justinefuller/Downloads/DESeq2_pipeline/sample_metadata.rds")
cat("✓ Saved sample metadata to: sample_metadata.rds\n")

write.csv(sample_metadata, 
          file = "/Users/justinefuller/Downloads/DESeq2_pipeline/sample_metadata.csv")
cat("✓ Saved sample metadata (CSV) to: sample_metadata.csv\n")

cat("\n=== SETUP COMPLETE ===\n")
cat("Ready to proceed with DESeq2 analyses\n")
cat("\nLoad these files in analysis scripts with:\n")
cat("  # For transcript-level analysis:\n")
cat("  count_matrix <- readRDS('/Users/justinefuller/Downloads/DESeq2_pipeline/count_matrix_transcripts.rds')\n")
cat("  sample_metadata <- readRDS('/Users/justinefuller/Downloads/DESeq2_pipeline/sample_metadata.rds')\n")
cat("\n  # For non-coding feature analysis:\n")
cat("  count_matrix <- readRDS('/Users/justinefuller/Downloads/DESeq2_pipeline/count_matrix_ncfeatures.rds')\n")
cat("  sample_metadata <- readRDS('/Users/justinefuller/Downloads/DESeq2_pipeline/sample_metadata.rds')\n")
