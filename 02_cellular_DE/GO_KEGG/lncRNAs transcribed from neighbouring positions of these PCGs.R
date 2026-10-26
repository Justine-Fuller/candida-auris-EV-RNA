#lncRNAs transcribed from neighbouring positions of these PCGs
#Run only if you don't already have these installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
a
BiocManager::install(c(
  "rtracklayer",
  "GenomicRanges"
))
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(readr)

# ============================================================
# 1. IMPORT YOUR GTF
# ============================================================

gtf_file <- "/Users/justinefuller/Library/Mobile Documents/com~apple~CloudDocs/University documents/Honours/Caro/CandidaEV/DESeq2_pipeline/final.gtf"

gtf <- import(gtf_file)

gtf

#inspect what annotation fields your GTF actually contains:
mcols(gtf) |> names()
head(as.data.frame(gtf))

#extract only the transcript rows, then split them into lncRNAs and PCGs.
# Keep transcript-level entries only
transcripts <- gtf[gtf$type == "transcript"]

# Check how many transcripts there are
length(transcripts)

# Inspect some transcript IDs
head(transcripts$transcript_id, 20)

#Then split them by ID pattern:
# lncRNAs
lncrna_transcripts <- transcripts[
  grepl("^LNCRNA", transcripts$gene_id)
]

# Protein-coding genes
pcg_transcripts <- transcripts[
  grepl("^B9J08_", transcripts$gene_id)
]

#check how many of each I have:
length(lncrna_transcripts)
length(pcg_transcripts)

#inspect them:
head(as.data.frame(lncrna_transcripts))
head(as.data.frame(pcg_transcripts))

# Find nearest protein-coding gene for each lncRNA
nearest_hits <- distanceToNearest(
  lncrna_transcripts,
  pcg_transcripts,
  ignore.strand = TRUE
)

# Build a table of lncRNA-PCG neighbour pairs
lncrna_pcg_pairs <- data.frame(
  lncRNA = lncrna_transcripts$gene_id[
    queryHits(nearest_hits)
  ],
  PCG = pcg_transcripts$gene_id[
    subjectHits(nearest_hits)
  ],
  distance_bp = mcols(nearest_hits)$distance
)

# Inspect the first few rows
head(lncrna_pcg_pairs)

#check the distance distribution:
summary(lncrna_pcg_pairs$distance_bp)

table(
  cut(
    lncrna_pcg_pairs$distance_bp,
    breaks = c(-1, 0, 1000, 5000, 10000, Inf),
    labels = c(
      "Overlapping",
      "<=1 kb",
      "1-5 kb",
      "5-10 kb",
      ">10 kb"
    )
  )
)

#distanceToNearest() reports: distance = 0 when the lncRNA and PCG 
# genomic ranges overlap. It is therefore useful to keep those 
# separate from merely nearby genes.