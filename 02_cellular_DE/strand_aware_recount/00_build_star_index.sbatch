#!/bin/bash
#SBATCH --job-name=cellmRNA_star_index
#SBATCH --partition=ada
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/flljus002/CandidaEV_Project/08_Cell_mRNA_corrected_2026-09-25/logs/%x_%j.out
# Builds the STAR index from 04_STAR_and_counting_PROPOSAL.md (VALIDATION_2026-09-25_categoryA).
# Skips if the index already exists. Writes only inside 08_Cell_mRNA_corrected_2026-09-25/.
set -euo pipefail
B=/scratch/flljus002/CandidaEV_Project/08_Cell_mRNA_corrected_2026-09-25
REF=/scratch/flljus002/CandidaEV_Project/_Archive_2026-08-31/01_Differential_expression_analysis/01_Reference_files
IDX=$B/02_STAR/index
mkdir -p $B/logs
if [ -s $IDX/SA ]; then echo "Index already exists at $IDX, nothing to do"; exit 0; fi
# Same reference files as the provenance log (md5 checked so nothing silently changes)
echo "e07db201a43d5d2ec0db905b27f02e10  $REF/C_auris_B8441_version_s02-m01-r01_chromosomes.fasta" | md5sum -c -
echo "6ef80f048cebceaa8b7b0fd1b33bcd13  $REF/final.gtf" | md5sum -c -
mkdir -p $IDX
module load software/STAR-2.7.11b
STAR --runMode genomeGenerate --runThreadN 6 \
  --genomeDir $IDX \
  --genomeFastaFiles $REF/C_auris_B8441_version_s02-m01-r01_chromosomes.fasta \
  --sjdbGTFfile $REF/final.gtf \
  --sjdbOverhang 49 --genomeSAindexNbases 10 \
  --outFileNamePrefix $B/02_STAR/index_build_
echo INDEX_DONE_OK
