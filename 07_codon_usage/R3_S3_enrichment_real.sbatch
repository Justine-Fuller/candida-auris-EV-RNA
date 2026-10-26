#!/bin/bash
#SBATCH -p ada -c 2 --mem=16G -t 02:00:00 -J c6r3_S3 -o /scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6C_codon/logs/%x_%j.out
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6C_codon; V1=/scratch/flljus002/CandidaEV_Project/05_EV_coverage_analysis/11_Full_length_biological_annotation/output/codon_usage_analysis
export R_LIBS_USER=/home/flljus002/Rlibs
cd $B/run_r3; [ -e S3_enrichment_real ] && { echo exists; exit 1; }; mkdir S3_enrichment_real
# same KEGG snapshot as r1/r2; r2 significant terms added to the tracked list (comparison only)
Rscript $B/scripts_r3/run_enrichment_real.R S2_real_scoring S3_enrichment_real $V1/human_enrichment/output $B/scripts_r3 $B/run/S3_enrichment_real/gson_KEGG_hsa.rds $B/run_r2/S3_enrichment_real
