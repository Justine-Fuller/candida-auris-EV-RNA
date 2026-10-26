#!/bin/bash
#SBATCH -p ada -c 1 --mem=8G -t 04:00:00 -J c6r3_S5 --array=1-48 -o /scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6C_codon/logs/%x_%A_%a.out
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6C_codon; export R_LIBS_USER=/home/flljus002/Rlibs
cd $B/run_r3; mkdir -p S5_null
Rscript $B/scripts_r3/run_null_chunk.R S4_null_sets/null_foregrounds.tsv S3_enrichment_real S5_null $SLURM_ARRAY_TASK_ID 48 $B/scripts_r3
