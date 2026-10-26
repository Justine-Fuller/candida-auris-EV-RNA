#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=s6f_F7
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6F_whole_RNA_tag/logs/%x_%j.out
# Step 6F F7: direct GC test (EV-enriched vs comparison genes, abundance + length covariates, permutation p), per group x region.
set -eo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6F_whole_RNA_tag
cd $B; (cd scripts && md5sum -c --quiet SCRIPTS_r7.md5)
source /opt/exp_soft/miniconda3-py3.12/etc/profile.d/conda.sh; conda activate /home/flljus002/.conda/envs/drimseq_env
set +e; Rscript scripts/s6f_gc_test.R run/F4_sets run/F4_sim run/F7_gc 20261014 100 100000 10; rc=$?; set -e
(cd run/F7_gc && md5sum * > F7_outputs.md5); chmod -R a-w run/F7_gc; exit $rc
