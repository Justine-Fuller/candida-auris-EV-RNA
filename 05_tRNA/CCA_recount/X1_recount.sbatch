#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=cca2_X1
#SBATCH --array=1-36%8
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_tRNA_2026-10-09/cca_recount_2026-10-13/logs/%x_%A_%a.out
# CCA recount per library (EV, Cell sRNA, Cell polyA). manifest_all36.tsv HAS a header.
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_tRNA_2026-10-09; C=$B/cca_recount_2026-10-13; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3
read -r LIB BAM KIND _ < <(awk -F'\t' -v i=$((SLURM_ARRAY_TASK_ID+1)) 'NR==i' $B/manifest_all36.tsv)
[ -n "$LIB" ] || exit 1; [ ! -e $C/run/R1_cca2/$LIB.cca2.json ] || { echo STOP exists; exit 1; }
K=sRNA; [ "$KIND" = polyA ] && K=polyA
(cd $C/scripts && md5sum -c --quiet SCRIPTS_r1.md5); (cd $B/inputs && md5sum -c --quiet inputs.md5)
mkdir -p $C/run/R1_cca2
$PY $C/scripts/trna_cca_recount.py $BAM $LIB $B/inputs/candidate_trna_reference_168.bed \
  /scratch/flljus002/CandidaEV_Project/_Archive_2026-08-31/01_Differential_expression_analysis/01_Reference_files/C_auris_B8441_version_s02-m01-r01_chromosomes.fasta $C/run/R1_cca2 $K
chmod a-w $C/run/R1_cca2/$LIB.cca2.json
