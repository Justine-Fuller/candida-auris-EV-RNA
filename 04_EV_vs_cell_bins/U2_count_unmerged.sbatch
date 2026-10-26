#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=bin100_U2_count
#SBATCH --array=1-24%8
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%A_%a.out
#SBATCH --error=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%A_%a.err
# U2: SENSITIVITY recount, UNMERGED 100-nt bins (last bin at true length), registered "half" 3-prime overhang rule (2026-10-09).
# manifest_all24.tsv HAS a header: row = task + 1. Refuses to overwrite.
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1; R=$B/run/U2_counts_unmerged_3prime; ANN=$B/run/B1_annotation/outputs/bin100_master_annotation.tsv
read -r LIB BAM KIND COMP _ < <(awk -F'\t' -v i=$((SLURM_ARRAY_TASK_ID+1)) 'NR==i' $B/manifest_all24.tsv)
[ -n "$LIB" ] && { [ "$KIND" = EV ] || [ "$KIND" = sRNA ]; } || { echo "STOP: bad row"; exit 1; }
[ ! -e $R/$LIB.accounting.json ] || { echo "STOP: exists"; exit 1; }
(cd $B/run/B1_annotation/outputs && md5sum -c FROZEN_BIN100_ANNOT.md5 --quiet) || { echo "STOP: annotation changed"; exit 1; }
grep "  $BAM\$" /scratch/flljus002/Phase3D_forensic_audit/v2_protein_coding/provenance/J00/$LIB.md5 | md5sum -c -
mkdir -p $R; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3
$PY $B/scripts/count_bins100.py $ANN $BAM $LIB $KIND $R
$PY $B/scripts/count_bins100_c1.py $ANN $BAM $LIB $KIND $R
echo U2_DONE_OK
