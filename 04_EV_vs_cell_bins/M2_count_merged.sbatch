#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=bin100_M2_count
#SBATCH --array=1-24%8
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%A_%a.out
#SBATCH --error=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%A_%a.err
# M2: recount 12 EV + 12 Cell sRNA libraries on the merged-final-bin annotation (D1), C0 + independent C1 (r2 counters).
# manifest_all24.tsv HAS a header: row = task + 1. Refuses to overwrite.
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1; R=$B/run/M2_counts_merged; ANN=$B/run/M1_annotation_merged/bin100m_master_annotation.tsv
read -r LIB BAM KIND COMP _ < <(awk -F'\t' -v i=$((SLURM_ARRAY_TASK_ID+1)) 'NR==i' $B/manifest_all24.tsv)
[ -n "$LIB" ] && { [ "$KIND" = EV ] || [ "$KIND" = sRNA ]; } || { echo "STOP: bad row"; exit 1; }
[ ! -e $R/$LIB.accounting.json ] || { echo "STOP: exists"; exit 1; }
(cd $B/run/M1_annotation_merged && md5sum -c FROZEN_BIN100M_ANNOT.md5 --quiet) || { echo "STOP: annotation changed"; exit 1; }
grep "  $BAM\$" /scratch/flljus002/Phase3D_forensic_audit/v2_protein_coding/provenance/J00/$LIB.md5 | md5sum -c -
mkdir -p $R; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3
$PY $B/scripts/count_bins100.py $ANN $BAM $LIB $KIND $R
$PY $B/scripts/count_bins100_c1.py $ANN $BAM $LIB $KIND $R
echo M2_DONE_OK
