#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=trna_T1
#SBATCH --array=1-24%8
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_tRNA_2026-10-09/logs/%x_%A_%a.out
#SBATCH --error=/scratch/flljus002/Phase3D_forensic_audit/v2_tRNA_2026-10-09/logs/%x_%A_%a.err
# tRNA step T1: family counts + positions + PCG/lncRNA bin coverage per library (EV + Cell sRNA). manifest.tsv HAS a header.
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_tRNA_2026-10-09; P=/scratch/flljus002/Phase3D_forensic_audit; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3
read -r LIB BAM KIND COMP _ < <(awk -F'\t' -v i=$((SLURM_ARRAY_TASK_ID+1)) 'NR==i' $B/manifest.tsv)
[ -n "$LIB" ] || { echo "STOP: empty row"; exit 1; }
[ ! -e $B/run/T1_counts/$LIB.trna.json ] || { echo "STOP: exists"; exit 1; }
grep "  $BAM\$" $P/v2_protein_coding/provenance/J00/$LIB.md5 | md5sum -c -
(cd $B/inputs && md5sum -c --quiet inputs.md5)
mkdir -p $B/run/T1_counts
$PY $B/scripts/trna_count.py $BAM $LIB $B/inputs/bin100_master_annotation.tsv $B/inputs/candidate_trna_reference_168.bed $B/run/T1_counts
echo T1_LIB_DONE
