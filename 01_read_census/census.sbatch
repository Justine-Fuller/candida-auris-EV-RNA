#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=census_count
#SBATCH --array=1-24%8
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_read_census_2026-10-09/logs/%x_%A_%a.out
#SBATCH --error=/scratch/flljus002/Phase3D_forensic_audit/v2_read_census_2026-10-09/logs/%x_%A_%a.err
# READ-ONLY read census per library (EV + Cell sRNA, frozen v2 T1b BAMs). manifest.tsv HAS a header: row = task + 1.
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_read_census_2026-10-09; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3; P=/scratch/flljus002/Phase3D_forensic_audit
read -r LIB BAM KIND COMP _ < <(awk -F'\t' -v i=$((SLURM_ARRAY_TASK_ID+1)) 'NR==i' $B/manifest.tsv)
[ -n "$LIB" ] || { echo "STOP: empty row"; exit 1; }
[ ! -e $B/run/counts/$LIB.census.json ] || { echo "STOP: exists"; exit 1; }
grep "  $BAM\$" $P/v2_protein_coding/provenance/J00/$LIB.md5 | md5sum -c -
(cd $B/inputs && md5sum -c --quiet inputs.md5)
mkdir -p $B/run/counts
$PY $B/scripts/census.py $BAM $LIB $B/inputs/intervals.tsv $B/inputs/candidate_trna_reference_168.bed $B/run/counts
echo CENSUS_LIB_DONE
