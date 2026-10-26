#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=s6h_T1
#SBATCH --array=1-24%8
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6H_DRACH_position/logs/%x_%A_%a.out
# Step 6H T1: 6E PCG_CDS reads with position class, per library (manifest_6E.tsv HAS a header).
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6H_DRACH_position; cd $B; source inputs/paths.sh
read -r LIB COMP GROUP BAM < <(awk -F'\t' -v i=$((SLURM_ARRAY_TASK_ID+1)) 'NR==i' $MAN6E); [ -n "$LIB" ] || { echo "STOP: empty row"; exit 1; }
mkdir -p run/H1_positions; [ ! -e run/H1_positions/$LIB.h1.json.gz ] || { echo "STOP: exists"; exit 1; }
grep "  $BAM\$" $P/v2_protein_coding/provenance/J00/$LIB.md5 | md5sum -c -
md5sum -c --quiet inputs/inputs.md5; (cd scripts && md5sum -c --quiet SCRIPTS_r1.md5)
$PY scripts/h1_positions6h.py $BAM $LIB $A100 $N100 $AM run/H1_positions > run/H1_positions/$LIB.accounting.json
chmod a-w run/H1_positions/$LIB.*; echo T1_DONE
