#!/bin/bash
#SBATCH --partition=ada
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.out
#SBATCH --error=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.err
#SBATCH --job-name=bin100_B3_tables
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
# B3: validation gates + per-replicate tables (no statistics). Refuses to overwrite.
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1; O=$B/run/B3_tables
[ ! -e $O ] || { echo "STOP: exists"; exit 1; }
/opt/exp_soft/miniconda3-py3.12/bin/python3 $B/scripts/build_tables100.py $B/run/B1_annotation/outputs/bin100_master_annotation.tsv $B/run/B2_counts $B/manifest_EV.tsv $O | tee $B/run/B3_log.txt
(cd $O && md5sum * > B3_outputs.md5)
echo B3_DONE_OK
