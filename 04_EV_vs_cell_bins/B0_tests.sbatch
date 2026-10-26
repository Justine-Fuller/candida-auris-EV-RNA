#!/bin/bash
#SBATCH --partition=ada
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.out
#SBATCH --error=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.err
#SBATCH --job-name=bin100_B0_tests
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
# B0: known-answer tests of both 100-nt counters (C0, C1). Refuses to overwrite.
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1
[ ! -e $B/run/B0_counter_test ] || { echo "STOP: exists"; exit 1; }
/opt/exp_soft/miniconda3-py3.12/bin/python3 $B/scripts/test_counters100.py $B/scripts $B/run/B0_counter_test | tee $B/run/B0_counter_test_result.txt
echo B0_DONE_OK
