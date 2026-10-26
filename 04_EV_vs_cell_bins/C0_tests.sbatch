#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=bin100_C0_tests
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.out
#SBATCH --error=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.err
# C0: known-answer test r3 of the r2 counters (EV answers unchanged; kind sRNA identical to EV; polyA refused).
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1; [ ! -e $B/run/C0_counter_test_r3 ] || { echo "STOP: exists"; exit 1; }
/opt/exp_soft/miniconda3-py3.12/bin/python3 $B/scripts/test_counters100.py $B/scripts $B/run/C0_counter_test_r3 | tee $B/run/C0_counter_test_r3_result.txt
grep -q ALL_COUNTER100_TESTS_PASS $B/run/C0_counter_test_r3_result.txt
echo C0_DONE_OK
