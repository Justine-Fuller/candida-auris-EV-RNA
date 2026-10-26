#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=bin100_S0_tests
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.out
#SBATCH --error=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.err
# S0: r4 counters (optional 'any' overhang rule) -- original known-answer test and the 3' test incl. the 'any' rule.
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3; O=$B/run/S0_tests
[ ! -e $O ] || { echo "STOP: exists"; exit 1; }; mkdir -p $O
$PY $B/scripts/test_counters100.py $B/scripts $O/counters | tee $O/counters_result.txt; grep -q ALL_COUNTER100_TESTS_PASS $O/counters_result.txt
$PY $B/scripts/test_counters100_3prime.py $B/scripts $O/three_prime | tee $O/three_prime_result.txt; grep -q ALL_3PRIME_TESTS_PASS $O/three_prime_result.txt
echo S0_DONE_OK
