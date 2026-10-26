#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=bin100_N0_tests
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.out
#SBATCH --error=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.err
# N0: r3 counters (3'-overhang exclusion) -- original known-answer test (unchanged expectations), new 3' test, tie-rule test
# with r3 table scripts, binomial test. All must pass.
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3; R=$B/run; O=$R/N0_tests
[ ! -e $O ] || { echo "STOP: exists"; exit 1; }; mkdir -p $O
$PY $B/scripts/test_counters100.py $B/scripts $O/counters | tee $O/counters_result.txt; grep -q ALL_COUNTER100_TESTS_PASS $O/counters_result.txt
$PY $B/scripts/test_counters100_3prime.py $B/scripts $O/three_prime | tee $O/three_prime_result.txt; grep -q ALL_3PRIME_TESTS_PASS $O/three_prime_result.txt
$PY $B/scripts/test_tie_rule100.py $B/scripts $O/tie_rule | tee $O/tie_rule_result.txt; grep -q ALL_TIE_RULE_TESTS_PASS $O/tie_rule_result.txt
source /opt/exp_soft/miniconda3-py3.12/etc/profile.d/conda.sh; set +u; conda activate /home/flljus002/.conda/envs/drimseq_env; set -u
$PY $B/scripts/test_binomial100.py $B/scripts $O/binomial $(command -v Rscript) | tee $O/binomial_result.txt; grep -q ALL_BINOMIAL_TESTS_PASS $O/binomial_result.txt
echo N0_DONE_OK
