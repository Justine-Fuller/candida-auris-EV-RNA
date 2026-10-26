#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=bin100_U0_tests
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.out
#SBATCH --error=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.err
# U0: known-answer test of binomial_test100.R r2 (layout-agnostic; final_bin flag).
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3; O=$B/run/U0_tests
[ ! -e $O ] || { echo "STOP: exists"; exit 1; }; mkdir -p $O
source /opt/exp_soft/miniconda3-py3.12/etc/profile.d/conda.sh; set +u; conda activate /home/flljus002/.conda/envs/drimseq_env; set -u
$PY $B/scripts/test_binomial100.py $B/scripts $O/binomial $(command -v Rscript) | tee $O/binomial_result.txt; grep -q ALL_BINOMIAL_TESTS_PASS $O/binomial_result.txt
echo U0_DONE_OK
