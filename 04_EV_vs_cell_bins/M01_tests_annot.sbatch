#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=bin100_M01
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.out
#SBATCH --error=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.err
# M0: known-answer test of binomial_test100.R; M1: merged-final-bin annotation (D1) from the frozen 100-nt annotation.
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3; R=$B/run
[ ! -e $R/M0_binomial_test ] && [ ! -e $R/M1_annotation_merged ] || { echo "STOP: exists"; exit 1; }
source /opt/exp_soft/miniconda3-py3.12/etc/profile.d/conda.sh; set +u; conda activate /home/flljus002/.conda/envs/drimseq_env; set -u
$PY $B/scripts/test_binomial100.py $B/scripts $R/M0_binomial_test $(command -v Rscript) | tee $R/M0_binomial_test_result.txt
grep -q ALL_BINOMIAL_TESTS_PASS $R/M0_binomial_test_result.txt
(cd $R/B1_annotation/outputs && md5sum -c FROZEN_BIN100_ANNOT.md5 --quiet)
mkdir -p $R/M1_annotation_merged
$PY $B/scripts/make_merged_annotation100m.py $R/B1_annotation/outputs/bin100_master_annotation.tsv $R/M1_annotation_merged | tee $R/M1_log.txt
(cd $R/M1_annotation_merged && md5sum * > FROZEN_BIN100M_ANNOT.md5 && chmod a-w *)
echo M01_DONE_OK
