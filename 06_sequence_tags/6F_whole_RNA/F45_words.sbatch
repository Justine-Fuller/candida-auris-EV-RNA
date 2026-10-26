#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=s6f_F45
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6F_whole_RNA_tag/logs/%x_%j.out
# Step 6F F4 + F5 (Justine "V2"): self-test, synthetic gates (planted word, GC shift), real sets, word regression + null.
# Stops at the first failure (synthetic gates before the real run).
set -eo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6F_whole_RNA_tag; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3
cd $B; md5sum -c --quiet inputs/inputs.md5; (cd scripts && md5sum -c --quiet SCRIPTS_r5.md5)
for d in F4_sim F5_sim F4_sets F5_words; do [ ! -e run/$d ] || { echo "STOP: run/$d exists"; exit 1; }; done
$PY scripts/s6f_sets_kmers.py selftest | tee run/F4_selftest.txt
source /opt/exp_soft/miniconda3-py3.12/etc/profile.d/conda.sh; conda activate /home/flljus002/.conda/envs/drimseq_env
$PY scripts/s6f_sets_kmers.py sim run/F3_gene_tests run/F1_sequences/gene_sequences.tsv.gz run/F4_sim 20261013
Rscript scripts/s6f_word_regression.R run/F4_sim run/F5_sim sim 20261013 0
$PY scripts/s6f_sets_kmers.py real run/F3_gene_tests run/F1_sequences/gene_sequences.tsv.gz run/F4_sets
Rscript scripts/s6f_word_regression.R run/F4_sets run/F5_words real 20261013 100
for d in F4_sim F5_sim F4_sets F5_words; do (cd run/$d && find . -type f ! -name '*.md5' -exec md5sum {} + > ${d}_outputs.md5); chmod -R a-w run/$d; done
chmod a-w run/F4_selftest.txt; echo F45_DONE
