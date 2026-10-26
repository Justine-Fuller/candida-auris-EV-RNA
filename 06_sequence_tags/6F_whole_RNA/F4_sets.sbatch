#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=s6f_F4
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6F_whole_RNA_tag/logs/%x_%j.out
# Step 6F F4 (r6): counter self-test, synthetic sets (planted / GC shift), real sets (FX + FE); task list for the F5 array.
set -eo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6F_whole_RNA_tag; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3
cd $B; md5sum -c --quiet inputs/inputs.md5; (cd scripts && md5sum -c --quiet SCRIPTS_r7.md5)
for d in F4_sim F4_sets F4_selftest.txt F5_tasks.tsv; do [ ! -e run/$d ] || { echo "STOP: run/$d exists"; exit 1; }; done
$PY scripts/s6f_sets_kmers.py selftest > run/F4_selftest.txt; cat run/F4_selftest.txt
$PY scripts/s6f_sets_kmers.py sim run/F3_gene_tests run/F1_sequences/gene_sequences.tsv.gz run/F4_sim 20261013
$PY scripts/s6f_sets_kmers.py real run/F3_gene_tests run/F1_sequences/gene_sequences.tsv.gz run/F4_sets
( awk -F'\t' 'NR>1{print "run/F4_sim\t"$1}' run/F4_sim/sets_summary.tsv; awk -F'\t' 'NR>1{print "run/F4_sets\t"$1}' run/F4_sets/sets_summary.tsv ) > run/F5_tasks.tsv
wc -l run/F5_tasks.tsv
for d in F4_sim F4_sets; do (cd run/$d && find . -type f ! -name '*.md5' -exec md5sum {} + > ${d}_outputs.md5); chmod -R a-w run/$d; done
chmod a-w run/F4_selftest.txt run/F5_tasks.tsv; echo F4_DONE
