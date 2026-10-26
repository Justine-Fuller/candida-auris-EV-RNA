#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=s6g_S3
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6G_GGAC_DRACH/logs/%x_%j.out
# Step 6G gene level, sets: synthetic (planted GGACT / GC shift) + real FX sets (Step 6F F3 gene tests + F1 sequences); task list.
set -eo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6G_GGAC_DRACH; cd $B; source inputs/paths.sh
md5sum -c --quiet inputs/inputs.md5; (cd scripts && md5sum -c --quiet SCRIPTS_r1.md5)
for d in G5_sim G5_sets G6_tasks.tsv; do [ ! -e run/$d ] || { echo "STOP: run/$d exists"; exit 1; }; done
$PY scripts/g5_sets6g.py sim $F3 $SEQ run/G5_sim 20261014
$PY scripts/g5_sets6g.py real $F3 $SEQ run/G5_sets
( awk -F'\t' 'NR>1{print "run/G5_sim\t"$1}' run/G5_sim/sets_summary.tsv; awk -F'\t' 'NR>1{print "run/G5_sets\t"$1}' run/G5_sets/sets_summary.tsv ) > run/G6_tasks.tsv
wc -l run/G6_tasks.tsv
for d in G5_sim G5_sets; do (cd run/$d && find . -type f ! -name '*.md5' -exec md5sum {} + > ${d}_outputs.md5); chmod -R a-w run/$d; done
chmod a-w run/G6_tasks.tsv; echo S3_DONE
