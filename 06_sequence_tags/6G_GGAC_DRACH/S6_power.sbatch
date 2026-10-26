#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=s6g_S6
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6G_GGAC_DRACH/logs/%x_%j.out
# Step 6G DRACH power check, stage 1: dose-series sets (g8) + task list for S7.
set -eo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6G_GGAC_DRACH; cd $B; source inputs/paths.sh
md5sum -c --quiet inputs/inputs.md5; (cd scripts && md5sum -c --quiet SCRIPTS_r3.md5)
for d in G8_power_sets G8_tasks.tsv; do [ ! -e run/$d ] || { echo "STOP: run/$d exists"; exit 1; }; done
$PY scripts/g8_power_sets6g.py $F3 $SEQ run/G8_power_sets 20261015
awk -F'\t' 'NR>1{print "run/G8_power_sets\t"$1}' run/G8_power_sets/sets_summary.tsv > run/G8_tasks.tsv; wc -l run/G8_tasks.tsv
(cd run/G8_power_sets && find . -type f ! -name '*.md5' -exec md5sum {} + > G8_outputs.md5); chmod -R a-w run/G8_power_sets run/G8_tasks.tsv; echo S6_DONE
