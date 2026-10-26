#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=s6f_F5agg
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6F_whole_RNA_tag/logs/%x_%j.out
# Step 6F F5 r7 aggregation (usage: sbatch F5_aggregate.sbatch <GC|noGC>): gates (incl. synthetic + random-label null), significant words, >= 3/4 candidates. Freezes F5_perm.
set -eo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6F_whole_RNA_tag
COV=$1; [ "$COV" = GC ] || [ "$COV" = noGC ] || { echo "STOP: COV"; exit 1; }
cd $B; (cd scripts && md5sum -c --quiet SCRIPTS_r7.md5)
source /opt/exp_soft/miniconda3-py3.12/etc/profile.d/conda.sh; conda activate /home/flljus002/.conda/envs/drimseq_env
set +e; Rscript scripts/s6f_word_aggregate.R run/F4_sim run/F4_sets run/F5_perm_$COV run/F5_words_$COV $COV; rc=$?; set -e
(cd run/F5_perm_$COV && find . -type f ! -name '*.md5' -exec md5sum {} + > F5_perm_outputs.md5); (cd run/F5_words_$COV && md5sum * > F5_words_outputs.md5)
chmod -R a-w run/F5_perm_$COV run/F5_words_$COV; exit $rc
