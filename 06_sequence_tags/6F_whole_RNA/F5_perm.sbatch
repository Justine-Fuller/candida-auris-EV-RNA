#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=s6f_F5
#SBATCH --array=1-42%14
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6F_whole_RNA_tag/logs/%x_%A_%a.out
# Step 6F F5 r7: usage sbatch F5_perm.sbatch <GC|noGC>. One set per task (run/F5_tasks.tsv has NO header: <F4 dir> <set>). Permutation word test + random-label null.
set -eo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6F_whole_RNA_tag
COV=$1; [ "$COV" = GC ] || [ "$COV" = noGC ] || { echo "STOP: COV"; exit 1; }
cd $B; (cd scripts && md5sum -c --quiet SCRIPTS_r7.md5)
read -r F4 SET < <(awk -F'\t' -v i=$SLURM_ARRAY_TASK_ID 'NR==i' run/F5_tasks.tsv); [ -n "$SET" ] || { echo "STOP: empty row"; exit 1; }
source /opt/exp_soft/miniconda3-py3.12/etc/profile.d/conda.sh; conda activate /home/flljus002/.conda/envs/drimseq_env
Rscript scripts/s6f_word_perm.R $F4 $SET run/F5_perm_$COV $((20261014 + SLURM_ARRAY_TASK_ID)) 100 100000 10 $COV
chmod a-w run/F5_perm_$COV/$SET.*; echo F5_TASK_DONE
