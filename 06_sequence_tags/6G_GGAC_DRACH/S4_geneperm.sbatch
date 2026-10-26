#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=s6g_S4
#SBATCH --array=1-22%11
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6G_GGAC_DRACH/logs/%x_%A_%a.out
# Step 6G gene level: one set per task (run/G6_tasks.tsv has NO header). Usage: sbatch S4_geneperm.sbatch <GC|noGC>
set -eo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6G_GGAC_DRACH; cd $B
COV=$1; [ "$COV" = GC ] || [ "$COV" = noGC ] || { echo "STOP: COV"; exit 1; }
(cd scripts && md5sum -c --quiet SCRIPTS_r1.md5)
read -r F4 SET < <(awk -F'\t' -v i=$SLURM_ARRAY_TASK_ID 'NR==i' run/G6_tasks.tsv); [ -n "$SET" ] || { echo "STOP: empty row"; exit 1; }
source /opt/exp_soft/miniconda3-py3.12/etc/profile.d/conda.sh; conda activate /home/flljus002/.conda/envs/drimseq_env
Rscript scripts/g6_gene_perm6g.R $F4 $SET run/G6_perm_$COV $((20261014 + SLURM_ARRAY_TASK_ID)) 100 100000 10 $COV
chmod a-w run/G6_perm_$COV/$SET.*; echo S4_TASK_DONE
