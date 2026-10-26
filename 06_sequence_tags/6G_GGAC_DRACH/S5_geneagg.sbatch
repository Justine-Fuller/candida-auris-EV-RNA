#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=s6g_S5
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6G_GGAC_DRACH/logs/%x_%j.out
# Step 6G gene level aggregation. Usage: sbatch S5_geneagg.sbatch <GC|noGC>
set -eo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6G_GGAC_DRACH; cd $B
COV=$1; [ "$COV" = GC ] || [ "$COV" = noGC ] || { echo "STOP: COV"; exit 1; }
(cd scripts && md5sum -c --quiet SCRIPTS_r1.md5)
source /opt/exp_soft/miniconda3-py3.12/etc/profile.d/conda.sh; conda activate /home/flljus002/.conda/envs/drimseq_env
set +e; Rscript scripts/g7_gene_aggregate6g.R run/G5_sim run/G5_sets run/G6_perm_$COV run/G7_genes_$COV $COV; rc=$?; set -e
(cd run/G6_perm_$COV && find . -type f ! -name '*.md5' -exec md5sum {} + > G6_outputs.md5); (cd run/G7_genes_$COV && md5sum * > G7_outputs.md5)
chmod -R a-w run/G6_perm_$COV run/G7_genes_$COV; exit $rc
