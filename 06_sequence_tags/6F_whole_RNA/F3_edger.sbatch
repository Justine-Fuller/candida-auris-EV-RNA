#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=s6f_F3
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6F_whole_RNA_tag/logs/%x_%j.out
# Step 6F F3: gene-level EV vs Cell sRNA edgeR per group; enriched / comparison-pool counts per region subset. STOP after this.
set -eo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6F_whole_RNA_tag
source $B/inputs/paths.sh
md5sum -c --quiet $B/inputs/inputs.md5; (cd $B/scripts && md5sum -c --quiet SCRIPTS_r1.md5)
[ $(grep -l $'F2_GATE\tPASS' $B/run/F2_counts/*.gate.tsv | wc -l) -eq 24 ] || { echo "STOP: F2 gates not all PASS"; exit 1; }
source /opt/exp_soft/miniconda3-py3.12/etc/profile.d/conda.sh; conda activate /home/flljus002/.conda/envs/drimseq_env
D=$B/run/F3_gene_tests
Rscript $B/scripts/s6f_edger_genes.R $MAN6E $B/run/F2_counts $B/run/F1_sequences/gene_sequences.tsv.gz $TPMTAB $MAN36 $D
(cd $D && md5sum * > F3_outputs.md5); chmod -R a-w $D
