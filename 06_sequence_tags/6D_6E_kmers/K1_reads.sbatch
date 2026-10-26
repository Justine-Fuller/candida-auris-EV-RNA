#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=s6e_K1
#SBATCH --array=1-24%8
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6DE/logs/%x_%A_%a.out
# Step 6E K1: RNA-oriented read sequences per RNA type per library (12 EV + 12 Cell sRNA). manifest_6E.tsv HAS a header.
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6DE; P=/scratch/flljus002/Phase3D_forensic_audit; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3
read -r LIB COMP GROUP BAM < <(awk -F'\t' -v i=$((SLURM_ARRAY_TASK_ID+1)) 'NR==i' $B/inputs/manifest_6E.tsv)
[ -n "$LIB" ] || { echo "STOP: empty row"; exit 1; }
[ ! -e $B/run/K1_reads/$LIB.k1.json.gz ] || { echo "STOP: exists"; exit 1; }
grep "  $BAM\$" $P/v2_protein_coding/provenance/J00/$LIB.md5 | md5sum -c -
md5sum -c --quiet $B/inputs/inputs.md5
(cd $B/scripts && md5sum -c --quiet SCRIPTS_r1.md5)
mkdir -p $B/run/K1_reads
$PY $B/scripts/kmer_reads6e.py $BAM $LIB $P/v2_tRNA_2026-10-09/inputs/bin100_master_annotation.tsv $P/v2_rRNA_2026-10-09/FL3_annotation/new_bins100.tsv \
    $P/v2_tRNA_2026-10-09/inputs/candidate_trna_reference_168.bed $P/v2_rRNA_2026-10-09/FL3_annotation/rRNA_5S_loci_V3.bed $B/run/K1_reads > $B/run/K1_reads/$LIB.accounting.json
chmod a-w $B/run/K1_reads/$LIB.k1.json.gz $B/run/K1_reads/$LIB.accounting.json
echo K1_LIB_DONE
