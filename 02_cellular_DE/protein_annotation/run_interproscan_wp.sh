#!/bin/bash
#SBATCH --job-name=ips_wp
#SBATCH --partition=ada
#SBATCH --output=/scratch/flljus002/CandidaEV_Project/06_PCG_functional_annotation/09_whole_proteome_2026-10-14/02_InterProScan/ips_wp_%j.log
#SBATCH --error=/scratch/flljus002/CandidaEV_Project/06_PCG_functional_annotation/09_whole_proteome_2026-10-14/02_InterProScan/ips_wp_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
# Same settings as run_interproscan_up/down.sh (InterProScan 5.78-109.0), input = whole proteome.
set -euo pipefail
IPS_DIR="/scratch/flljus002/software/interproscan/interproscan-5.78-109.0"
INPUT="/scratch/flljus002/CandidaEV_Project/06_PCG_functional_annotation/09_whole_proteome_2026-10-14/input/proteome_5417.fasta"
OUTPUT_DIR="/scratch/flljus002/CandidaEV_Project/06_PCG_functional_annotation/09_whole_proteome_2026-10-14/02_InterProScan/results"
TEMP_DIR="/scratch/flljus002/interproscan_tmp/wp_${SLURM_JOB_ID}"
if [ -e "$OUTPUT_DIR/wp_interproscan.tsv" ]; then echo "output exists, refusing to overwrite"; exit 1; fi
mkdir -p "$OUTPUT_DIR" "$TEMP_DIR"
cd "$IPS_DIR"
./interproscan.sh -i "$INPUT" -b "$OUTPUT_DIR/wp_interproscan" -f tsv,json,gff3 -goterms -pa -dp -dra -cpu 8 -T "$TEMP_DIR"
rm -rf "$TEMP_DIR"
chmod a-w "$OUTPUT_DIR"/wp_interproscan.*
