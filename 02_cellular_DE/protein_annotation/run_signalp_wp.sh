#!/bin/bash
#SBATCH --job-name=signalp_wp
#SBATCH --partition=ada
#SBATCH --output=/scratch/flljus002/CandidaEV_Project/06_PCG_functional_annotation/09_whole_proteome_2026-10-14/04_SignalP/signalp_wp_%j.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=08:00:00
# Same settings as run_signalp_up/down.sh (SignalP 6.0), input = whole proteome.
set -eo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate signalp6
OUT=/scratch/flljus002/CandidaEV_Project/06_PCG_functional_annotation/09_whole_proteome_2026-10-14/04_SignalP/results
if [ -e "$OUT/prediction_results.txt" ]; then echo "output exists, refusing to overwrite"; exit 1; fi
mkdir -p "$OUT"
signalp6 --fastafile /scratch/flljus002/CandidaEV_Project/06_PCG_functional_annotation/09_whole_proteome_2026-10-14/input/proteome_5417.fasta --organism eukarya --output_dir "$OUT" --format none --mode fast --torch_num_threads 8 --write_procs 4
chmod a-w "$OUT"/*
