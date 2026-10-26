#!/bin/bash
#SBATCH --job-name=predgpi_wp
#SBATCH --partition=ada
#SBATCH --output=/scratch/flljus002/CandidaEV_Project/06_PCG_functional_annotation/09_whole_proteome_2026-10-14/06_PredGPI/predgpi_wp_%j.log
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=06:00:00
# Same settings as run_predgpi_up/down.sh, input = whole proteome.
set -eo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate predgpi
export PREDGPI_HOME="/scratch/flljus002/software/predgpi"
OUT=/scratch/flljus002/CandidaEV_Project/06_PCG_functional_annotation/09_whole_proteome_2026-10-14/06_PredGPI/results
if [ -e "$OUT/wp_predgpi.json" ]; then echo "output exists, refusing to overwrite"; exit 1; fi
mkdir -p "$OUT"
cd /scratch/flljus002/software/predgpi
./predgpi.py -f /scratch/flljus002/CandidaEV_Project/06_PCG_functional_annotation/09_whole_proteome_2026-10-14/input/proteome_5417.fasta -o "$OUT/wp_predgpi.json" -m json
chmod a-w "$OUT/wp_predgpi.json"
