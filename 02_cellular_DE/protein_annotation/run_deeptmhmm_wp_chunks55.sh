#!/bin/bash
#SBATCH --job-name=deeptmhmm_wp
#SBATCH --partition=ada
#SBATCH --array=1-55%4
#SBATCH --output=/scratch/flljus002/CandidaEV_Project/06_PCG_functional_annotation/09_whole_proteome_2026-10-14/05_DeepTMHMM/deeptmhmm_wp_chunk55_%a_%A.log
#SBATCH --error=/scratch/flljus002/CandidaEV_Project/06_PCG_functional_annotation/09_whole_proteome_2026-10-14/05_DeepTMHMM/deeptmhmm_wp_chunk55_%a_%A.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
# Same settings as run_deeptmhmm_up/down.sh (DeepTMHMM 1.0.24 via biolib), input = whole proteome in 55 chunks of <=99 proteins
# (the full 2.85 MB file was rejected by biolib as too large; job 1477378).
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate deeptmhmm
OUT=/scratch/flljus002/CandidaEV_Project/06_PCG_functional_annotation/09_whole_proteome_2026-10-14/05_DeepTMHMM/results/chunk${SLURM_ARRAY_TASK_ID}
if [ -e "$OUT/biolib_results" ]; then echo "output exists, refusing to overwrite"; exit 1; fi
mkdir -p "$OUT"; cd "$OUT"
biolib run DTU/DeepTMHMM:1.0.24 --fasta /scratch/flljus002/CandidaEV_Project/06_PCG_functional_annotation/09_whole_proteome_2026-10-14/input/chunks55/proteome_chunk${SLURM_ARRAY_TASK_ID}.fasta
chmod -R a-w "$OUT"/biolib_results
