#!/bin/bash
#SBATCH --job-name=deeptmhmm_retry
#SBATCH --partition=ada
#SBATCH --output=/scratch/flljus002/CandidaEV_Project/06_PCG_functional_annotation/09_whole_proteome_2026-10-14/05_DeepTMHMM/retry_chunk%a_%A.log
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=06:00:00
# Retry of the 37 chunks that biolib refused anonymously in job 1477570 ("compute_limit_exceeded_please_signin")
# or that timed out. Same DeepTMHMM 1.0.24 call; up to 8 attempts, 15 min apart.
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate deeptmhmm
OUT=/scratch/flljus002/CandidaEV_Project/06_PCG_functional_annotation/09_whole_proteome_2026-10-14/05_DeepTMHMM/results/chunk${SLURM_ARRAY_TASK_ID}
if [ -e "$OUT/biolib_results" ]; then echo "output exists, refusing to overwrite"; exit 1; fi
mkdir -p "$OUT"
for attempt in 1 2 3 4 5 6 7 8; do
  echo "attempt $attempt $(date)"
  T=$(mktemp -d "$OUT/../.tmp_chunk${SLURM_ARRAY_TASK_ID}_XXXX"); cd "$T"
  timeout 40m biolib run DTU/DeepTMHMM:1.0.24 --fasta /scratch/flljus002/CandidaEV_Project/06_PCG_functional_annotation/09_whole_proteome_2026-10-14/input/chunks55/proteome_chunk${SLURM_ARRAY_TASK_ID}.fasta
  if [ -s "$T/biolib_results/TMRs.gff3" ] && [ -s "$T/biolib_results/predicted_topologies.3line" ]; then
    mv "$T/biolib_results" "$OUT/"; cd "$OUT"; rmdir "$T"; chmod -R a-w "$OUT"/biolib_results; echo DONE; exit 0
  fi
  cd "$OUT"; rm -rf "$T"; sleep 900
done
echo "GAVE UP"; exit 1
