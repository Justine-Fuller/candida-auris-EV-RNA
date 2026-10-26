#!/bin/bash
#SBATCH --job-name=lc_ds_caur
#SBATCH --output=/scratch/flljus002/Candida_Host_Splicing_GSE154911/logs/lc_ds_caur_%j.log
#SBATCH --error=/scratch/flljus002/Candida_Host_Splicing_GSE154911/logs/lc_ds_caur_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00

set -eo pipefail
source ~/.bashrc
conda activate leafcutter_env
export R_LIBS_USER=/home/flljus002/Rlibs

PREFIX=/scratch/flljus002/Candida_Host_Splicing_GSE154911/03_LeafCutter/ds_results/RPMI_vs_caur_live

Rscript /scratch/flljus002/Candida_Host_Splicing_GSE154911/environment_setup/leafcutter/scripts/leafcutter_ds.R \
  -o $PREFIX \
  -e /scratch/flljus002/Candida_Host_Splicing_GSE154911/01_Reference/gtf/gencode.v50.exons.txt.gz \
  -p 4 \
  /scratch/flljus002/Candida_Host_Splicing_GSE154911/03_LeafCutter/clustering/gse154911_perind_numers.counts.gz \
  /scratch/flljus002/Candida_Host_Splicing_GSE154911/03_LeafCutter/groups/RPMI_vs_caur_live.txt

for f in ${PREFIX}_cluster_significance.txt ${PREFIX}_effect_sizes.txt; do
  if [ ! -s "$f" ]; then
    echo "ERROR: $f is missing or empty" >&2
    exit 1
  fi
done

echo "caur DS test complete ($(wc -l < ${PREFIX}_cluster_significance.txt) cluster rows)"
