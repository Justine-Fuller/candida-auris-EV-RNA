#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=s6f_F6
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6F_whole_RNA_tag/logs/%x_%j.out
# Step 6F F6: STREME (MEME 5.5.7, defaults except --rna and --seed), enriched vs comparison sequences, FX sets (primary) + the
# synthetic planted set (does STREME find the planted word?). STREME cannot use covariates: abundance / length / GC are NOT adjusted.
set -eo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6F_whole_RNA_tag
cd $B; (cd scripts && md5sum -c --quiet SCRIPTS_r7.md5); [ ! -e run/F6_streme ] || { echo "STOP: exists"; exit 1; }; mkdir run/F6_streme
source /opt/exp_soft/miniconda3-py3.12/etc/profile.d/conda.sh; conda activate /home/flljus002/.conda/envs/meme_env
echo -e "set\tenriched_seqs\tcomparison_seqs\tstatus\tmotif_line" > run/F6_streme/summary.tsv
for d in run/F4_sets/sets/FX.* run/F4_sim/sets/FX.*.planted; do
  s=$(basename $d); ne=$(grep -c '>' $d/enriched.fa || true); nc=$(grep -c '>' $d/comparison.fa || true)
  if [ "$ne" -lt 1 ] || [ "$nc" -lt 1 ]; then echo -e "$s\t$ne\t$nc\tskipped_empty_side\t" >> run/F6_streme/summary.tsv; continue; fi
  if streme --p $d/enriched.fa --n $d/comparison.fa --rna --seed 20261013 --oc run/F6_streme/$s > run/F6_streme/$s.log 2>&1; then
    n=$(grep -c '^MOTIF' run/F6_streme/$s/streme.txt || true)
    if [ "$n" -eq 0 ]; then echo -e "$s\t$ne\t$nc\tok_0_motifs\t" >> run/F6_streme/summary.tsv
    else grep '^MOTIF' run/F6_streme/$s/streme.txt | while read -r l; do echo -e "$s\t$ne\t$nc\tok\t$l" >> run/F6_streme/summary.tsv; done; fi
  else echo -e "$s\t$ne\t$nc\tstreme_error\t$(tail -1 run/F6_streme/$s.log)" >> run/F6_streme/summary.tsv; fi
done
cat run/F6_streme/summary.tsv
(cd run/F6_streme && find . -type f ! -name '*.md5' -exec md5sum {} + > F6_outputs.md5); chmod -R a-w run/F6_streme; echo F6_DONE
