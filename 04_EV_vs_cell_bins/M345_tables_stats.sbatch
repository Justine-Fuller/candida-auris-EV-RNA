#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=bin100_M345
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.out
#SBATCH --error=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.err
# M3 gates/tables, M4 tie rule, merge consistency vs the unmerged counts, M5 one-sided binomial test (D1-D8) + independent
# Python re-check -- for EV and Cell sRNA separately. Refuses to overwrite.
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3; R=$B/run; ANN=$R/M1_annotation_merged/bin100m_master_annotation.tsv; C=$R/M2_counts_merged
for d in M3_tables_EV M3_tables_cellsRNA M4_tie_rule_EV M4_tie_rule_cellsRNA M5_binomial_EV M5_binomial_cellsRNA; do [ ! -e $R/$d ] || { echo "STOP: $d exists"; exit 1; }; done
(cd $C && md5sum *.gz *.json > M2_outputs.md5 && chmod a-w *)
source /opt/exp_soft/miniconda3-py3.12/etc/profile.d/conda.sh; set +u; conda activate /home/flljus002/.conda/envs/drimseq_env; set -u
for X in "EV EV manifest_EV.tsv B2_counts" "cellsRNA sRNA manifest_Cell_sRNA.tsv C2_counts_cellsRNA"; do
  set -- $X; T=$1; K=$2; MAN=$B/$3; OLD=$R/$4
  $PY $B/scripts/build_tables100.py $ANN $C $MAN $R/M3_tables_$T $K | tee $R/M3_${T}_log.txt
  $PY $B/scripts/apply_tie_rule100.py $ANN $C $MAN $R/M4_tie_rule_$T | tee $R/M4_${T}_log.txt
  $PY $B/scripts/compare_merge100.py $OLD $C $MAN $ANN $R/M4_tie_rule_$T/merge_consistency.tsv | tee $R/M4_${T}_merge_log.txt
  Rscript $B/scripts/binomial_test100.R $R/M4_tie_rule_$T $ANN $MAN $R/M5_binomial_$T 20261008 | tee $R/M5_${T}_log.txt
  $PY $B/scripts/check_binomial100.py $R/M4_tie_rule_$T $ANN $MAN $R/M5_binomial_$T 20261008 | tee $R/M5_${T}_check_log.txt
  for d in M3_tables_$T M4_tie_rule_$T M5_binomial_$T; do (cd $R/$d && md5sum * > ../${d}.md5 && mv ../${d}.md5 ${d}_outputs.md5 && chmod a-w *); done
done
echo M345_DONE_OK
