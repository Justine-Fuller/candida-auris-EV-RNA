#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=bin100_S345
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.out
#SBATCH --error=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.err
# S3-S5 SENSITIVITY ("any" 3-prime overhang rule): tables, tie rule, consistency vs M2, binomial + independent check.
# for EV and Cell sRNA separately, after the 3'-overhang exclusion. Refuses to overwrite.
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3; R=$B/run; ANN=$R/M1_annotation_merged/bin100m_master_annotation.tsv; C=$R/S2_counts_anyoverhang
for d in S3_tables_anyoverhang_EV S3_tables_anyoverhang_cellsRNA S4_tie_rule_anyoverhang_EV S4_tie_rule_anyoverhang_cellsRNA S5_binomial_anyoverhang_EV S5_binomial_anyoverhang_cellsRNA; do [ ! -e $R/$d ] || { echo "STOP: $d exists"; exit 1; }; done
(cd $C && md5sum *.gz *.json > S2_outputs.md5 && chmod a-w *)
source /opt/exp_soft/miniconda3-py3.12/etc/profile.d/conda.sh; set +u; conda activate /home/flljus002/.conda/envs/drimseq_env; set -u
for X in "EV EV manifest_EV.tsv" "cellsRNA sRNA manifest_Cell_sRNA.tsv"; do
  set -- $X; T=$1; K=$2; MAN=$B/$3
  $PY $B/scripts/build_tables100.py $ANN $C $MAN $R/S3_tables_anyoverhang_$T $K | tee $R/S3_${T}_log.txt
  $PY $B/scripts/apply_tie_rule100.py $ANN $C $MAN $R/S4_tie_rule_anyoverhang_$T | tee $R/S4_${T}_log.txt
  $PY $B/scripts/compare_3prime100.py $R/M2_counts_merged $C $MAN $ANN $R/S4_tie_rule_anyoverhang_$T/three_prime_consistency.tsv | tee $R/S4_${T}_3prime_log.txt
  Rscript $B/scripts/binomial_test100.R $R/S4_tie_rule_anyoverhang_$T $ANN $MAN $R/S5_binomial_anyoverhang_$T 20261008 | tee $R/S5_${T}_log.txt
  $PY $B/scripts/check_binomial100.py $R/S4_tie_rule_anyoverhang_$T $ANN $MAN $R/S5_binomial_anyoverhang_$T 20261008 | tee $R/S5_${T}_check_log.txt
  for d in S3_tables_anyoverhang_$T S4_tie_rule_anyoverhang_$T S5_binomial_anyoverhang_$T; do (cd $R/$d && md5sum * > ../${d}.md5 && mv ../${d}.md5 ${d}_outputs.md5 && chmod a-w *); done
done
echo S345_DONE_OK
