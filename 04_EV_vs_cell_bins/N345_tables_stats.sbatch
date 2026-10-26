#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=bin100_N345
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.out
#SBATCH --error=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.err
# N3 gates/tables, N4 tie rule, consistency vs M2 (only excluded reads differ), N5 binomial (D1-D8) + independent check,
# for EV and Cell sRNA separately, after the 3'-overhang exclusion. Refuses to overwrite.
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3; R=$B/run; ANN=$R/M1_annotation_merged/bin100m_master_annotation.tsv; C=$R/N2_counts_3prime
for d in N3_tables_EV N3_tables_cellsRNA N4_tie_rule_EV N4_tie_rule_cellsRNA N5_binomial_EV N5_binomial_cellsRNA; do [ ! -e $R/$d ] || { echo "STOP: $d exists"; exit 1; }; done
(cd $C && md5sum *.gz *.json > N2_outputs.md5 && chmod a-w *)
source /opt/exp_soft/miniconda3-py3.12/etc/profile.d/conda.sh; set +u; conda activate /home/flljus002/.conda/envs/drimseq_env; set -u
for X in "EV EV manifest_EV.tsv" "cellsRNA sRNA manifest_Cell_sRNA.tsv"; do
  set -- $X; T=$1; K=$2; MAN=$B/$3
  $PY $B/scripts/build_tables100.py $ANN $C $MAN $R/N3_tables_$T $K | tee $R/N3_${T}_log.txt
  $PY $B/scripts/apply_tie_rule100.py $ANN $C $MAN $R/N4_tie_rule_$T | tee $R/N4_${T}_log.txt
  $PY $B/scripts/compare_3prime100.py $R/M2_counts_merged $C $MAN $ANN $R/N4_tie_rule_$T/three_prime_consistency.tsv | tee $R/N4_${T}_3prime_log.txt
  Rscript $B/scripts/binomial_test100.R $R/N4_tie_rule_$T $ANN $MAN $R/N5_binomial_$T 20261008 | tee $R/N5_${T}_log.txt
  $PY $B/scripts/check_binomial100.py $R/N4_tie_rule_$T $ANN $MAN $R/N5_binomial_$T 20261008 | tee $R/N5_${T}_check_log.txt
  for d in N3_tables_$T N4_tie_rule_$T N5_binomial_$T; do (cd $R/$d && md5sum * > ../${d}.md5 && mv ../${d}.md5 ${d}_outputs.md5 && chmod a-w *); done
done
echo N345_DONE_OK
