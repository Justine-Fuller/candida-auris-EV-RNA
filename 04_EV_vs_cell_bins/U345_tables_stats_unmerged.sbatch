#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=bin100_U345
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.out
#SBATCH --error=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.err
# U3-U5 SENSITIVITY (unmerged bins, "half" rule): tables, tie rule, consistency vs unmerged B2/C2 counts, binomial + check.
# for EV and Cell sRNA separately, after the 3'-overhang exclusion. Refuses to overwrite.
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3; R=$B/run; ANN=$R/B1_annotation/outputs/bin100_master_annotation.tsv; C=$R/U2_counts_unmerged_3prime
for d in U3_tables_unmerged_EV U3_tables_unmerged_cellsRNA U4_tie_rule_unmerged_EV U4_tie_rule_unmerged_cellsRNA U5_binomial_unmerged_EV U5_binomial_unmerged_cellsRNA; do [ ! -e $R/$d ] || { echo "STOP: $d exists"; exit 1; }; done
(cd $C && md5sum *.gz *.json > U2_outputs.md5 && chmod a-w *)
source /opt/exp_soft/miniconda3-py3.12/etc/profile.d/conda.sh; set +u; conda activate /home/flljus002/.conda/envs/drimseq_env; set -u
for X in "EV EV manifest_EV.tsv B2_counts" "cellsRNA sRNA manifest_Cell_sRNA.tsv C2_counts_cellsRNA"; do
  set -- $X; T=$1; K=$2; MAN=$B/$3; OLD=$R/$4
  $PY $B/scripts/build_tables100.py $ANN $C $MAN $R/U3_tables_unmerged_$T $K | tee $R/U3_${T}_log.txt
  $PY $B/scripts/apply_tie_rule100.py $ANN $C $MAN $R/U4_tie_rule_unmerged_$T | tee $R/U4_${T}_log.txt
  $PY $B/scripts/compare_3prime100.py $OLD $C $MAN $ANN $R/U4_tie_rule_unmerged_$T/three_prime_consistency.tsv | tee $R/U4_${T}_3prime_log.txt
  Rscript $B/scripts/binomial_test100.R $R/U4_tie_rule_unmerged_$T $ANN $MAN $R/U5_binomial_unmerged_$T 20261008 | tee $R/U5_${T}_log.txt
  $PY $B/scripts/check_binomial100.py $R/U4_tie_rule_unmerged_$T $ANN $MAN $R/U5_binomial_unmerged_$T 20261008 | tee $R/U5_${T}_check_log.txt
  for d in U3_tables_unmerged_$T U4_tie_rule_unmerged_$T U5_binomial_unmerged_$T; do (cd $R/$d && md5sum * > ../${d}.md5 && mv ../${d}.md5 ${d}_outputs.md5 && chmod a-w *); done
done
echo U345_DONE_OK
