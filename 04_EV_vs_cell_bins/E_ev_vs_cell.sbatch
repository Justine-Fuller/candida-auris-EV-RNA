#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=bin100_E_evcell
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.out
#SBATCH --error=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.err
# E: step 5/6 EV-vs-Cell sRNA bin-share comparison (edgeR diffSplice.DGEGLM, QL) + combination with hotspot calls.
# E0 known-answer test; E1 main (merged layout, N4 tables + N5 hotspots); E2 sensitivity (unmerged, U4 + U5);
# each with an informational null-simulation run and an independent Python check. Refuses to overwrite.
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1; R=$B/run; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3
for d in E0_test E1_main_merged E1_null_merged E2_sens_unmerged E2_null_unmerged; do [ ! -e $R/$d ] || { echo "STOP: $d exists"; exit 1; }; done
source /opt/exp_soft/miniconda3-py3.12/etc/profile.d/conda.sh; set +u; conda activate /home/flljus002/.conda/envs/drimseq_env; set -u
RS=$(which Rscript)
for d in N4_tie_rule_EV N4_tie_rule_cellsRNA N5_binomial_EV N5_binomial_cellsRNA U4_tie_rule_unmerged_EV U4_tie_rule_unmerged_cellsRNA U5_binomial_unmerged_EV U5_binomial_unmerged_cellsRNA; do
  (cd $R/$d && md5sum -c --quiet *_outputs.md5) || { echo "STOP: $d changed"; exit 1; }; done
(cd $R/M1_annotation_merged && md5sum -c --quiet FROZEN_BIN100M_ANNOT.md5); (cd $R/B1_annotation/outputs && md5sum -c --quiet FROZEN_BIN100_ANNOT.md5)
$PY $B/scripts/test_ev_vs_cell100.py $B/scripts $R/E0_test $RS | tee $R/E0_test_result.txt
grep -q ALL_EV_VS_CELL_TESTS_PASS $R/E0_test_result.txt || { echo "STOP: E0 failed"; exit 1; }
run() {  # name tablesEV tablesCell annotation hotEV hotCell
  Rscript $B/scripts/ev_vs_cell100.R $R/$2 $R/$3 $4 $B/manifest_EV.tsv $B/manifest_Cell_sRNA.tsv $R/$5 $R/$6 $R/$1 real 20261009 | tee $R/${1}_log.txt
  $PY $B/scripts/check_ev_vs_cell100.py $R/$2 $R/$3 $4 $B/manifest_EV.tsv $B/manifest_Cell_sRNA.tsv $R/$5 $R/$6 $R/$1 | tee $R/${1}_check_log.txt
  Rscript $B/scripts/ev_vs_cell100.R $R/$2 $R/$3 $4 $B/manifest_EV.tsv $B/manifest_Cell_sRNA.tsv $R/$5 $R/$6 $R/${7} null_sim 20261009 | tee $R/${7}_log.txt
  for d in $1 $7; do (cd $R/$d && md5sum * > ../${d}.md5 && mv ../${d}.md5 ${d}_outputs.md5 && chmod a-w *); done
}
run E1_main_merged N4_tie_rule_EV N4_tie_rule_cellsRNA $R/M1_annotation_merged/bin100m_master_annotation.tsv N5_binomial_EV N5_binomial_cellsRNA E1_null_merged
run E2_sens_unmerged U4_tie_rule_unmerged_EV U4_tie_rule_unmerged_cellsRNA $R/B1_annotation/outputs/bin100_master_annotation.tsv U5_binomial_unmerged_EV U5_binomial_unmerged_cellsRNA E2_null_unmerged
echo E_DONE_OK
