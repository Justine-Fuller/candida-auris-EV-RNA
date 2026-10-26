#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=bin100_C34_cell
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.out
#SBATCH --error=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.err
# C3 + C4 (Cell sRNA): validation gates/tables with ties unresolved (C3), then the same tie rule as EV (C4: neighbouring
# exact ties +1 to both bins; non-neighbouring ties and ambiguous reads unassigned). No statistics. Refuses to overwrite.
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3; ANN=$B/run/B1_annotation/outputs/bin100_master_annotation.tsv
C2=$B/run/C2_counts_cellsRNA; C3=$B/run/C3_tables_cellsRNA; C4=$B/run/C4_tie_rule_cellsRNA
[ ! -e $C3 ] && [ ! -e $C4 ] || { echo "STOP: exists"; exit 1; }
cd $C2 && md5sum *.gz *.json > C2_outputs.md5 && chmod a-w *
$PY $B/scripts/build_tables100.py $ANN $C2 $B/manifest_Cell_sRNA.tsv $C3 sRNA | tee $B/run/C3_log.txt
(cd $C3 && md5sum * > C3_outputs.md5 && chmod a-w *)
$PY $B/scripts/apply_tie_rule100.py $ANN $C2 $B/manifest_Cell_sRNA.tsv $C4 | tee $B/run/C4_log.txt
(cd $C4 && md5sum * > C4_outputs.md5 && chmod a-w *)
echo C34_DONE_OK
