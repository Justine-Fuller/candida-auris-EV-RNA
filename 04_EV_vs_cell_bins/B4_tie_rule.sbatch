#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=bin100_B4_tierule
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.out
#SBATCH --error=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.err
# B4: apply Caro's confirmed tie rule (+1 to both bins of an exact 2-bin greatest-overlap tie; ambiguous_model reads stay
# unassigned) to the validated B2 counts. No BAM recount, no statistics. Freezes B2/B3 first. Refuses to overwrite.
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1; O=$B/run/B4_tie_rule; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3
[ ! -e $O ] && [ ! -e $B/run/B4_tie_rule_test ] || { echo "STOP: exists"; exit 1; }
cd $B/run/B2_counts; [ -e B2_outputs.md5 ] || md5sum *.gz *.json > B2_outputs.md5; md5sum -c B2_outputs.md5 --quiet; chmod a-w *
cd $B/run/B3_tables; md5sum -c B3_outputs.md5 --quiet; chmod a-w *
$PY $B/scripts/test_tie_rule100.py $B/scripts $B/run/B4_tie_rule_test | tee $B/run/B4_tie_rule_test_result.txt
$PY $B/scripts/apply_tie_rule100.py $B/run/B1_annotation/outputs/bin100_master_annotation.tsv $B/run/B2_counts $B/manifest_EV.tsv $O | tee $B/run/B4_log.txt
(cd $O && md5sum * > B4_outputs.md5 && chmod a-w *)
echo B4_DONE_OK
