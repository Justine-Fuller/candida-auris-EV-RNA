#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=s6e_K4b
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6DE/logs/%x_%j.out
# Step 6E K4b (2026-10-14, Justine): Q6 stretch test on the no-GC-control (naive) tag candidates. Gate first: the r2 script with
# "primary" must reproduce the frozen K4 outputs byte for byte.
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6DE; P=/scratch/flljus002/Phase3D_forensic_audit; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3; R=$B/run
(cd $B/scripts && md5sum -c --quiet SCRIPTS_r4.md5)
(cd $R/K4_candidates_q6 && md5sum -c --quiet K4_candidates_q6_outputs.md5); (cd $R/K3_tests && md5sum -c --quiet K3_tests_outputs.md5)
D=$R/K4b_q6_naive; [ ! -e $D ] || { echo "STOP: exists"; exit 1; }; mkdir $D
ARGS="$R/K3_tests $P/v2_bin100_option1/run/E1_main_merged $P/v2_bin100_option1/run/M1_annotation_merged/bin100m_master_annotation.tsv /scratch/flljus002/CandidaEV_Project/_Archive_2026-08-31/01_Differential_expression_analysis/01_Reference_files/C_auris_B8441_version_s02-m01-r01_chromosomes.fasta"
$PY $B/scripts/k4_candidates_q6_r2_select.py $ARGS $D/repro_primary primary | tee $D/repro_primary.log
for f in tag_candidates.tsv q6_fisher_ev_selective_bins.tsv significant_words_per_group.tsv; do
  cmp $D/repro_primary/$f $R/K4_candidates_q6/$f && echo "GATE_REPRO_PASS $f" || { echo "GATE_REPRO_FAIL $f"; exit 1; }; done
$PY $B/scripts/k4_candidates_q6_r2_select.py $ARGS $D/naive naive | tee $D/naive.log
(cd $D && find . -type f -exec md5sum {} + > ../K4b_outputs.md5 && mv ../K4b_outputs.md5 .); chmod -R a-w $D; echo K4b_DONE
