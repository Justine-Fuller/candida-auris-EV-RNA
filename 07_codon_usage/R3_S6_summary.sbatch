#!/bin/bash
#SBATCH -p ada -c 1 --mem=32G -t 02:00:00 -J c6r3_S6 -o /scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6C_codon/logs/%x_%j.out
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6C_codon; export R_LIBS_USER=/home/flljus002/Rlibs
cd $B/run_r3; [ -e S6_summary ] && { echo exists; exit 1; }
N=$(($(wc -l < S4_null_sets/null_sets.tsv)-1))
Rscript $B/scripts_r3/summarize_null.R S3_enrichment_real S5_null S4_null_sets/null_sets.tsv S6_summary $N
/opt/exp_soft/miniconda3-py3.12/bin/python3 $B/scripts_r3/compare_runs.py $B/run_r2 $B/run_r3 S7_compare_r2_r3
(cd $B/run_r3 && find . -type f ! -path "./S5_null/*" ! -name "null_foregrounds.tsv" -exec md5sum {} + > $B/provenance/r3_outputs_md5.txt; find . -type f -exec chmod a-w {} +)
echo S6_R3_DONE
