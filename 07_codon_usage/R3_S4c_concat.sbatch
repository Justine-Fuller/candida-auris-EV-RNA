#!/bin/bash
#SBATCH -p ada -c 1 --mem=4G -t 00:30:00 -J c6r3_S4c -o /scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6C_codon/logs/%x_%j.out
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6C_codon; cd $B/run_r3/S4_null_sets
[ $(ls null_generation_QC_part*.txt | wc -l) -eq 24 ] && ! grep -L "NULL GENERATION: PASS" null_generation_QC_part*.txt | grep -q . || { echo "STOP part QC"; exit 1; }
{ head -1 null_sets_part001.tsv; cat null_sets_part*.tsv | grep -v "^set_id" | sort; } > null_sets.tsv
cat null_foregrounds_part*.tsv | sort > null_foregrounds.tsv
echo "sets: $(($(wc -l < null_sets.tsv)-1)); foregrounds: $(wc -l < null_foregrounds.tsv)"
rm null_foregrounds_part*.tsv null_sets_part*.tsv
