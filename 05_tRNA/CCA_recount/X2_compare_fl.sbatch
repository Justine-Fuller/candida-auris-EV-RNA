#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=cca2_X2
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_tRNA_2026-10-09/cca_recount_2026-10-13/logs/%x_%j.out
# Gates + old/new CCA table + FL r3 rebuilt with recounted CCA (incl / excl templated-ambiguous) + comparison. Refuses to overwrite.
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_tRNA_2026-10-09; C=$B/cca_recount_2026-10-13; R=$C/run; P=/scratch/flljus002/Phase3D_forensic_audit; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3
for d in R2_compare FL_r3_cca2_incl FL_r3_cca2_excl; do [ ! -e $R/$d ] || { echo "STOP $d exists"; exit 1; }; done
[ $(ls $R/R1_cca2/*.cca2.json | wc -l) -eq 36 ] || { echo STOP incomplete; exit 1; }
(cd $C/scripts && md5sum -c --quiet SCRIPTS_r1.md5)
$PY $C/scripts/cca2_compare.py $R/R1_cca2 $B/run/T1_all $B/manifest_all36.tsv $R/R2_compare | tee $R/R2_log.txt
grep -q CCA2_COMPARE_DONE $R/R2_log.txt || { echo "STOP gate"; exit 1; }
for v in incl excl; do
  $PY $P/v2_rRNA_2026-10-09/scripts/fl_classify_r3.py $R/R2_compare/T1_all_cca2_$v $P/v2_rRNA_2026-10-09/manifest_all36.tsv $B/inputs/bin100_master_annotation.tsv \
    $B/inputs/candidate_trna_reference_168.bed $P/v2_read_census_2026-10-09/run_r2/counts $P/v2_bin100_option1/run/B4_tie_rule $P/v2_bin100_option1/run/C4_tie_rule_cellsRNA \
    $R/FL_r3_cca2_$v $P/v2_bin150/run/R2_counts $P/v2_rRNA_2026-10-09/run/FL3_counts $P/v2_rRNA_2026-10-09/FL3_annotation/new_bins100.tsv > $R/FL_${v}_log.txt
done
$PY $C/scripts/fl_compare_cca2.py $P/v2_rRNA_2026-10-09/run/FL_r3/unit_table.tsv.gz $R/FL_r3_cca2_incl $R/FL_r3_cca2_excl $B/run/A_real $R/R2_compare/fl_and_codon_changes.tsv
for d in R2_compare FL_r3_cca2_incl FL_r3_cca2_excl; do (cd $R/$d && find . -type f ! -name '*_outputs.md5' -exec md5sum {} + > ../$d.md5 && mv ../$d.md5 ${d}_outputs.md5 && find . -type f -exec chmod a-w {} +); done
echo X2_DONE
