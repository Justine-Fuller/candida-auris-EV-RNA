#!/bin/bash
#SBATCH -p ada
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 01:00:00
#SBATCH -J s8_prio
#SBATCH -o /scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/logs/s8_%j.out
# Step 8 (v2) r2 (A3 + B2) prioritisation of EV hotspot regions: self-test, run, independent sequence check (samtools faidx), freeze.
set -euo pipefail
P=/scratch/flljus002/Phase3D_forensic_audit
D=$P/v2_downstream_2026-10-12/step8_prioritise
B=$P/v2_bin100_option1/run
GENOME=/scratch/flljus002/CandidaEV_Project/_Archive_2026-08-31/01_Differential_expression_analysis/01_Reference_files/C_auris_B8441_version_s02-m01-r01_chromosomes.fasta
PY=/opt/exp_soft/miniconda3-py3.12/bin/python3
module load software/samtools-1.20
for d in $D/selftest_r2 $D/run; do [ ! -e $d ] || { echo "STOP: $d exists"; exit 1; }; done
cd $D && md5sum -c scripts/SCRIPTS_r2.md5
md5sum $B/N5_binomial_EV/*.group_calls.tsv.gz $B/E1_main_merged/*.ev_vs_cell_bins.tsv.gz $B/N4_tie_rule_EV/*.bin100_final_table.tsv.gz \
  $B/N4_tie_rule_cellsRNA/*.bin100_final_table.tsv.gz $B/M1_annotation_merged/bin100m_master_annotation.tsv $GENOME \
  $P/v2_downstream_2026-10-12/step5A/run/hotspot_regions_q0.05.tsv.gz $P/v2_downstream_2026-10-12/step5B/prepare/hotspot_genes_annotated.tsv \
  > $D/provenance_s8_inputs_md5.txt
$PY $D/scripts/s8_prioritise.py selftest $D/selftest_r2
$PY $D/scripts/s8_prioritise.py run $B/N5_binomial_EV $B/E1_main_merged $B/N4_tie_rule_EV $B/N4_tie_rule_cellsRNA \
  $B/M1_annotation_merged/bin100m_master_annotation.tsv $GENOME \
  $P/v2_downstream_2026-10-12/step5A/run/hotspot_regions_q0.05.tsv.gz $P/v2_downstream_2026-10-12/step5B/prepare/hotspot_genes_annotated.tsv $D/run
# independent check: first + strand region, first - strand region, first region with >1 genomic block (spliced), via samtools faidx
T=$D/run/prioritised_regions_q0.05.tsv
H=$(head -1 $T | tr '\t' '\n' | grep -n -x -E 'region_id|chrom|strand|genomic_blocks_5to3|sequence' | cut -d: -f1 | paste -sd,)
{ awk -F'\t' 'NR>1 && $0 ~ /\t\+\t/' <(cut -f$H $T) | head -1 || true
  awk -F'\t' 'NR>1 && $0 ~ /\t-\t/' <(cut -f$H $T) | head -1 || true
  awk -F'\t' 'NR>1 && $4 ~ /;/' <(cut -f$H $T) | head -1; } > $D/run/handcheck_input.tsv
echo -e "region_id\tstrand\tblocks\tsamtools_seq_md5\ttable_seq_md5\tmatch\tfirst30\tlast30" > $D/run/handcheck_sequences.tsv
while IFS=$'\t' read rid chrom strand blocks seq; do
  s=""
  for blk in ${blocks//;/ }; do a=${blk%-*}; b=${blk#*-}
    if [ "$strand" = "-" ]; then x=$(samtools faidx -i $GENOME "$chrom:$((a+1))-$b" | grep -v '>' | tr -d '\n')
    else x=$(samtools faidx $GENOME "$chrom:$((a+1))-$b" | grep -v '>' | tr -d '\n'); fi
    s="$s$x"; done
  s=$(echo "$s" | tr a-z A-Z)
  m1=$(echo -n "$s" | md5sum | cut -c1-12); m2=$(echo -n "$seq" | md5sum | cut -c1-12)
  echo -e "$rid\t$strand\t$blocks\t$m1\t$m2\t$([ "$m1" = "$m2" ] && echo PASS || echo FAIL)\t${s:0:30}\t${s: -30}" >> $D/run/handcheck_sequences.tsv
done < $D/run/handcheck_input.tsv
cat $D/run/handcheck_sequences.tsv
[ $(grep -c PASS$'\t' $D/run/handcheck_sequences.tsv) -eq $(wc -l < $D/run/handcheck_input.tsv) ] && echo HANDCHECK_PASS || echo HANDCHECK_FAIL
cd $D/run && md5sum * > ../provenance_s8_outputs_md5.txt
chmod -R a-w $D/run $D/selftest_r2
echo S8_DONE
