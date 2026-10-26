#!/bin/bash
# independent samtools faidx check of 3 region sequences; usage: handcheck.sh <run_dir>
set -euo pipefail
module load software/samtools-1.20
GENOME=/scratch/flljus002/CandidaEV_Project/_Archive_2026-08-31/01_Differential_expression_analysis/01_Reference_files/C_auris_B8441_version_s02-m01-r01_chromosomes.fasta
R=$1; O=$2
# independent check: first + strand region, first - strand region, first region with >1 genomic block (spliced), via samtools faidx
T=$R/prioritised_regions_q0.05.tsv
H=$(head -1 $T | tr '\t' '\n' | grep -n -x -E 'region_id|chrom|strand|genomic_blocks_5to3|sequence' | cut -d: -f1 | paste -sd,)
{ awk -F'\t' 'NR>1 && $0 ~ /\t\+\t/' <(cut -f$H $T) | head -1 || true
  awk -F'\t' 'NR>1 && $0 ~ /\t-\t/' <(cut -f$H $T) | head -1 || true
  awk -F'\t' 'NR>1 && $4 ~ /;/' <(cut -f$H $T) | head -1; } > $O/handcheck_input.tsv
echo -e "region_id\tstrand\tblocks\tsamtools_seq_md5\ttable_seq_md5\tmatch\tfirst30\tlast30" > $O/handcheck_sequences.tsv
while IFS=$'\t' read rid chrom strand blocks seq; do
  s=""
  for blk in ${blocks//;/ }; do a=${blk%-*}; b=${blk#*-}
    if [ "$strand" = "-" ]; then x=$(samtools faidx -i $GENOME "$chrom:$((a+1))-$b" | grep -v '>' | tr -d '\n')
    else x=$(samtools faidx $GENOME "$chrom:$((a+1))-$b" | grep -v '>' | tr -d '\n'); fi
    s="$s$x"; done
  s=$(echo "$s" | tr a-z A-Z)
  m1=$(echo -n "$s" | md5sum | cut -c1-12); m2=$(echo -n "$seq" | md5sum | cut -c1-12)
  echo -e "$rid\t$strand\t$blocks\t$m1\t$m2\t$([ "$m1" = "$m2" ] && echo PASS || echo FAIL)\t${s:0:30}\t${s: -30}" >> $O/handcheck_sequences.tsv
done < $O/handcheck_input.tsv
cat $O/handcheck_sequences.tsv
[ $(grep -c PASS$'\t' $O/handcheck_sequences.tsv) -eq $(wc -l < $O/handcheck_input.tsv) ] && echo HANDCHECK_PASS || echo HANDCHECK_FAIL
