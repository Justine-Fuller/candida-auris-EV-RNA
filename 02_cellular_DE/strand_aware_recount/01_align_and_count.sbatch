#!/bin/bash
#SBATCH --job-name=cellmRNA_stranded
#SBATCH --partition=ada
#SBATCH --array=1-12%3
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/flljus002/CandidaEV_Project/08_Cell_mRNA_corrected_2026-09-25/logs/%x_%A_%a.out
# One cellular mRNA library per array task:
#   1. STAR alignment with ONE set of settings for all 12 libraries (04_STAR_and_counting_PROPOSAL.md),
#      replacing the old treatment-dependent settings.
#   2. bedtools coverage against the same final_with_introns.sorted.bed as before, three ways:
#        -S  = reads on the OPPOSITE strand to the feature  -> the real counts (libraries are reverse-stranded,
#              confirmed for all 12 in POLYA_VALIDATION_REPORT.md)
#        -s  = reads on the SAME strand (should be small; sanity check)
#        none = strand ignored (old convention, so the effect of strand alone can be separated from realignment)
# Input: the trimmed FASTQs from 2026-09-25 (correct TruSeq adapter). BAMs are KEPT this time.
set -euo pipefail
B=/scratch/flljus002/CandidaEV_Project/08_Cell_mRNA_corrected_2026-09-25
IDX=$B/02_STAR/index
BED=/scratch/flljus002/CandidaEV_Project/03_BEDTools_coverage/01_Input/final_with_introns.sorted.bed
OUT=$B/03_stranded_2026-10-14
SAMPLES=(B8441_mRNA_1 B8441_mRNA_2 B8441_mRNA_3 B8441_mRNA_CASP_1 B8441_mRNA_CASP_2 B8441_mRNA_CASP_3
         MMC1_mRNA_1 MMC1_mRNA_2 MMC1_mRNA_3 MMC1_mRNA_CASP_1 MMC1_mRNA_CASP_2 MMC1_mRNA_CASP_3)
S=${SAMPLES[$((SLURM_ARRAY_TASK_ID-1))]}
FQ=$B/01_cutadapt/${S}_trimmed.fastq.gz

[ -s $IDX/SA ] || { echo "STOP: STAR index missing, run 00_build_star_index.sbatch first"; exit 1; }
[ -s $FQ ] || { echo "STOP: $FQ missing"; exit 1; }
mkdir -p $OUT/bam $OUT/counts $OUT/logs
[ ! -e $OUT/counts/${S}.stranded_S.txt ] || { echo "STOP: $S already counted"; exit 1; }
free=$(quota -w | awk '/hpc-scratch/{printf "%d", ($3-$2)/1048576}'); echo "free_GB=$free"
[ "$free" -ge 5 ] || { echo "STOP: scratch free < 5 GB"; exit 1; }
echo "host $(hostname)  sample $S"

module load software/STAR-2.7.11b
STAR --runThreadN 6 --genomeDir $IDX \
  --readFilesIn $FQ --readFilesCommand zcat \
  --outFileNamePrefix $OUT/bam/${S}_ \
  --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 8000000000 \
  --outSAMattributes NH HI AS nM NM MD \
  --outSAMunmapped Within \
  --outSAMprimaryFlag OneBestScore \
  --outFilterType BySJout \
  --outFilterMultimapNmax 10 \
  --outFilterMismatchNmax 3 \
  --outFilterScoreMinOverLread 0.66 --outFilterMatchNminOverLread 0.66 \
  --outFilterMatchNmin 25 \
  --alignIntronMin 20 --alignIntronMax 2000 \
  --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 \
  --alignEndsType Local \
  --twopassMode None
rm -rf $OUT/bam/${S}__STARtmp
cp $OUT/bam/${S}_Log.final.out $OUT/logs/

module load software/samtools-1.20
BAM=$OUT/bam/${S}_Aligned.sortedByCoord.out.bam
samtools index $BAM
# Count mapped reads only (unmapped reads are kept in the BAM but must not be counted)
samtools view -b -F 4 $BAM > $OUT/bam/${S}_mapped.bam && samtools index $OUT/bam/${S}_mapped.bam

module load software/bedtools-2.31.0
bedtools coverage -split -sorted -S -a $BED -b $OUT/bam/${S}_mapped.bam -counts > $OUT/counts/${S}.stranded_S.txt
bedtools coverage -split -sorted -s -a $BED -b $OUT/bam/${S}_mapped.bam -counts > $OUT/counts/${S}.samestrand_s.txt
bedtools coverage -split -sorted    -a $BED -b $OUT/bam/${S}_mapped.bam -counts > $OUT/counts/${S}.unstranded.txt
rm -f $OUT/bam/${S}_mapped.bam $OUT/bam/${S}_mapped.bam.bai

# Sanity check: on protein-coding _transcript rows, opposite-strand (-S) reads should far outnumber same-strand (-s)
paste <(awk -F'\t' '$4 ~ /^B9J08_.*_transcript$/ {s+=$7} END{print s}' $OUT/counts/${S}.stranded_S.txt) \
      <(awk -F'\t' '$4 ~ /^B9J08_.*_transcript$/ {s+=$7} END{print s}' $OUT/counts/${S}.samestrand_s.txt) |
  awk -v s=$S '{printf "%s\tPCG reads -S=%d  -s=%d  fraction_correct_strand=%.3f\n", s, $1, $2, $1/($1+$2)}' | tee $OUT/logs/${S}_strand_check.txt
md5sum $OUT/counts/${S}.* > $OUT/logs/${S}.counts.md5
echo DONE_OK
