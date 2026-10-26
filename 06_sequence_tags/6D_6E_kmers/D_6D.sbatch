#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=s6d
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6DE/logs/%x_%j.out
# Step 6D: prepare (universe = PCGs detected in B8441 caspofungin EVs, >= 2/3 replicates) -> annotation table (goslim_yeast words)
# -> s5b_enrich.R r5 (unchanged, md5 176c2953...) standard + read-depth-corrected, covariate mean EV TPM and mean EV reads.
set -euo pipefail
export R_LIBS_USER=/home/flljus002/Rlibs
D=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12; B=$D/step6DE; R=$B/run; P=/scratch/flljus002/Phase3D_forensic_audit; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3
for d in D1_prepare D2_table D3_enrich_TPM D3_enrich_reads; do [ ! -e $R/$d ] || { echo "STOP: $d exists"; exit 1; }; done
echo "176c295374884e1923e814d065d1d03a  $D/scripts/s5b_enrich.R" | md5sum -c -
(cd $B/scripts && md5sum -c --quiet SCRIPTS_r3.md5)
INT=/scratch/flljus002/CandidaEV_Project/06_PCG_functional_annotation/08_Integrated_annotation/PCG_functional_annotation/PCG_functional_annotation.tsv
V3=$P/v2_rRNA_2026-10-09/refs/raw/genome_annotations/GCA_002759435.3_Cand_auris_B8441_V3_genomic.gff.gz
md5sum $P/v2_rRNA_2026-10-09/run/FL_r3/unit_table.tsv.gz $P/v2_rRNA_2026-10-09/run/TPM_abundance/tpm_with_rRNA_per_library.tsv.gz $INT $V3 \
  $D/refs/kegg_pathway_ko_list.tsv $B/inputs/6D/* $P/v2_bin100_option1/run/B4_tie_rule/B8441_EV_CASP_*.gz $P/v2_rRNA_2026-10-09/run/FL3_counts/B8441_EV_CASP_*.json > $B/provenance/s6d_inputs_md5.txt
$PY $B/scripts/s6d_prepare.py $P/v2_rRNA_2026-10-09/run/FL_r3/unit_table.tsv.gz $P/v2_rRNA_2026-10-09/run/TPM_abundance/tpm_with_rRNA_per_library.tsv.gz \
    $P/v2_bin100_option1/run/B4_tie_rule $P/v2_rRNA_2026-10-09/run/FL3_counts $INT $V3 $R/D1_prepare | tee $R/D1_log.txt
Rscript $B/scripts/s6d_table.R $R/D1_prepare/fl_genes.tsv $B/inputs/6D/goslim_yeast.obo $B/inputs/6D/v1_go_to_slim_mapping.tsv \
    $B/inputs/6D/full_length_mRNA_go_annotation_table.tsv $B/inputs/6D/eggnog_15_genes.tsv $R/D2_table | tee $R/D2_log.txt
grep -q S6D_TABLE_GATES_PASS $R/D2_log.txt || { echo "STOP: 6D table gate failed"; exit 1; }
Rscript $D/scripts/s5b_enrich.R --selftest | tee $R/D3_selftest_log.txt
grep -q S5B_SELFTEST_PASS $R/D3_selftest_log.txt || { echo "STOP: selftest failed"; exit 1; }
for c in TPM reads; do Rscript $D/scripts/s5b_enrich.R $R/D1_prepare/prepare_$c $D/refs/kegg_pathway_ko_list.tsv $R/D3_enrich_$c 10000 | tee $R/D3_${c}_log.txt; done
for d in D1_prepare D2_table D3_enrich_TPM D3_enrich_reads; do (cd $R/$d && find . -type f ! -name "*_outputs.md5" -exec md5sum {} + > ../${d}.md5 && mv ../${d}.md5 ${d}_outputs.md5 && find . -type f -exec chmod a-w {} +); done
echo S6D_DONE
