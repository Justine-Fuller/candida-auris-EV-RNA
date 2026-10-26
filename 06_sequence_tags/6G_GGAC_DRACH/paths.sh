# Step 6G input paths (read-only; md5s in inputs.md5)
P=/scratch/flljus002/Phase3D_forensic_audit
V=$P/v2_downstream_2026-10-12
MAN6E=$V/step6DE/inputs/manifest_6E.tsv
K1=$V/step6DE/run/K1_reads
K2FROZEN=$V/step6DE/run/K2_counts
SIMK1=$V/step6DE/run/G0_tests/sim/K1
SIMMAN=$V/step6DE/run/G0_tests/sim/manifest.tsv
E1=$P/v2_bin100_option1/run/E1_main_merged
ANN=$P/v2_bin100_option1/run/M1_annotation_merged/bin100m_master_annotation.tsv
GENOME=/scratch/flljus002/CandidaEV_Project/_Archive_2026-08-31/01_Differential_expression_analysis/01_Reference_files/C_auris_B8441_version_s02-m01-r01_chromosomes.fasta
F3=$V/step6F_whole_RNA_tag/run/F3_gene_tests
SEQ=$V/step6F_whole_RNA_tag/run/F1_sequences/gene_sequences.tsv.gz
PY=/opt/exp_soft/miniconda3-py3.12/bin/python3
