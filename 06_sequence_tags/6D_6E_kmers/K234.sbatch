#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=s6e_K234
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6DE/logs/%x_%j.out
# Step 6E K2 (GC matching + 6-mer counts) -> K3 (edgeR QL word-share tests, cross-check, label-swap null) -> K4 (Q5 candidates + Q6).
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6DE; P=/scratch/flljus002/Phase3D_forensic_audit; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3
R=$B/run; for d in K2_counts K3_tests K4_candidates_q6; do [ ! -e $R/$d ] || { echo "STOP: $d exists"; exit 1; }; done
(cd $B/scripts && md5sum -c --quiet SCRIPTS_r2.md5)
[ $(ls $R/K1_reads/*.k1.json.gz | wc -l) -eq 24 ] || { echo "STOP: K1 incomplete"; exit 1; }
(cd $P/v2_bin100_option1/run/E1_main_merged && md5sum -c --quiet E1_main_merged_outputs.md5)
(cd $P/v2_bin100_option1/run/M1_annotation_merged && md5sum -c --quiet FROZEN_BIN100M_ANNOT.md5)
$PY $B/scripts/kmer_match6e.py $B/inputs/manifest_6E.tsv $R/K1_reads $R/K2_counts 3 > $R/K2_log.txt
source /opt/exp_soft/miniconda3-py3.12/etc/profile.d/conda.sh; set +u; conda activate /home/flljus002/.conda/envs/drimseq_env; set -u
Rscript $B/scripts/kmer_edger6e.R $B/inputs/manifest_6E.tsv $R/K2_counts $R/K3_tests | tee $R/K3_log.txt
$PY $B/scripts/k4_candidates_q6.py $R/K3_tests $P/v2_bin100_option1/run/E1_main_merged $P/v2_bin100_option1/run/M1_annotation_merged/bin100m_master_annotation.tsv \
   /scratch/flljus002/CandidaEV_Project/_Archive_2026-08-31/01_Differential_expression_analysis/01_Reference_files/C_auris_B8441_version_s02-m01-r01_chromosomes.fasta \
   $R/K4_candidates_q6 | tee $R/K4_log.txt
for d in K2_counts K3_tests K4_candidates_q6; do (cd $R/$d && md5sum * > ../${d}.md5 && mv ../${d}.md5 ${d}_outputs.md5 && chmod a-w *); done
echo K234_DONE
