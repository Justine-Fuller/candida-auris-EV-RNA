#!/bin/bash
#SBATCH -p ada -c 1 --mem=8G -t 01:00:00 -J c6r3_S1 -o /scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6C_codon/logs/%x_%j.out
# r3 (Justine 2026-10-13, Q2 a): 9 EV tRNA families after the CCA recount. Gates: known-answer; r3 code with the r2 set reproduces r2.
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6C_codon; S=$B/scripts_r3; R2=$B/run_r2; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3
V1=/scratch/flljus002/CandidaEV_Project/05_EV_coverage_analysis/11_Full_length_biological_annotation/output/codon_usage_analysis
export R_LIBS_USER=/home/flljus002/Rlibs
FA=$B/refs/Homo_sapiens.GRCh38.cds.all.fa.gz
REAL9="Tyr:UAC;Arg:CGG,CGU;Ala:GCU,GCG;His:CAC;Phe:UUC;Leu:UUA;Thr:ACG"
(cd $S && md5sum -c --quiet SCRIPTS_r3.md5)
mkdir -p $B/run_r3; cd $B/run_r3
for d in S1_known_answer S1_r2_reproduction S1_codon_preference S2_real_scoring; do [ -e $d ] && { echo "$d exists"; exit 1; }; done
mkdir S1_known_answer S1_r2_reproduction S1_codon_preference
$PY $S/test_codon_scorer.py > S1_known_answer/known_answer_test.txt 2>&1 || { cat S1_known_answer/known_answer_test.txt; exit 2; }
$PY $S/codon_scorer.py --fasta $FA --codons "Tyr:UAC;Arg:CGG;Ala:GCU;His:CAC;Phe:UUC;Leu:UUA" --outdir S1_r2_reproduction/real_set --primary_only
{ for f in human_gene_codon_scores.tsv foreground_top10pct.txt background_all_scored.txt; do
    cmp -s S1_r2_reproduction/real_set/$f $R2/S2_real_scoring/$f && r=IDENTICAL || r=DIFFERENT; echo -e "$f\t$r"; done; } > S1_r2_reproduction/gate.tsv
cat S1_r2_reproduction/gate.tsv; grep -q DIFFERENT S1_r2_reproduction/gate.tsv && { echo "GATE r2 reproduction FAIL"; exit 2; }
$PY $S/codon_preference.py --catalogue $B/inputs/family_catalogue_43.tsv --cauris $V1/c_auris_codon_usage.tsv --outdir S1_codon_preference --chosen r3
Rscript $S/fisher_preference.R S1_codon_preference/fisher_counts.tsv S1_codon_preference/fisher_results.tsv
$PY $S/codon_scorer.py --fasta $FA --codons "$REAL9" --outdir S2_real_scoring --primary_only
echo S1_R3_DONE
