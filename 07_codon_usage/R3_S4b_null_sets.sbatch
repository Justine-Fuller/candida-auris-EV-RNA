#!/bin/bash
#SBATCH -p ada -c 1 --mem=8G -t 03:00:00 -J c6r3_S4b --array=1-24 -o /scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6C_codon/logs/%x_%A_%a.out
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6C_codon; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3; cd $B/run_r3; mkdir -p S4_null_sets
REAL9="Tyr:UAC;Arg:CGG,CGU;Ala:GCU,GCG;His:CAC;Phe:UUC;Leu:UUA;Thr:ACG"
$PY $B/scripts_r3/make_null_foregrounds.py --fasta $B/refs/Homo_sapiens.GRCh38.cds.all.fa.gz --primary_only --real "$REAL9"   --real_foreground S2_real_scoring/foreground_top10pct.txt --real_background S2_real_scoring/background_all_scored.txt --outdir S4_null_sets   --task $SLURM_ARRAY_TASK_ID --ntasks 24
