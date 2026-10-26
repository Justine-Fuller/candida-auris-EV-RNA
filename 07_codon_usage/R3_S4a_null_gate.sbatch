#!/bin/bash
#SBATCH -p ada -c 1 --mem=8G -t 01:00:00 -J c6r3_S4a -o /scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6C_codon/logs/%x_%j.out
# gate: the r3 null generator with the r2 real set reproduces the r2 null files byte for byte
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6C_codon; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3; cd $B/run_r3; [ -e S4_gate_r2 ] && exit 1
$PY $B/scripts_r3/make_null_foregrounds.py --fasta $B/refs/Homo_sapiens.GRCh38.cds.all.fa.gz --primary_only   --real_foreground $B/run_r2/S2_real_scoring/foreground_top10pct.txt --real_background $B/run_r2/S2_real_scoring/background_all_scored.txt --outdir S4_gate_r2
{ for f in null_sets.tsv null_foregrounds.tsv; do cmp -s S4_gate_r2/$f $B/run_r2/S4_null_sets/$f && r=IDENTICAL || r=DIFFERENT; echo -e "$f\t$r"; done; } > S4_gate_r2/gate.tsv
cat S4_gate_r2/gate.tsv; grep -q DIFFERENT S4_gate_r2/gate.tsv && { echo GATE FAIL; exit 2; }
rm S4_gate_r2/null_foregrounds.tsv; echo "GATE r3 null code == r2: PASS"
