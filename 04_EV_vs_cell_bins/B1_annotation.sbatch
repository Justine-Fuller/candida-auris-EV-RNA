#!/bin/bash
#SBATCH --partition=ada
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.out
#SBATCH --error=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1/logs/%x_%j.err
#SBATCH --job-name=bin100_B1_annot
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
# B1: 100-nt bin annotation from the md5-verified frozen B1 inputs (copied read-only from v2_bin150 R1); checks vs frozen
# 150-nt B1 outputs; then freeze read-only. Refuses to overwrite.
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_bin100_option1; F=/scratch/flljus002/Phase3D_forensic_audit/v2_bin150/run/R1_annotation; R=$B/run/B1_annotation
[ ! -e $R ] || { echo "STOP: exists"; exit 1; }
(cd $F/outputs && md5sum -c FROZEN_B1.md5 --quiet) || { echo "STOP: frozen B1 changed"; exit 1; }
mkdir -p $R/inputs $R/scripts
cp $F/inputs/* $R/inputs/
grep inputs/ $F/inputs.expected.md5 > $R/inputs.expected.md5
(cd $R && md5sum -c inputs.expected.md5)
cp $B/scripts/build_bin100.py $R/scripts/
PY=/opt/exp_soft/miniconda3-py3.12/bin/python3; $PY --version > $R/python_version.txt 2>&1
$PY $R/scripts/build_bin100.py | tee $R/build_log.txt
$PY $B/scripts/check_annotation100.py $R/outputs $F/outputs | tee $R/check_log.txt
(cd $R/outputs && md5sum * > FROZEN_BIN100_ANNOT.md5)
chmod a-w $R/outputs/* $R/outputs $R/inputs/* $R/scripts/* $R/annotation_checks.tsv
echo B1_DONE_OK_FROZEN
