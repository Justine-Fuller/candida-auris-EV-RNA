#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=cca2_X0
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_tRNA_2026-10-09/cca_recount_2026-10-13/logs/%x_%j.out
set -euo pipefail
C=/scratch/flljus002/Phase3D_forensic_audit/v2_tRNA_2026-10-09/cca_recount_2026-10-13; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3
(cd $C/scripts && md5sum -c --quiet SCRIPTS_r1.md5)
$PY $C/scripts/test_cca_recount.py $C/scripts /scratch/flljus002/Phase3D_forensic_audit/v2_tRNA_2026-10-09/scripts/trna_count.py $C/run/X0_test | tee $C/run/X0_test_result.txt
grep -q ALL_CCA_RECOUNT_TESTS_PASS $C/run/X0_test_result.txt
