#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=s6e_G0
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6DE/logs/%x_%j.out
# Step 6E gates before the real run: known-answer counter/matching tests + synthetic planted / GC-shift / null through K2 + K3.
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6DE; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3
[ ! -e $B/run/G0_tests ] || { echo "STOP: G0_tests exists"; exit 1; }
mkdir -p $B/run/G0_tests
(cd $B/scripts && md5sum -c --quiet SCRIPTS_r1.md5)
$PY $B/scripts/test_kmer6e.py $B/scripts $B/run/G0_tests/known_answer | tee $B/run/G0_tests/known_answer_result.txt
grep -q ALL_6E_COUNTER_TESTS_PASS $B/run/G0_tests/known_answer_result.txt || { echo "STOP: known-answer failed"; exit 1; }
$PY $B/scripts/sim6e.py $B/run/G0_tests/sim
$PY $B/scripts/kmer_match6e.py $B/run/G0_tests/sim/manifest.tsv $B/run/G0_tests/sim/K1 $B/run/G0_tests/sim/K2 3
source /opt/exp_soft/miniconda3-py3.12/etc/profile.d/conda.sh; set +u; conda activate /home/flljus002/.conda/envs/drimseq_env; set -u
Rscript $B/scripts/kmer_edger6e.R $B/run/G0_tests/sim/manifest.tsv $B/run/G0_tests/sim/K2 $B/run/G0_tests/sim/K3 PCG_CDS
$PY $B/scripts/gate6e_sim.py $B/run/G0_tests/sim/K3 | tee $B/run/G0_tests/sim_gates_result.txt
echo G0_DONE
