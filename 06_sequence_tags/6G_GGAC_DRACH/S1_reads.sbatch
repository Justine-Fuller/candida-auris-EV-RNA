#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=s6g_S1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6G_GGAC_DRACH/logs/%x_%j.out
# Step 6G read level: known-answer tests -> 6E K2 reproduction gate (k6) -> synthetic gates (6E SIM libraries) -> real k4 + DRACH.
set -eo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6G_GGAC_DRACH; cd $B; source inputs/paths.sh
md5sum -c --quiet inputs/inputs.md5; (cd scripts && md5sum -c --quiet SCRIPTS_r1.md5)
R=run; for d in G0_known_answer G1_k6_repro G1_sim G2_sim G3_sim G1_k4 G1_k5drach G2_k4 G2_k5drach G3_reads; do [ ! -e $R/$d ] || { echo "STOP: $R/$d exists"; exit 1; }; done
$PY scripts/test_g1_6g.py scripts $R/G0_known_answer | tee $R/G0_known_answer.txt; grep -q ALL_G1_COUNTER_TESTS_PASS $R/G0_known_answer.txt
$PY scripts/g5_sets6g.py selftest | tee -a $R/G0_known_answer.txt; grep -q F4_SELFTEST_PASS $R/G0_known_answer.txt
# --- gate: unit k6 must reproduce the frozen 6E K2 counts exactly (same GC-matched draws)
$PY scripts/g1_count6g.py $MAN6E $K1 $R/G1_k6_repro 3 k6 > $R/G1_k6_repro.log
n=0; for f in $K2FROZEN/*.counts.tsv.gz; do b=$(basename $f); cmp <(zcat $f) <(zcat $R/G1_k6_repro/$b) || { echo "REPRO_FAIL $b"; exit 1; }; n=$((n+1)); done
for f in $K2FROZEN/*.gc.tsv; do cmp $f $R/G1_k6_repro/$(basename $f) || { echo "REPRO_FAIL $(basename $f)"; exit 1; }; done
echo "GATE_K6_REPRO_PASS $n count files + gc tables identical to 6E K2" | tee $R/G1_k6_repro_gate.txt
source /opt/exp_soft/miniconda3-py3.12/etc/profile.d/conda.sh; conda activate /home/flljus002/.conda/envs/drimseq_env
# --- synthetic gates on the frozen 6E SIM libraries (planted GGACTA, GC shift, null)
mkdir $R/G1_sim $R/G2_sim
$PY scripts/g1_count6g.py $SIMMAN $SIMK1 $R/G1_sim/k4 3 k4 > /dev/null; $PY scripts/g1_count6g.py $SIMMAN $SIMK1 $R/G1_sim/k5drach 3 k5drach > /dev/null
Rscript scripts/g2_edger6g.R $SIMMAN $R/G1_sim/k4 $R/G2_sim/k4 PCG_CDS 256 > $R/G2_sim_k4.log
Rscript scripts/g2_edger6g.R $SIMMAN $R/G1_sim/k5drach $R/G2_sim/k5drach PCG_CDS 1007 > $R/G2_sim_k5drach.log
$PY scripts/g3_reads_summary6g.py sim $R/G2_sim/k4 $R/G2_sim/k5drach $R/G3_sim
# --- real
$PY scripts/g1_count6g.py $MAN6E $K1 $R/G1_k4 3 k4 > $R/G1_k4.log; $PY scripts/g1_count6g.py $MAN6E $K1 $R/G1_k5drach 3 k5drach > $R/G1_k5drach.log
Rscript scripts/g2_edger6g.R $MAN6E $R/G1_k4 $R/G2_k4 PCG_CDS,lncRNA,tRNA,5S_rRNA 256 > $R/G2_k4.log
Rscript scripts/g2_edger6g.R $MAN6E $R/G1_k5drach $R/G2_k5drach PCG_CDS,lncRNA,tRNA,5S_rRNA 1007 > $R/G2_k5drach.log
$PY scripts/g3_reads_summary6g.py real $R/G2_k4 $R/G2_k5drach $R/G3_reads
for d in G0_known_answer G1_k6_repro G1_sim G2_sim G3_sim G1_k4 G1_k5drach G2_k4 G2_k5drach G3_reads; do
  (cd $R/$d && find . -type f ! -name '*.md5' -exec md5sum {} + > ${d}_outputs.md5); chmod -R a-w $R/$d; done
chmod a-w $R/*.log $R/*.txt; echo S1_DONE
