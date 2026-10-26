#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=s6h_T2
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6H_DRACH_position/logs/%x_%j.out
# Step 6H T2 r2 (Justine option a): synthetic gates incl. the two null checks -> genome check -> real gc / pos / posgc -> table.
set -eo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6H_DRACH_position; cd $B; source inputs/paths.sh
md5sum -c --quiet inputs/inputs.md5; (cd scripts && md5sum -c --quiet SCRIPTS_r2.md5); R=run
for d in H0_sim H0_sim_nojitter H2_sim H3_genome H2_gc H2_pos H2_posgc G2_gc G2_pos G2_posgc H4_table; do [ ! -e $R/$d ] || { echo "STOP: $R/$d exists"; exit 1; }; done
[ $(grep -c $'\tPASS\t' $R/H1_gate.tsv) -eq 24 ] || { echo "STOP: H1 gate not 24 PASS"; exit 1; }
source /opt/exp_soft/miniconda3-py3.12/etc/profile.d/conda.sh; conda activate /home/flljus002/.conda/envs/drimseq_env
$PY scripts/h0_sim6h.py $R/H0_sim; $PY scripts/h0_sim6h.py $R/H0_sim_nojitter nojitter; mkdir $R/H2_sim
$PY scripts/h2_match6h.py $R/H0_sim/manifest.tsv $R/H0_sim/H1 $R/H2_sim/pos 1 pos > /dev/null
$PY scripts/h2_match6h.py $R/H0_sim/manifest.tsv $R/H0_sim/H1 $R/H2_sim/posgc 1 posgc > /dev/null
$PY scripts/h2_match6h.py $R/H0_sim_nojitter/manifest.tsv $R/H0_sim_nojitter/H1 $R/H2_sim/pos_nojitter 1 pos > /dev/null
for m in pos posgc pos_nojitter; do Rscript scripts/g2_edger6g.R $([ $m = pos_nojitter ] && echo $R/H0_sim_nojitter/manifest.tsv || echo $R/H0_sim/manifest.tsv) $R/H2_sim/$m $R/H2_sim/G2_$m PCG_CDS 1007 > $R/H2_sim/G2_$m.log; done
$PY scripts/h4_summary6h.py sim $R/H2_sim/G2_pos $R/H2_sim/G2_posgc $R/H2_sim/G2_pos_nojitter $R/H2_sim/sim_gates.tsv
$PY scripts/h3_genome6h.py $AM $GENOME $R/H3_genome
for m in gc pos posgc; do
  $PY scripts/h2_match6h.py $MAN6E $R/H1_positions $R/H2_$m 3 $m > $R/H2_$m.log
  Rscript scripts/g2_edger6g.R $MAN6E $R/H2_$m $R/G2_$m PCG_CDS 1007 > $R/G2_$m.log; done
$PY scripts/h4_summary6h.py real $G6K5 $R/G2_gc $R/G2_pos $R/G2_posgc $R/H4_table
for d in H0_sim H0_sim_nojitter H2_sim H3_genome H2_gc H2_pos H2_posgc G2_gc G2_pos G2_posgc H4_table; do (cd $R/$d && find . -type f ! -name '*.md5' -exec md5sum {} + > ${d}_outputs.md5); chmod -R a-w $R/$d; done
chmod a-w $R/*.log; echo T2_DONE
