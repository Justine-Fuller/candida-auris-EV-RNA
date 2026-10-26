#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=s6g_S2
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6G_GGAC_DRACH/logs/%x_%j.out
# Step 6G stretch level: synthetic gates first (stop on failure), then real + random-label null.
set -eo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6G_GGAC_DRACH; cd $B; source inputs/paths.sh
md5sum -c --quiet inputs/inputs.md5; (cd scripts && md5sum -c --quiet SCRIPTS_r2.md5)
for d in G4_sim G4_stretch; do [ ! -e run/$d ] || { echo "STOP: run/$d exists"; exit 1; }; done
$PY scripts/g4_stretch6g.py sim $E1 $ANN $GENOME run/G4_sim 20261014 10 100000 0
$PY scripts/g4_stretch6g.py real $E1 $ANN $GENOME run/G4_stretch 20261014 10 100000 100
for d in G4_sim G4_stretch; do (cd run/$d && md5sum * > ${d}_outputs.md5); chmod -R a-w run/$d; done; echo S2_DONE
