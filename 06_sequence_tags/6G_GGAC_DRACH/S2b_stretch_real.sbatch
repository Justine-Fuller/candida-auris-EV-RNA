#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=s6g_S2b
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6G_GGAC_DRACH/logs/%x_%j.out
# Step 6G stretch level, REAL data (Justine 2026-10-15 07:04 UTC, option a): run with the r2 GC correction (5 GC classes) although
# the synthetic GC-shift gate let 1/257 units through (CTCA; GGAC/DRACH unaffected) -- reported as a caveat. Random-label null gate
# per set still applies.
set -eo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6G_GGAC_DRACH; cd $B; source inputs/paths.sh
md5sum -c --quiet inputs/inputs.md5; (cd scripts && md5sum -c --quiet SCRIPTS_r3.md5)
[ ! -e run/G4_stretch ] || { echo "STOP: run/G4_stretch exists"; exit 1; }
set +e; $PY scripts/g4_stretch6g.py real $E1 $ANN $GENOME run/G4_stretch 20261014 10 100000 100; rc=$?; set -e
(cd run/G4_stretch && md5sum * > G4_stretch_outputs.md5); chmod -R a-w run/G4_stretch; echo "G4 real exit $rc"; exit $rc
