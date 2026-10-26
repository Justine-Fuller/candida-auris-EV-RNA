#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=s6g_S8
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6G_GGAC_DRACH/logs/%x_%j.out
# Step 6G DRACH power check, stage 3: dose table vs real folds.
set -eo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6G_GGAC_DRACH; cd $B; source inputs/paths.sh
(cd scripts && md5sum -c --quiet SCRIPTS_r3.md5)
$PY scripts/g9_power_table6g.py run/G8_power_sets run/G10_power_perm_GC run/G10_power_perm_noGC run/G7_genes_GC run/G7_genes_noGC run/G9_power_table
for d in G10_power_perm_GC G10_power_perm_noGC; do (cd run/$d && md5sum * > ${d}_outputs.md5); chmod -R a-w run/$d; done
(cd run/G9_power_table && md5sum * > G9_outputs.md5); chmod -R a-w run/G9_power_table; echo S8_DONE
