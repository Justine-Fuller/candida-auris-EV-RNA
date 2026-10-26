#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=s6f_F2
#SBATCH --array=1-24%8
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6F_whole_RNA_tag/logs/%x_%A_%a.out
# Step 6F F2: reads per gene model per library (12 EV + 12 Cell sRNA), Step 6E read rule. manifest_6E.tsv HAS a header.
# Gate: accounting must equal the 6E K1 accounting of the same library (same rule) for every shared key.
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6F_whole_RNA_tag; P=/scratch/flljus002/Phase3D_forensic_audit; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3
source $B/inputs/paths.sh
read -r LIB COMP GROUP BAM < <(awk -F'\t' -v i=$((SLURM_ARRAY_TASK_ID+1)) 'NR==i' $MAN6E)
[ -n "$LIB" ] || { echo "STOP: empty row"; exit 1; }
D=$B/run/F2_counts; mkdir -p $D; [ ! -e $D/$LIB.gene_counts.tsv ] || { echo "STOP: exists"; exit 1; }
grep "  $BAM\$" $P/v2_protein_coding/provenance/J00/$LIB.md5 | md5sum -c -
md5sum -c --quiet $B/inputs/inputs.md5; (cd $B/scripts && md5sum -c --quiet SCRIPTS_r1.md5)
$PY $B/scripts/s6f_gene_counts.py $BAM $LIB $A100 $N100 $D
$PY - $D/$LIB.accounting.json $K1/$LIB.accounting.json $D/$LIB.gene_counts.tsv > $D/$LIB.gate.tsv <<'PYEOF'
import sys, json, csv
a = json.load(open(sys.argv[1]))["accounting"]; k = json.load(open(sys.argv[2]))
keys = ["records", "primary_mapped", "eligible_NH1_ge18", "model_none", "model_ambiguous", "PCG_CDS_reads", "lncRNA_reads"]
ok = True
for key in keys + ["model_competitor_only", "model_no_SEQ"]:
    x, y = a.get(key, 0), k.get(key, 0)
    if key == "model_no_SEQ": good = y == 0          # 6E dropped reads without SEQ; none expected
    else: good = x == y
    ok &= good; print(key, x, y, "PASS" if good else "FAIL", sep="\t")
c = list(csv.DictReader(open(sys.argv[3]), delimiter="\t"))
for cl in ("PCG_CDS", "lncRNA"):
    s = sum(int(r["reads"]) for r in c if r["rna_class"] == cl); good = s == a.get(cl + "_reads", 0); ok &= good
    print("sum_" + cl, s, a.get(cl + "_reads", 0), "PASS" if good else "FAIL", sep="\t")
print("F2_GATE", "PASS" if ok else "FAIL", sep="\t")
PYEOF
cat $D/$LIB.gate.tsv
chmod a-w $D/$LIB.gene_counts.tsv $D/$LIB.accounting.json $D/$LIB.gate.tsv
grep -q $'F2_GATE\tPASS' $D/$LIB.gate.tsv && echo F2_LIB_DONE_PASS || { echo F2_LIB_GATE_FAIL; exit 1; }
