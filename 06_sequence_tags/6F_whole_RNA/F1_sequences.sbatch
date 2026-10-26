#!/bin/bash
#SBATCH --partition=ada
#SBATCH --job-name=s6f_F1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6F_whole_RNA_tag/logs/%x_%j.out
# Step 6F F1: per-gene sequences (UTR5 / CDS / UTR3 / whole; lncRNA whole) + independent samtools faidx check of 3 models.
set -euo pipefail
B=/scratch/flljus002/Phase3D_forensic_audit/v2_downstream_2026-10-12/step6F_whole_RNA_tag; PY=/opt/exp_soft/miniconda3-py3.12/bin/python3
source $B/inputs/paths.sh
md5sum -c --quiet $B/inputs/inputs.md5; (cd $B/scripts && md5sum -c --quiet SCRIPTS_r1.md5)
D=$B/run/F1_sequences; [ ! -e $D ] || { echo "STOP: $D exists"; exit 1; }; mkdir -p $D
$PY $B/scripts/s6f_sequences.py $A100 $N100 $GTF $GENOME $D
module load software/samtools-1.20
# independent check: rebuild "whole" for the first + strand PCG with both UTRs, first - strand PCG with both UTRs, first
# multi-exon lncRNA, from GTF UTR rows + bin blocks via samtools faidx (-i = reverse complement)
$PY - $D/gene_sequences.tsv.gz $GTF $A100 > $D/handcheck_input.tsv <<'PYEOF'
import sys, gzip, csv, re
seqs = list(csv.DictReader(gzip.open(sys.argv[1], "rt"), delimiter="\t"))
pick = []
for st in "+-":
    pick.append(next(r for r in seqs if r["rna_class"] == "PCG_CDS" and r["strand"] == st and int(r["len_UTR5"]) > 0 and int(r["len_UTR3"]) > 0))
blk = {}
for r in csv.DictReader(open(sys.argv[3]), delimiter="\t"):
    blk.setdefault(r["model_id"], []).extend(r["genomic_blocks_5to3"].split(";"))
def spliced(m):
    iv = sorted(tuple(map(int, b.split("-"))) for b in blk[m])
    return any(b[0] != a[1] for a, b in zip(iv, iv[1:]))
lnc = [r for r in seqs if r["rna_class"] == "lncRNA"]
# r2: a spliced lncRNA if one exists (r1 crashed: none found by its test), else the first - strand lncRNA
pick.append(next((r for r in lnc if spliced(r["model_id"])), None) or next(r for r in lnc if r["strand"] == "-"))
sys.stderr.write("spliced lncRNA models: %d\n" % sum(spliced(r["model_id"]) for r in lnc))
utr = {}
for l in open(sys.argv[2]):
    x = l.split("\t")
    if len(x) > 8 and x[2] in ("5_prime_UTR", "3_prime_UTR"):
        t = re.search(r'transcript_id "([^"]+)"', x[8]).group(1); utr.setdefault(t, []).append(f"{int(x[3])-1}-{x[4]}")
for r in pick:
    ivs = blk[r["model_id"]] + utr.get(r["model_id"], [])
    print(r["model_id"], r["chrom"], r["strand"], ";".join(ivs), r["seq_whole"], sep="\t")
PYEOF
: > $D/handcheck_sequences.tsv
while IFS=$'\t' read -r mid chrom strand ivs sw; do
  s=""; for iv in $(echo $ivs | tr ';' '\n' | sort -t- -k1,1n); do a=${iv%-*}; b=${iv#*-}
    s="$s$(samtools faidx $GENOME "$chrom:$((a+1))-$b" | grep -v '>' | tr -d '\n' || true)"; done
  # merge adjacent bins: samtools pieces are concatenated in coordinate order; reverse complement for - strand
  [ "$strand" = "-" ] && s=$(echo "$s" | rev | tr ACGTacgt TGCAtgca)
  s=$(echo "$s" | tr a-z A-Z)
  echo -e "$mid\t$strand\t${#s}\t${#sw}\t$([ "$s" = "$sw" ] && echo PASS || echo FAIL)" >> $D/handcheck_sequences.tsv
done < $D/handcheck_input.tsv
cat $D/handcheck_sequences.tsv
[ $(grep -c $'\tPASS$' $D/handcheck_sequences.tsv || true) -eq 3 ] && echo HANDCHECK_PASS || echo HANDCHECK_FAIL
(cd $D && md5sum gene_sequences.tsv.gz checks.tsv handcheck_input.tsv handcheck_sequences.tsv > F1_outputs.md5); chmod -R a-w $D
