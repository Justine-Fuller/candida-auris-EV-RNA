#!/bin/bash
# P0: read-only diagnostic. Do the downloaded curated rRNA sequences match the B8441 rDNA? Not the annotation itself.
set -euo pipefail
D=/scratch/flljus002/Phase3D_forensic_audit/v2_rRNA_2026-10-09; B=/scratch/flljus002/.conda/envs/orthofinder/bin
A=/scratch/flljus002/Phase3D_forensic_audit/lncRNA1946_annotation_audit/results
O=$D/P0_reference_match; mkdir -p $O
python3 - "$D/refs/raw/PDB_9SPJ.fasta" "$O/PDB_9SPJ_rRNA.fa" <<"PY"
import sys
out=open(sys.argv[2],"w"); keep=False
for l in open(sys.argv[1]):
    if l.startswith(">"):
        keep=("rRNA" in l or "5.8S RNA" in l)
        if keep: out.write(">PDB9SPJ_"+l.split("|")[2].replace(" ","_")+"\n")
    elif keep: out.write(l.strip().replace("U","T")+"\n")
PY
cat $O/PDB_9SPJ_rRNA.fa $D/refs/raw/Cauris_NG_073563.1.fasta $D/refs/raw/Calbicans_SC5314_rDNA_RefSeq.fasta $D/refs/raw/Scerevisiae_S288C_rDNA_RefSeq.fasta > $O/queries.fa
cat $A/region_Chr2_2280000_2400000.fa $A/region_Chr3_1_20000.fa $A/region_Chr5_1_20000.fa > $O/subject.fa
FMT="6 qseqid sseqid qstart qend qlen sstart send pident length mismatch gapopen evalue bitscore sstrand"
$B/blastn -task blastn -dust no -query $O/queries.fa -subject $O/subject.fa -evalue 1e-10 -outfmt "$FMT" > $O/blast_all.tsv
$B/blastn -version | head -1 > $O/blast_version.txt
python3 - $O <<"PY"
import sys,collections
o=sys.argv[1]; best=collections.defaultdict(list)
for l in open(o+"/blast_all.tsv"):
    f=l.rstrip().split("\t"); best[(f[0],f[1])].append(f)
with open(o+"/summary.tsv","w") as w:
    w.write("query\tqlen\tsubject\tn_hits_covering>=90%qlen\tbest_pident\tbest_qcov\tbest_qstart\tbest_qend\tbest_strand\n")
    for (q,s),hs in sorted(best.items()):
        ql=int(hs[0][4]); full=[h for h in hs if (int(h[3])-int(h[2])+1)>=0.9*ql]
        b=max(hs,key=lambda h:float(h[12]))
        w.write(f"{q}\t{ql}\t{s}\t{len(full)}\t{b[7]}\t{(int(b[3])-int(b[2])+1)/ql:.3f}\t{b[2]}\t{b[3]}\t{b[13]}\n")
PY
md5sum $O/* > $D/provenance/P0_outputs_md5.txt
