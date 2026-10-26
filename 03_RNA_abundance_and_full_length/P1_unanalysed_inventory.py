# P1 r2 (read-only; r2 = lncRNA units matched by gene id, V3 also matched by 5-digit locus_tag)
# P1 (read-only): which final.gtf transcripts have NOT been analysed (not in the PCG/lncRNA/tRNA full-length units),
# and what FungiDB / NCBI V3 call them. Descriptive only; no decisions.
import sys, gzip, re, collections
P="/scratch/flljus002/Phase3D_forensic_audit"; D=P+"/v2_rRNA_2026-10-09"; O=D+"/P1_unanalysed_inventory"
GTF="/scratch/flljus002/CandidaEV_Project/_Archive_2026-08-31/01_Differential_expression_analysis/01_Reference_files/final.gtf"
FDB="/scratch/flljus002/CandidaEV_Project/05_EV_coverage_analysis/11_Full_length_biological_annotation/output/fungidb_gene_annotation_source.gff"
V3=D+"/refs/raw/genome_annotations/GCA_002759435.3_Cand_auris_B8441_V3_genomic.gff.gz"
UNIT=P+"/v2_tRNA_2026-10-09/run/FL_r2/unit_table.tsv.gz"
INTV=P+"/v2_read_census_2026-10-09/inputs/intervals.tsv"
TRNA=P+"/v2_tRNA_2026-10-09/inputs/candidate_trna_reference_168.bed"
LRD=P+"/v2_bin150/run/R1_annotation/inputs/lncRNA_rDNA_provenance.tsv"
SRD=P+"/v2_bin150/run/R1_annotation/inputs/StringTie_rDNA_provenance.tsv"
# transcripts
tx={}
for l in open(GTF):
    f=l.rstrip("\n").split("\t")
    if len(f)<9 or f[2]!="transcript": continue
    t=re.search(r"transcript_id \"([^\"]+)\"",f[8]).group(1); g=re.search(r"gene_id \"([^\"]+)\"",f[8]).group(1)
    tx[t]=dict(t=t,g=g,src=f[1],chrom=f[0],s=int(f[3]),e=int(f[4]),strand=f[6])
analysed=set()
with gzip.open(UNIT,"rt") as h:
    next(h)
    for l in h: analysed.add(l.split("\t")[1])
cls=collections.defaultdict(set)
for l in open(INTV):
    f=l.rstrip("\n").split("\t")
    if f[0]=="chrom": continue
    cls[f[5]].add(f[3])
# FungiDB descriptions (gene level)
fdesc={}
for l in open(FDB):
    f=l.rstrip("\n").split("\t")
    if len(f)<9: continue
    m=re.search(r"ID=([^;]+)",f[8]); 
    if not m: continue
    gid=m.group(1).split("-t")[0]
    d=re.search(r"description=([^;]+)",f[8]); b=re.search(r"ebi_biotype=([^;]+)",f[8])
    if gid not in fdesc or f[2].endswith("gene"):
        fdesc[gid]=( (d.group(1) if d else ""), (b.group(1) if b else f[2]) )
# V3 products keyed by old v2 ID (Note=B8441v2:<id>T0)
v3={}
with gzip.open(V3,"rt") as h:
    for l in h:
        f=l.rstrip("\n").split("\t")
        if len(f)<9 or f[2] in ("gene","exon","CDS","region"): continue
        m=re.search(r"B8441v2:([A-Za-z0-9_]+?)T0",f[8]); p=re.search(r"product=([^;]+)",f[8]); lt=re.search(r"locus_tag=([^;]+)",f[8])
        val=(f[2],p.group(1) if p else "",lt.group(1) if lt else "")
        if m: v3[m.group(1)]=val
        if lt: v3.setdefault(lt.group(1),val)
# tRNA bed and rDNA flags
trna=[l.split("\t") for l in open(TRNA) if not l.startswith("#")]
rdna={}
for fn in (LRD,SRD):
    for i,l in enumerate(open(fn)):
        f=l.rstrip("\n").split("\t")
        if i==0: continue
        rdna[f[0]]=f[-1]
def trna_hit(r):
    for b in trna:
        if b[0]==r["chrom"] and int(b[1])<r["e"] and int(b[2])>r["s"]-1: return b[3]
    return ""
rows=[]
for t,r in tx.items():
    if t in analysed or r["g"] in analysed: continue
    g=r["g"]; fd=fdesc.get(g,("","")); v=v3.get(g,("","",""))
    rows.append([t,g,r["src"],r["chrom"],r["s"]-1,r["e"],r["strand"],r["e"]-r["s"]+1,";".join(sorted(cls.get(t,[]))) or "none",
                 fd[1],fd[0],v[0],urllib_unq(v[1]) if False else v[1].replace("%2C",",").replace("%20"," "),v[2],trna_hit(r),rdna.get(g,"")])
rows.sort(key=lambda x:(x[2],x[8],x[3],x[4]))
with open(O+"/unanalysed_transcripts.tsv","w") as w:
    w.write("transcript\tgene\tgtf_source\tchrom\tstart0\tend\tstrand\tlength\tcensus_interval_class\tfungidb_biotype\tfungidb_description\tV3_feature\tV3_product\tV3_locus_tag\toverlaps_tRNAscan168\trDNA_flag\n")
    for x in rows: w.write("\t".join(map(str,x))+"\n")
c=collections.Counter((x[2],x[8],x[9] or x[11] or "-") for x in rows)
with open(O+"/summary.tsv","w") as w:
    w.write(f"# transcripts in final.gtf: {len(tx)}; analysed (FL units): {len([t for t in tx if t in analysed or tx[t]["g"] in analysed])}; not analysed: {len(rows)}\n")
    w.write("gtf_source\tcensus_class\tfungidb_or_V3_type\tn\n")
    for k,n in sorted(c.items()): w.write("\t".join(k)+f"\t{n}\n")
