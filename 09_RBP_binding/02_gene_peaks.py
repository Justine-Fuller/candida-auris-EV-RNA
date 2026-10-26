# Step 1-2: eCLIP peaks (hg19) on DDX21, SMC1A, RPS3, filtered as in LncLOOM (same strand, log2FC>=2, -log10p>=2)
import pyBigWig, csv, collections
OUT='./rbp/out/'
GENES={'DDX21':'NM_004728','SMC1A':'NM_006306','RPS3':'NM_001005'}
tx={}
for l in open('./rbp/ref/refGene.txt'):
    f=l.rstrip('\n').split('\t')
    if f[12] in GENES and f[1]==GENES[f[12]] and '_' not in f[2]:
        tx[f[12]]=dict(chrom=f[2],strand=f[3],txs=int(f[4]),txe=int(f[5]),cds=int(f[6]),cde=int(f[7]),
            es=[int(x) for x in f[9].strip(',').split(',')],ee=[int(x) for x in f[10].strip(',').split(',')])
# LeafCutter junctions (GRCh38) -> hg19 by a per-gene constant shift that puts junction ends on RefSeq exon edges
lc=collections.defaultdict(list)
for l in open('/mnt/user-data/uploads/CandidaEV/Differential splicing analysis/01_LeafCutter_ds_caur/03_Output/RPMI_vs_caur_live_effect_sizes.txt'):
    f=l.split('\t')
    if f[0].count(':')!=3: continue
    c,s,e,clu=f[0].split(':')
    g={'clu_149_+':'DDX21','clu_1183_-':'SMC1A','clu_218_+':'RPS3'}.get(clu)
    if g: lc[g].append((int(s),int(e),float(f[4])))
junc=[]
for g,J in lc.items():
    t=tx[g]; edges=set(t['ee'])|set(t['es'])
    cand=collections.Counter()
    for s,e,_ in J:
        for x in edges:
            for y in (s,e):
                for k in (-1,0,1): cand[x-y+k]+=0  # register
    best=max(cand,key=lambda d:sum(((s+d) in edges or (s+d-1) in edges)+((e+d) in edges or (e+d-1) in edges) for s,e,_ in J))
    for s,e,dpsi in J: junc.append(dict(gene=g,chrom=t['chrom'],intron_start_hg19=s+best+1,intron_end_hg19=e+best-1,dPSI=dpsi,shift=best))
with open(OUT+'leafcutter_junctions_hg19.tsv','w') as f:
    w=csv.DictWriter(f,fieldnames=list(junc[0]),delimiter='\t'); w.writeheader(); w.writerows(junc)
def region(g,p):
    t=tx[g]
    for s,e in zip(t['es'],t['ee']):
        if s<=p<e:
            if p<t['cds']: return '5UTR' if t['strand']=='+' else '3UTR'
            if p>=t['cde']: return '3UTR' if t['strand']=='+' else '5UTR'
            return 'CDS'
    return 'intron'
cat=[r for r in csv.DictReader(open(OUT+'eclip_file_catalogue.tsv'),delimiter='\t') if r['used']=='True']
peaks=[]
for r in cat:
    bb=pyBigWig.open(f"./rbp/eclip_bb/{r['acc']}.bigBed")
    for g,t in tx.items():
        for s,e,rest in bb.entries(t['chrom'],t['txs'],t['txe']) or []:
            name,score,strand,fc,p=rest.split('\t')[:5]
            if strand!=t['strand']: continue
            fc,p=float(fc),float(p)
            if not (fc>=2 and p>=2): continue
            mid=(s+e)//2
            ss=[(j['intron_start_hg19'],j['intron_end_hg19'],j['dPSI']) for j in junc if j['gene']==g]
            dist=min(min(abs(mid-a),abs(mid-b)) for a,b,_ in ss) if ss else ''
            peaks.append(dict(gene=g,chrom=t['chrom'],start=s,end=e,strand=strand,rbp=r['rbp'],cell=r['cell'],rep=r['rep'],
                log2FC=round(fc,2),neglog10p=round(p,2),region=region(g,mid),dist_to_LeafCutter_splice_site=dist,
                encode_strict=int(fc>=3 and p>=3),acc=r['acc']))
    bb.close()
# reproducible = overlapping peak for the same RBP + cell line in the other replicate
idx=collections.defaultdict(list)
for p in peaks: idx[(p['rbp'],p['cell'],p['rep'],p['gene'])].append(p)
for p in peaks:
    other='rep02' if p['rep']=='rep01' else 'rep01'
    p['reproducible']=int(any(q['start']<p['end'] and p['start']<q['end'] for q in idx[(p['rbp'],p['cell'],other,p['gene'])]))
peaks.sort(key=lambda p:(p['gene'],p['start'],p['rbp']))
with open(OUT+'gene_eclip_peaks_hg19.tsv','w') as f:
    w=csv.DictWriter(f,fieldnames=list(peaks[0]),delimiter='\t'); w.writeheader(); w.writerows(peaks)
# per gene x RBP summary
S=collections.defaultdict(lambda:collections.defaultdict(list))
for p in peaks: S[(p['gene'],p['rbp'])]['p'].append(p)
rows=[]
for (g,rbp),d in S.items():
    P=d['p']; R=[p for p in P if p['reproducible']]
    rows.append(dict(gene=g,rbp=rbp,cells=','.join(sorted({p['cell'] for p in P})),n_peaks=len(P),n_reproducible=len(R),
        n_encode_strict=sum(p['encode_strict'] for p in P),max_log2FC=max(p['log2FC'] for p in P),
        regions=','.join(f'{k}:{v}' for k,v in collections.Counter(p['region'] for p in P).most_common()),
        min_dist_to_splice_site=min((p['dist_to_LeafCutter_splice_site'] for p in P if p['dist_to_LeafCutter_splice_site']!=''),default='')))
rows.sort(key=lambda r:(r['gene'],-r['n_reproducible'],-r['max_log2FC']))
with open(OUT+'gene_rbp_summary.tsv','w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)
for g in GENES:
    G=[r for r in rows if r['gene']==g]
    print(g,tx[g]['chrom'],tx[g]['txs'],tx[g]['txe'],tx[g]['strand'],'RBPs',len(G),'with reproducible peak',sum(r['n_reproducible']>0 for r in G))
print({g:(J[0]['shift'] if (J:=[j for j in junc if j['gene']==g]) else None) for g in GENES})
