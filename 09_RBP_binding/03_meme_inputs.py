# Step 3a: for each RBP with a reproducible peak on DDX21/SMC1A/RPS3, collect its top genome-wide bound sequences (hg19)
import pyBigWig, csv, collections, subprocess, os
OUT='./rbp/out/'; TOPN=2000; MAXLEN=100   # 2000 ranked candidates; 03b drops repeats and keeps the top 500
P=list(csv.DictReader(open(OUT+'gene_eclip_peaks_hg19.tsv'),delimiter='\t'))
rbps=sorted({p['rbp'] for p in P if p['reproducible']=='1'})
cat=[r for r in csv.DictReader(open(OUT+'eclip_file_catalogue.tsv'),delimiter='\t') if r['used']=='True' and r['rbp'] in rbps]
os.makedirs(OUT+'meme_input',exist_ok=True)
chroms=[f'chr{i}' for i in range(1,23)]+['chrX']
log=[]
for rbp in rbps:
    peaks=[]
    for cell in sorted({r['cell'] for r in cat if r['rbp']==rbp}):
        reps={}
        for r in cat:
            if r['rbp']==rbp and r['cell']==cell:
                bb=pyBigWig.open(f"./rbp/eclip_bb/{r['acc']}.bigBed"); L=[]
                for c in chroms:
                    for s,e,rest in bb.entries(c,0,bb.chroms(c)) or []:
                        f=rest.split('\t'); fc,p=float(f[3]),float(f[4])
                        if fc>=2 and p>=2: L.append((c,s,e,f[2],fc))
                bb.close(); reps[r['rep']]=L
        if len(reps)<2: continue
        # reproducible: overlaps a same-strand peak in the other replicate
        idx=collections.defaultdict(list)
        for c,s,e,st,fc in reps['rep02']: idx[(c,st,s//1000)].append((s,e))
        for c,s,e,st,fc in reps['rep01']:
            if any(a<e and s<b for k in (s//1000-1,s//1000,e//1000) for a,b in idx[(c,st,k)]): peaks.append((fc,c,s,e,st,cell))
    peaks.sort(reverse=True)
    seen=set(); top=[]
    for fc,c,s,e,st,cell in peaks:
        if st=='+': e=min(e,s+MAXLEN)
        else: s=max(s,e-MAXLEN)        # keep the 5' (crosslink) end of the peak
        key=(c,st,s//50)
        if key in seen: continue     # avoid the same site twice (two cell lines)
        if e-s<30: m=(s+e)//2; s,e=m-15,m+15   # pad very short peaks to 30 nt
        seen.add(key); top.append((c,s,e,f'{rbp}|{cell}|{c}:{s}-{e}({st})|{fc:.2f}',0,st))
        if len(top)==TOPN: break
    with open(OUT+f'meme_input/{rbp}.bed','w') as f:
        for t in top: f.write('\t'.join(map(str,t))+'\n')
    log.append((rbp,len(peaks),len(top)))
    print(rbp,len(peaks),len(top),flush=True)
with open(OUT+'meme_input/summary.tsv','w') as f:
    f.write('rbp\treproducible_peaks_genomewide\tsequences_used\n')
    for l in log: f.write('\t'.join(map(str,l))+'\n')
