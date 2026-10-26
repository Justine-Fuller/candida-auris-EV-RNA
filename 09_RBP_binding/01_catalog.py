# Catalogue hg19 eCLIP bigBeds: RBP, cell line, replicate, ENCODE accession, upload date
import pyBigWig, os, csv, collections
rows=[]
for k,name in (l.rstrip('\n').split('\t') for l in open('./rbp/eclip_hg19_files.tsv')):
    acc=os.path.basename(k).split('.')[0]
    bb=pyBigWig.open(f'./rbp/eclip_bb/{acc}.bigBed')
    assert bb.chroms('chr1')==249250621
    e=bb.entries('chr10',0,135534747)
    f=e[0][2].split('\t') if e else []
    bb.close()
    rbp,cell,rep=name.split('_')[0],name.split('_')[1],name.split('_')[2][:5]
    rows.append(dict(acc=acc,date=k[:10],rbp=rbp,cell=cell,rep=rep,label=name,nfields=len(f),example='|'.join(f)))
# keep the most recent file per RBP/cell/rep (reprocessed releases supersede older ones)
best={}
for r in sorted(rows,key=lambda r:r['date']):
    best[(r['rbp'],r['cell'],r['rep'])]=r
for r in rows: r['used']=best[(r['rbp'],r['cell'],r['rep'])]['acc']==r['acc']
with open('./rbp/out/eclip_file_catalogue.tsv','w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(sorted(rows,key=lambda r:(r['rbp'],r['cell'],r['rep'],r['date'])))
u=[r for r in rows if r['used']]
print('files',len(rows),'used',len(u),'RBPs',len({r['rbp'] for r in u}),'RBP-cell',len({(r['rbp'],r['cell']) for r in u}))
print(collections.Counter(r['nfields'] for r in rows)); print(rows[0]['example'])
