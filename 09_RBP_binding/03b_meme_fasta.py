# Step 3b: sequences for MEME. Drop peaks that are mostly repeat (>50% soft-masked in hg19), keep the top 500 per RBP
import glob, os, subprocess
O='./rbp/out/'; os.makedirs(O+'meme_fa',exist_ok=True); log=[]
for b in sorted(glob.glob(O+'meme_input/*.bed')):
    r=os.path.basename(b)[:-4]
    fa=subprocess.run(['./mm/env/bin/bedtools','getfasta','-s','-name','-tab','-fi','./rbp/ref/hg19.fa','-bed',b],capture_output=True,text=True).stdout.splitlines()
    keep=[]; rep=0
    for l in fa:                       # bed is already ranked by log2FC
        n,s=l.split('\t')
        if sum(c.islower() for c in s)/len(s)>0.5: rep+=1; continue
        keep.append((n,s.upper()))
        if len(keep)==500: break
    with open(O+f'meme_fa/{r}.fa','w') as f:
        for n,s in keep: f.write(f'>{n}\n{s}\n')
    log.append(f'{r}\t{len(fa)}\t{rep}\t{len(keep)}')
open(O+'meme_fa/summary.tsv','w').write('rbp\tcandidates\trepeat_dropped_before_top500\tsequences_used\n'+'\n'.join(log)+'\n')
