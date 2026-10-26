#!/bin/bash
# Step 4: collect MEME motifs (E < 0.05), scan EV-enriched and comparison C. auris RNAs with FIMO (given strand only)
set -euo pipefail
E=./mm/env/bin; O=./rbp/out; EV=./rbp/ev
cd $O
# keep motifs with MEME E-value < 0.05, renamed RBP_motifN
python3 - <<'PY'
import glob,re,os
keep=[]
for d in sorted(glob.glob('meme_out/*/meme.txt')):
    rbp=d.split('/')[1]
    for m in re.finditer(r'^MOTIF (\S+) MEME-(\d+)\s+width =\s+(\d+)\s+sites =\s+(\d+)\s+llr =\s+\d+\s+E-value =\s+(\S+)',open(d).read(),re.M):
        if float(m.group(5))<0.05: keep.append((rbp,m.group(2),m.group(1),m.group(5)))
open('motifs_kept.tsv','w').write('rbp\tmeme_n\tconsensus\tE\n'+''.join('\t'.join(k)+'\n' for k in keep))
print(len(keep),'motifs from',len({k[0] for k in keep}),'RBPs')
PY
rm -f all_motifs.meme
for d in meme_out/*/; do r=$(basename $d); n=$(awk -F'\t' -v r=$r '$1==r{printf "%s,",$2}' motifs_kept.tsv); [ -z "$n" ] && continue
  for i in ${n//,/ }; do $E/meme-get-motif -id $(awk -F'\t' -v r=$r -v i=$i '$1==r && $2==i{print $3}' motifs_kept.tsv) $d/meme.txt | sed "s/^MOTIF \(\S*\) MEME-$i/MOTIF ${r}_m$i \1/" > tmp_${r}_${i}.meme
    if [ ! -s all_motifs.meme ]; then cat tmp_${r}_${i}.meme > all_motifs.meme; else sed -n '/^MOTIF/,$p' tmp_${r}_${i}.meme >> all_motifs.meme; fi; rm tmp_${r}_${i}.meme; done; done
grep -c ^MOTIF all_motifs.meme
# C. auris background (0-order) from all scanned RNAs
cat $EV/ev_enriched.fa $EV/comparison.fa > $EV/all_scanned.fa
$E/fasta-get-markov -m 0 -norc $EV/all_scanned.fa $EV/cauris_bg0.txt >/dev/null 2>&1
# dinucleotide-shuffled EV RNAs = "how many hits by chance"
$E/fasta-shuffle-letters -kmer 2 -seed 1 -tag _shuf $EV/ev_enriched.fa > $EV/ev_enriched_shuffled.fa
for s in ev_enriched comparison ev_enriched_shuffled; do
  $E/fimo --norc --bfile $EV/cauris_bg0.txt --thresh 1e-4 --max-stored-scores 10000000 --text all_motifs.meme $EV/$s.fa > fimo_$s.tsv 2>/dev/null
  echo $s $(($(wc -l < fimo_$s.tsv)-1)) hits
done
