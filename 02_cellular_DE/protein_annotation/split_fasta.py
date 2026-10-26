"""Split the proteome FASTA into N chunks of consecutive records (biolib DeepTMHMM rejects ~2.8 MB inputs;
the old up/down runs were 0.8 / 1.3 MB)."""
import sys
src, outprefix, n = sys.argv[1], sys.argv[2], int(sys.argv[3])
recs, cur = [], None
for line in open(src):
    if line.startswith(">"): cur = [line]; recs.append(cur)
    else: cur.append(line)
size = -(-len(recs) // n)
for i in range(n):
    part = recs[i*size:(i+1)*size]
    with open(f"{outprefix}{i+1}.fasta", "w") as fh:
        for r in part: fh.writelines(r)
    print(f"chunk {i+1}: {len(part)} proteins")
