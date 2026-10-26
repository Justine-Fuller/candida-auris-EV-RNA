"""Make the whole-proteome input FASTA with transcript IDs (header transcript= field), as in up/down.fasta,
and check that every up/down protein is identical to the proteome copy."""
import sys, hashlib
src, out, up, down, report = sys.argv[1:6]
def read(p):
    seqs, order, cur = {}, [], None
    for line in open(p):
        line = line.rstrip("\n")
        if line.startswith(">"):
            f = [x.strip() for x in line[1:].split("|")]
            tr = [x[11:] for x in f if x.startswith("transcript=")]
            cur = tr[0] if tr else f[0]   # use transcript ID, as in up/down.fasta
            assert cur not in seqs, cur
            order.append(cur); seqs[cur] = []
        elif line:
            seqs[cur].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}, order
prot, order = read(src)
assert len(prot) == len(order), "duplicate IDs"
with open(out, "w") as fh:
    for k in order:
        s = prot[k]
        fh.write(f">{k}\n")
        for i in range(0, len(s), 60): fh.write(s[i:i+60] + "\n")
lines = [f"proteome {src}: {len(order)} proteins",
         f"proteins containing '*': {sum('*' in s for s in prot.values())}",
         f"proteins with non-standard letters (not ACDEFGHIKLMNPQRSTVWY*): {sum(bool(set(s)-set('ACDEFGHIKLMNPQRSTVWY*')) for s in prot.values())}"]
for name, p in (("up", up), ("down", down)):
    d, _ = read(p)
    missing = [k for k in d if k not in prot]
    diff = [k for k in d if k in prot and d[k] != prot[k]]
    lines.append(f"{name}.fasta: {len(d)} proteins; not in proteome {len(missing)}; sequence differs {len(diff)} {diff[:5]}")
open(report, "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
