#!/usr/bin/env python3
"""trna_cca_recount.py -- v2 tRNA branch, CCA recount (2026-10-13; Justine chose option b after the no-CCA check, see
../cca_check_2026-10-13/README.md). READ-ONLY. Does NOT replace trna_count.py / T1 outputs (kept as provenance).
Usage: trna_cca_recount.py <bam> <library> <trna168.bed> <genome.fasta> <outdir> [sRNA|polyA]
Family assignment and the OLD CCA rule are copied unchanged from trna_count.py r2 (reads >= 18 nt; records grouped per QNAME
until NH; a read belongs to family F if all its same-strand tRNA hits are in F; reads with an N gap > 1000 are not checked):
  old: an alignment ends EXACTLY at the locus 3' end -> reach; its 3' soft clip -> clipped; clip starting with CCA -> CCA.
NEW rule (2026-10-13): for every alignment of the read that hits F and every F locus it hits, if the alignment covers the
locus' last mature base (the 3' end base, aligned as M/=/X), the TAIL = all read bases after that base in RNA orientation,
whether aligned to the genome past the locus end or soft-clipped. The read reaches the 3' end (reach_new) if any such
alignment exists; its tail class is the best over those alignments: CCA (tail starts with CCA) > CC > C > other > none (empty).
CCA reads are split: templated_ambiguous = in EVERY alignment giving a CCA tail, the 3 CCA bases are aligned to genome
bases that are CCA themselves (the tail could be genomic trailer); otherwise non_templated.
Invariant (gate): for reads with no aligned base past the locus end in any F alignment, new CCA == old CCA.
Output: <lib>.cca2.json (per family: old fields as in T1 + new fields)."""
import sys, json, collections, pysam
BAM, LIB, TB, FA, OUT = sys.argv[1:6]; KIND = sys.argv[6] if len(sys.argv) > 6 else "sRNA"
assert KIND in ("sRNA", "polyA")
B = 10000; RC = str.maketrans("ACGTN", "TGCAN"); rc = lambda s: s.translate(RC)[::-1]
G = {}; nm = None; buf = []
for l in open(FA):
    if l.startswith(">"):
        if nm: G[nm] = "".join(buf)
        nm = l[1:].split()[0]; buf = []
    else: buf.append(l.strip().upper())
G[nm] = "".join(buf)
tidx = collections.defaultdict(list); loci = []
with open(TB) as f:
    next(f)
    for line in f:
        ch, s, e, name, sc, st, fam, fid, intron, tier = line.rstrip("\n").split("\t"); s, e = int(s), int(e)
        k = len(loci); loci.append(dict(ch=ch, s=s, e=e, st=st, fam=fam, name=name))
        for b in range(s // B, (e - 1) // B + 1): tidx[(ch, st, b)].append(k)
def rstrand(a):
    s = "-" if a.is_reverse else "+"
    return ({"+": "-", "-": "+"}[s]) if KIND == "polyA" else s
def trna_hits(a):
    hits = set()
    for s, e in a.get_blocks():
        for b in range(s // B, (e - 1) // B + 1):
            for k in tidx.get((a.reference_name, rstrand(a), b), ()):
                L = loci[k]
                if L["s"] < e and s < L["e"]: hits.add(k)
    return hits
RANK = {"CCA": 5, "CC": 4, "C": 3, "other": 2, "none": 1}
def tclass(t): return "none" if not t else "CCA" if t.startswith("CCA") else "CC" if t.startswith("CC") else "C" if t.startswith("C") else "other"
def new_tail(a, seq, L):
    """(tail, templated, aligned_past_end) for alignment a on locus L, or None if a does not cover the last mature base."""
    pairs = a.get_aligned_pairs(); ch = G[a.reference_name]
    if L["st"] == "+":
        qp = next((q for q, r in pairs if r == L["e"] - 1 and q is not None), None)
        if qp is None: return None
        tail = seq[qp + 1:]; q2r = {q: r for q, r in pairs if q is not None}
        tpos = [q2r.get(qp + i) for i in (1, 2, 3)]
        tmpl = len(tail) >= 3 and all(r is not None and ch[r] == seq[qp + i] for i, r in zip((1, 2, 3), tpos))
        past = a.reference_end > L["e"]
    else:
        qp = next((q for q, r in pairs if r == L["s"] and q is not None), None)
        if qp is None: return None
        tail = rc(seq[:qp]); q2r = {q: r for q, r in pairs if q is not None}
        tpos = [q2r.get(qp - i) for i in (1, 2, 3)]
        tmpl = len(tail) >= 3 and all(r is not None and ch[r] == seq[qp - i] for i, r in zip((1, 2, 3), tpos))
        past = a.reference_start < L["s"]
    return tail, tmpl, past
acc = collections.Counter(); fam = collections.defaultdict(collections.Counter); pending = {}
def flush(rs):
    L_ = rs[0].infer_read_length() or rs[0].query_length
    if L_ < 18: return
    hits = {i: trna_hits(a) for i, a in enumerate(rs)}
    fams = {loci[k]["fam"] for hs in hits.values() for k in hs}
    if not fams: return
    acc["trna_reads"] += 1
    if len(fams) > 1: acc["trna_reads_cross_family_excluded"] += 1; return
    F = fams.pop(); fam[F]["reads"] += 1
    if rs[0].get_tag("NH") == 1: fam[F]["reads_NH1"] += 1
    if any(op == 3 and ln > 1000 for a in rs for op, ln in a.cigartuples):
        fam[F]["reads_bigN_excluded_from_CCA_check"] += 1; return
    reach = clip = cca = False; best = None; tmpl_all = True; past_any = False; reach_new = False
    for i, a in enumerate(rs):
        if not hits[i]: continue
        ct = a.cigartuples; seq = a.query_sequence or ""
        if not seq:
            p0 = next((x for x in rs if not x.is_secondary and x.query_sequence), None)
            if p0 is not None:
                seq = p0.query_sequence if p0.is_reverse == a.is_reverse else p0.query_sequence.upper().translate(RC)[::-1]
        for k in hits[i]:
            L = loci[k]
            # ---- old rule (verbatim logic of trna_count.py r2)
            if L["st"] == "+" and a.reference_end == L["e"]:
                reach = True
                if ct[-1][0] == 4:
                    clip = True; cl = seq[len(seq) - ct[-1][1]:]
                    if cl.upper().startswith("CCA"): cca = True
            elif L["st"] == "-" and a.reference_start == L["s"]:
                reach = True
                if ct[0][0] == 4:
                    clip = True; cl = seq[:ct[0][1]].upper().translate(RC)[::-1]
                    if cl.startswith("CCA"): cca = True
            # ---- new rule
            if not seq: continue
            r = new_tail(a, seq.upper(), L)
            if r is None: continue
            tail, tmpl, past = r; reach_new = True; past_any = past_any or past; c = tclass(tail)
            if c == "CCA" and not tmpl: tmpl_all = False
            if best is None or RANK[c] > RANK[best]: best = c
    fam[F]["reach_3prime"] += reach; fam[F]["clipped_at_3prime"] += clip; fam[F]["CCA_reads"] += cca
    if reach_new:
        fam[F]["reach_3prime_new"] += 1; fam[F]["tail_" + best] += 1
        if best == "CCA": fam[F]["tail_CCA_templated_ambiguous" if tmpl_all else "tail_CCA_non_templated"] += 1
    if not past_any and (best == "CCA") != cca: acc["INVARIANT_MISMATCH"] += 1; fam[F]["invariant_mismatch"] += 1
with pysam.AlignmentFile(BAM, "rb", check_sq=False) as f:
    for r in f.fetch(until_eof=True):
        acc["records"] += 1
        if r.is_unmapped or r.is_supplementary: continue
        nh = r.get_tag("NH"); q = r.query_name; pending.setdefault(q, []).append(r)
        if len(pending[q]) == nh: flush(pending.pop(q))
    for rs in pending.values(): acc["reads_records_lt_NH"] += 1; flush(rs)
json.dump(dict(library=LIB, bam=BAM, kind=KIND, accounting=dict(acc), trna_families={F: dict(v) for F, v in sorted(fam.items())}),
          open(f"{OUT}/{LIB}.cca2.json", "w"), sort_keys=True, indent=0)
print(json.dumps(acc, sort_keys=True))
