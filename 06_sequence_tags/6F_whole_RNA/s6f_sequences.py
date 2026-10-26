#!/usr/bin/env python3
"""s6f_sequences.py -- Step 6F (v2, 2026-10-13), stage F1. One sequence set per gene model, transcript strand, 5'->3'.
Plan 00_PLAN_DRAFT_whole_RNA_tag_2026-10-13.md (Q3 b: UTRs = final.gtf 5_prime_UTR / 3_prime_UTR rows only; the
5_prime_extension / 3_prime_extension rows are NOT used).
Usage: s6f_sequences.py <bin100 master annotation (unmerged)> <FL r3 new_bins100.tsv> <final.gtf> <genome.fasta> <outdir>
Regions per model:
- PCG_CDS: CDS = the model's frozen 100-nt-bin genomic blocks (= CDS frame, FL r3 CDS for the 26 added models), in bin order;
  UTR5 / UTR3 = final.gtf rows of that feature with transcript_id == model_id (none -> region absent);
  whole = UTR5 + CDS + UTR3 (only the regions that exist).
- lncRNA: whole = the model's bin blocks (whole transcript); no UTR / CDS split.
Coordinates: bins are 0-based half-open; GTF is 1-based closed (start-1 converts).
Output: <outdir>/gene_sequences.tsv.gz (one row per model) + checks.tsv (gates + info)."""
import sys, csv, gzip, re, collections
A100, N100, GTF, FA, OUT = sys.argv[1:6]
RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")
rc = lambda s: s.translate(RC)[::-1]
G = {}; name = None; buf = []
for l in open(FA):
    if l.startswith(">"):
        if name: G[name] = "".join(buf).upper()
        name = l[1:].split()[0]; buf = []
    else: buf.append(l.strip())
if name: G[name] = "".join(buf).upper()
M = collections.OrderedDict()
for p in (A100, N100):
    rows = list(csv.DictReader(open(p), delimiter="\t"))
    rows.sort(key=lambda r: (r["model_id"], int(r["bin_number"])))
    for r in rows:
        m = M.setdefault(r["model_id"], dict(model_id=r["model_id"], gene_id=r["gene_id"], rna_class=r["rna_class"], role=r["role"],
                         chrom=r["chrom"], strand=r["strand"], transcript_length=int(r["transcript_length"]), blocks=[],
                         same_strand_overlap=r.get("same_strand_binned_model_overlap", ""), source=p.split("/")[-1]))
        for blk in r["genomic_blocks_5to3"].split(";"):
            s, e = map(int, blk.split("-")); b = m["blocks"]
            if b and ((m["strand"] == "+" and b[-1][1] == s) or (m["strand"] == "-" and b[-1][0] == e)):   # join adjacent bins
                b[-1] = (min(b[-1][0], s), max(b[-1][1], e))
            else: b.append((s, e))
U = collections.defaultdict(list)   # (tid, feature) -> blocks
EXT = collections.Counter()
for l in open(GTF):
    if l.startswith("#"): continue
    x = l.rstrip("\n").split("\t")
    if len(x) < 9: continue
    t = re.search(r'transcript_id "([^"]+)"', x[8])
    if not t: continue
    if x[2] in ("5_prime_UTR", "3_prime_UTR"): U[(t.group(1), x[2])].append((x[0], int(x[3]) - 1, int(x[4]), x[6], x[1]))
    elif x[2] in ("5_prime_extension", "3_prime_extension"): EXT[x[2]] += 1
def seq(ch, blocks, st):
    s = "".join(G[ch][a:b] for a, b in blocks)
    return rc(s) if st == "-" else s
gc = lambda s: (s.count("G") + s.count("C")) / len(s) if s else ""
nonacgt = lambda s: len(s) - sum(s.count(c) for c in "ACGT")
checks = []; ck = lambda n, ok, d="": checks.append((n, "PASS" if ok else "FAIL", str(d)))
info = collections.Counter(); rows = []; bad_len = []; utr_bad = []
for m in M.values():
    st, ch = m["strand"], m["chrom"]
    blocks = m["blocks"]   # 5'->3' order
    core = seq(ch, sorted(blocks), st)
    if len(core) != m["transcript_length"]: bad_len.append(m["model_id"])
    r = dict(model_id=m["model_id"], gene_id=m["gene_id"], rna_class=m["rna_class"], role=m["role"], chrom=ch, strand=st,
             same_strand_binned_model_overlap=m["same_strand_overlap"], annotation_source=m["source"])
    if m["rna_class"] == "PCG_CDS":
        reg = {"CDS": core}
        cs, ce = min(a for a, b in blocks), max(b for a, b in blocks)
        for feat, key in (("5_prime_UTR", "UTR5"), ("3_prime_UTR", "UTR3")):
            u = U.get((m["model_id"], feat), [])
            if not u: continue
            if any(z[0] != ch or z[3] != st for z in u): utr_bad.append((m["model_id"], feat, "chrom/strand")); continue
            ub = sorted((a, b) for _, a, b, _, _ in u)
            if any(min(b, ce) > max(a, cs) for a, b in ub): utr_bad.append((m["model_id"], feat, "overlaps_CDS")); continue
            up = (feat == "5_prime_UTR") == (st == "+")   # True -> UTR lies at lower coordinates than the CDS
            if (up and max(b for a, b in ub) > cs) or (not up and min(a for a, b in ub) < ce):
                utr_bad.append((m["model_id"], feat, "wrong_side")); continue
            adj = (max(b for a, b in ub) == cs) if up else (min(a for a, b in ub) == ce)
            info[f"{key}_adjacent_to_CDS"] += adj; info[f"{key}_gap_to_CDS"] += (not adj)
            info[f"{key}_rows_{len(u)}"] += 1
            info[f"{key}_source_" + ",".join(sorted({z[4] for z in u}))] += 1
            reg[key] = seq(ch, ub, st)
        reg["whole"] = reg.get("UTR5", "") + reg["CDS"] + reg.get("UTR3", "")
        info["CDS_starts_ATG"] += core[:3] == "ATG"; info["CDS_ends_stop"] += core[-3:] in ("TAA", "TAG", "TGA")
        info["CDS_len_multiple_of_3"] += len(core) % 3 == 0
    else:
        reg = {"whole": core}
    for k in ("UTR5", "CDS", "UTR3", "whole"):
        s = reg.get(k, "")
        r[f"len_{k}"] = len(s); r[f"gc_{k}"] = gc(s); r[f"nonACGT_{k}"] = nonacgt(s) if s else ""; r[f"seq_{k}"] = s
    rows.append(r)
P = [r for r in rows if r["rna_class"] == "PCG_CDS"]; L = [r for r in rows if r["rna_class"] == "lncRNA"]
ck("F1_all_models_read", len(rows) == len(M), f"{len(P)} PCG_CDS, {len(L)} lncRNA")
ck("F1_length_equals_transcript_length", not bad_len, ";".join(bad_len[:10]))
ck("F1_utr_rows_consistent (same chrom/strand, outside CDS, correct side)", not utr_bad, ";".join(map(str, utr_bad[:10])))
ck("F1_chromosomes_in_genome", all(r["chrom"] in G for r in rows))
for k in ("UTR5", "UTR3"):
    ck(f"INFO_PCG_with_{k}", True, sum(1 for r in P if r[f"len_{k}"] > 0))
ck("INFO_PCG_with_both_UTRs", True, sum(1 for r in P if r["len_UTR5"] > 0 and r["len_UTR3"] > 0))
ck("INFO_gtf_extension_rows_not_used", True, dict(EXT))
for k, v in sorted(info.items()): ck("INFO_" + k, True, v)
for k in ("UTR5", "CDS", "UTR3", "whole"):
    xs = sorted(r[f"len_{k}"] for r in P if r[f"len_{k}"] > 0)
    if xs: ck(f"INFO_PCG_{k}_length_min_median_max", True, f"{xs[0]}/{xs[len(xs)//2]}/{xs[-1]}")
    ck(f"INFO_PCG_{k}_nonACGT_bases", True, sum(r[f"nonACGT_{k}"] or 0 for r in P))
ck("INFO_lncRNA_whole_nonACGT_bases", True, sum(r["nonACGT_whole"] or 0 for r in L))
cols = list(rows[0].keys())
with gzip.open(f"{OUT}/gene_sequences.tsv.gz", "wt") as w:
    dw = csv.DictWriter(w, cols, delimiter="\t", lineterminator="\n"); dw.writeheader(); dw.writerows(rows)
with open(f"{OUT}/checks.tsv", "w") as w:
    w.write("check\tresult\tdetail\n"); [w.write("\t".join(c) + "\n") for c in checks]
for c in checks: print(*c, sep="\t")
ok = all(c[1] == "PASS" for c in checks); print("F1_ALL_GATES_PASS" if ok else "F1_GATE_FAILURE"); sys.exit(0 if ok else 1)
