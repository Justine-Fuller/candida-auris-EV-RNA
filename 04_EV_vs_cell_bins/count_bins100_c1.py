#!/usr/bin/env python3
"""count_bins100_c1.py -- Option-1 100-nt branch INDEPENDENT counter C1 (2026-10-08).
r4 (2026-10-09): optional 6th argument = overhang rule: 'half' (default; the registered rule) or 'any' (SENSITIVITY
check only, Justine: exclude a read with ANY aligned base beyond the 3' end). Default behaviour unchanged.
r3 (2026-10-09, Caro rule; definitions Justine): a read assigned to one transcript is EXCLUDED when HALF OR MORE of
its aligned bases (M/=/X) lie beyond that transcript's 3' end (PCG: past the CDS stop end, i.e. 3' UTR side;
lncRNA: past the model's 3' end). 5' overhang and intronic bases are ignored for this rule; soft clips never count.
Counted as excluded_3prime_overhang (before tie evaluation), listed in <lib>.excluded_3prime.tsv.gz, and per
transcript in transcript_totals column excluded_3prime_overhang_reads.
 r2: also accepts Cell sRNA (kind "sRNA"). Derived from frozen
count_bins_c1.py (150-nt branch, md5 4e5dc9cf...): per-(chrom, strand) numpy position->bin maps (positions shared by
> 1 bin kept in a side table) and pysam get_blocks() (gapless aligned blocks; D and N excluded) instead of CIGAR parsing
and a bucketed interval index. Same eligibility/strand/model rules and the same Option-1 tie handling as C0
(count_bins100.py): ties at the greatest overlap are never counted in bins or transcript totals; they are listed.
Writes the same files as C0 with suffix .c1; the R2 gate requires them to be byte-identical."""
import sys, json, gzip, collections, pysam, numpy as np
bins_tsv, bam, lib, kind, outp = sys.argv[1:6]
RULE = sys.argv[6] if len(sys.argv) > 6 else "half"   # r4
assert RULE in ("half", "any"), "overhang rule must be half or any"
assert kind in ("EV", "sRNA"), "small-RNA libraries only (EV, Cell sRNA: both sense-stranded); polyA refused"   # r2: Cell sRNA added
NEX = 200
bid, bmodel, brole, blocks = [], [], [], []
with open(bins_tsv) as f:
    h = f.readline().rstrip("\n").split("\t"); c = {k: i for i, k in enumerate(h)}
    for line in f:
        r = line.rstrip("\n").split("\t")
        bid.append(r[c["bin_id"]]); bmodel.append(r[c["model_id"]]); brole.append(r[c["role"]])
        blocks.append((r[c["chrom"]], r[c["strand"]], [tuple(map(int, b.split("-"))) for b in r[c["genomic_blocks_5to3"]].split(";")]))
with pysam.AlignmentFile(bam, "rb", check_sq=False) as f:
    clen = dict(zip(f.references, f.lengths))
maps, multi = {}, collections.defaultdict(list)
for k, (ch, st, bl) in enumerate(blocks):
    a = maps.setdefault((ch, st), np.zeros(clen[ch], dtype=np.int64))
    for s, e in bl:
        for p in range(s, e):
            v = a[p]
            if v == 0: a[p] = k + 1
            else:
                if v > 0: multi[(ch, st, p)].append(v - 1); a[p] = -1
                multi[(ch, st, p)].append(k)
acc = collections.Counter(); bc = collections.Counter(); phys = collections.Counter(); tie = collections.Counter()
mpos = collections.defaultdict(list)   # r3: all exonic genomic block coordinates per model -> 3' boundary
for k, (ch, st, bl) in enumerate(blocks): mpos[bmodel[k]].extend(x for se in bl for x in se)
mstrand = {bmodel[k]: blocks[k][1] for k in range(len(bmodel))}
bound3 = {m: (max(v) if mstrand[m] == "+" else min(v)) for m, v in mpos.items()}
x3 = collections.Counter(); x3rows = []
pat = collections.Counter(); ties = []; ex = []
with pysam.AlignmentFile(bam, "rb", check_sq=False) as f:
    for r in f.fetch(until_eof=True):
        acc["records"] += 1
        if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
        acc["primary_mapped"] += 1
        if r.get_tag("NH") != 1: continue
        rl = r.infer_read_length() or r.query_length
        if rl < 18: continue
        acc["eligible_NH1_ge18"] += 1
        st = "-" if r.is_reverse else "+"
        a = maps.get((r.reference_name, st)); ov = collections.Counter()
        if a is not None:
            for s, e in r.get_blocks():
                seg = a[s:e]
                for v, n in zip(*np.unique(seg[seg > 0], return_counts=True)): ov[int(v) - 1] += int(n)
                for p in (np.nonzero(seg < 0)[0] + s):
                    for k in multi[(r.reference_name, st, int(p))]: ov[k] += 1
        if not ov: acc["no_binned_model"] += 1; continue
        if len({bmodel[k] for k in ov}) > 1: acc["ambiguous_model"] += 1; continue
        if any(brole[k] == "competitor_only" for k in ov): acc["competitor_only"] += 1; continue
        mm = bmodel[next(iter(ov))]; gb = r.get_blocks(); al = sum(e - s for s, e in gb); bd = bound3[mm]
        bey = sum(len(range(max(s, bd), e)) for s, e in gb) if mstrand[mm] == "+" else sum(len(range(s, min(e, bd))) for s, e in gb)
        if (2 * bey >= al) if RULE == "half" else (bey > 0):
            acc["excluded_3prime_overhang"] += 1; x3[mm] += 1
            x3rows.append((lib, r.query_name, r.reference_name, r.reference_start, st, r.cigarstring, rl, mm, al, bey)); continue
        vals = sorted(((-v, bid[k], k) for k, v in ov.items()))
        mx = -vals[0][0]; top = [k for nv, _, k in vals if -nv == mx]; m = bmodel[vals[0][2]]
        desc = ";".join(f"{b}:{-nv}" for nv, b, _ in vals)
        if len(vals) >= 2: pat[(len(vals), tuple(-nv for nv, _, _ in vals), len(top) >= 2)] += 1
        row = (lib, r.query_name, r.reference_name, r.reference_start, st, r.cigarstring, rl, m, len(vals), len(top), mx, desc)
        if len(top) >= 2:
            acc["tie_unresolved_2bin" if len(top) == 2 else "tie_unresolved_3plus"] += 1; tie[m] += 1; ties.append(row); continue
        acc["assigned"] += 1; phys[m] += 1; bc[bid[top[0]]] += 1
        if len(vals) >= 2:
            acc["assigned_touching_2plus_bins"] += 1
            if len(ex) < NEX: ex.append(row + (bid[top[0]],))
H = "library\tread\tchrom\tstart0\trna_strand\tcigar\tread_length\tmodel_id\tn_bins_touched\tn_bins_at_max\tmax_overlap\tbin_overlaps_desc"
# transcript summed_bin_counts recomputed from the bin counter itself (independent of phys)
sb = collections.Counter()
for b, n in bc.items(): sb[b.split("|")[0]] += n
w = lambda name: gzip.open(f"{outp}/{lib}.{name}.c1.tsv.gz", "wt")
with w("bin_counts") as o:
    o.write("bin_id\tassigned_reads\n"); [o.write(f"{b}\t{bc[b]}\n") for b in sorted(bc)]
with w("transcript_totals") as o:
    o.write("model_id\tassigned_reads\tsummed_bin_counts\tunresolved_tie_reads\texcluded_3prime_overhang_reads\n")
    [o.write(f"{m}\t{phys[m]}\t{sb[m]}\t{tie[m]}\t{x3[m]}\n") for m in sorted(set(phys) | set(tie) | set(x3))]
with w("excluded_3prime") as o:
    o.write("library\tread\tchrom\tstart0\trna_strand\tcigar\tread_length\tmodel_id\taligned_bases\tbases_beyond_3prime_end\n")
    [o.write("\t".join(map(str, t)) + "\n") for t in x3rows]
with w("ties") as o:
    o.write(H + "\n"); [o.write("\t".join(map(str, t)) + "\n") for t in ties]
with w("multibin_examples") as o:
    o.write(H + "\tassigned_bin\n"); [o.write("\t".join(map(str, t)) + "\n") for t in ex]
with w("overlap_patterns") as o:
    o.write("n_bins_touched\toverlaps_desc\tis_tie\treads\n")
    [o.write(f"{n}\t{','.join(map(str, v))}\t{t}\t{cnt}\n") for (n, v, t), cnt in sorted(pat.items())]
acc["overhang_rule"] = RULE
json.dump(dict(acc), open(f"{outp}/{lib}.accounting.c1.json", "w"), indent=1, sort_keys=True); print(json.dumps(dict(acc), sort_keys=True))
