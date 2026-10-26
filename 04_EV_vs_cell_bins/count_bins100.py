#!/usr/bin/env python3
"""count_bins100.py -- Option-1 100-nt branch counter C0 (2026-10-08).
r4 (2026-10-09): optional 6th argument = overhang rule: 'half' (default; the registered rule) or 'any' (SENSITIVITY
check only, Justine: exclude a read with ANY aligned base beyond the 3' end). Default behaviour unchanged.
r3 (2026-10-09, Caro rule; definitions Justine): a read assigned to one transcript is EXCLUDED when HALF OR MORE of
its aligned bases (M/=/X) lie beyond that transcript's 3' end (PCG: past the CDS stop end, i.e. 3' UTR side;
lncRNA: past the model's 3' end). 5' overhang and intronic bases are ignored for this rule; soft clips never count.
Counted as excluded_3prime_overhang (before tie evaluation), listed in <lib>.excluded_3prime.tsv.gz, and per
transcript in transcript_totals column excluded_3prime_overhang_reads.
 r2 (2026-10-08): also accepts Cell sRNA libraries
(kind "sRNA"; TruSeq Small RNA, sense-stranded exactly like EV, frozen rule); nothing else changed. Derived from frozen count_bins.py r2 (150-nt
branch, md5 b1fcb232...). Read eligibility, strand and model rules are UNCHANGED:
  primary, mapped, NH == 1, trimmed read length (incl. soft clips) >= 18; RNA strand = read strand (EV small RNA, sense);
  aligned overlap = CIGAR M/=/X reference bases on bin exonic positions (D, N, I, S contribute 0); read touching > 1
  same-strand binned model = ambiguous_model (unassigned); touching only a competitor-only model = competitor_only
  (unassigned); 10-kb index buckets de-duplicated per read block (r2 fix).
CHANGED for Option 1 (no tie rule decided yet): a read is assigned (+1) only to the single bin with the strictly
greatest aligned overlap. A read whose greatest overlap is shared EXACTLY by >= 2 bins is NOT counted in any bin and
NOT in the transcript's assigned total; it is written to <lib>.ties.tsv.gz (all details) and counted per transcript as
unresolved_tie_reads, so a tie rule can be applied later without recounting.
Outputs (outdir): <lib>.bin_counts.tsv.gz (bin_id, assigned_reads; non-zero only);
  <lib>.transcript_totals.tsv.gz (model_id, assigned_reads, summed_bin_counts, unresolved_tie_reads);
  <lib>.ties.tsv.gz (one row per tied read); <lib>.multibin_examples.tsv.gz (first 200 resolved reads touching >= 2 bins
  of one model, BAM order); <lib>.overlap_patterns.tsv.gz (histogram of per-read overlap vectors for reads touching
  >= 2 bins of one model); <lib>.accounting.json."""
import sys, json, gzip, collections, pysam
bins_tsv, bam, lib, kind, outp = sys.argv[1:6]
RULE = sys.argv[6] if len(sys.argv) > 6 else "half"   # r4
assert RULE in ("half", "any"), "overhang rule must be half or any"
assert kind in ("EV", "sRNA"), "small-RNA libraries only (EV, Cell sRNA: both sense-stranded); polyA refused"   # r2: Cell sRNA added
B = 10000; NEX = 200
idx = collections.defaultdict(list); binid = []; binmodel = []; binrole = []
end3 = {}   # r3: model -> (strand, 3' boundary): "+" positions >= boundary, "-" positions < boundary are beyond the 3' end
with open(bins_tsv) as f:
    h = f.readline().rstrip("\n").split("\t"); c = {k: i for i, k in enumerate(h)}
    for line in f:
        r = line.rstrip("\n").split("\t"); k = len(binid)
        binid.append(r[c["bin_id"]]); binmodel.append(r[c["model_id"]]); binrole.append(r[c["role"]])
        for blk in r[c["genomic_blocks_5to3"]].split(";"):
            s, e = map(int, blk.split("-"))
            for b in range(s // B, (e - 1) // B + 1):
                idx[(r[c["chrom"]], r[c["strand"]], b)].append((s, e, k))
            m_, st_ = r[c["model_id"]], r[c["strand"]]; old = end3.get(m_, (st_, e if st_ == "+" else s))[1]
            end3[m_] = (st_, max(old, e) if st_ == "+" else min(old, s))
acc = collections.Counter(); bc = collections.Counter(); phys = collections.Counter(); summ = collections.Counter()
tie = collections.Counter(); pat = collections.Counter(); ties = []; ex = []; x3 = collections.Counter(); x3rows = []
fmt = lambda ks, ov: ";".join(f"{binid[k]}:{ov[k]}" for k in ks)
with pysam.AlignmentFile(bam, "rb", check_sq=False) as f:
    for r in f.fetch(until_eof=True):
        acc["records"] += 1
        if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
        acc["primary_mapped"] += 1
        if r.get_tag("NH") != 1: continue
        rl = r.infer_read_length() or r.query_length
        if rl < 18: continue
        acc["eligible_NH1_ge18"] += 1
        rs = "-" if r.is_reverse else "+"
        ov = collections.Counter(); pos = r.reference_start; mblk = []
        for op, ln in r.cigartuples:
            if op in (0, 7, 8):
                s, e = pos, pos + ln; seen = set(); mblk.append((s, e))
                for b in range(s // B, (e - 1) // B + 1):
                    for (bs, be, k) in idx.get((r.reference_name, rs, b), ()):
                        if (bs, be, k) in seen: continue
                        seen.add((bs, be, k))
                        o = min(e, be) - max(s, bs)
                        if o > 0: ov[k] += o
                pos += ln
            elif op in (2, 3): pos += ln
        if not ov: acc["no_binned_model"] += 1; continue
        models = {binmodel[k] for k in ov}
        if len(models) > 1: acc["ambiguous_model"] += 1; continue
        if binrole[next(iter(ov))] == "competitor_only": acc["competitor_only"] += 1; continue
        mm = binmodel[next(iter(ov))]; st3, bd = end3[mm]; al = sum(e - s for s, e in mblk)   # r3: 3' overhang rule
        bey = sum(max(0, e - max(s, bd)) for s, e in mblk) if st3 == "+" else sum(max(0, min(e, bd) - s) for s, e in mblk)
        if (2 * bey >= al) if RULE == "half" else (bey > 0):
            acc["excluded_3prime_overhang"] += 1; x3[mm] += 1
            x3rows.append((lib, r.query_name, r.reference_name, r.reference_start, rs, r.cigarstring, rl, mm, al, bey)); continue
        ks = sorted(ov, key=lambda k: (-ov[k], binid[k]))            # bins by overlap desc, then bin id (5'->3')
        mx = ov[ks[0]]; top = [k for k in ks if ov[k] == mx]; m = binmodel[ks[0]]
        if len(ks) >= 2: pat[(len(ks), tuple(ov[k] for k in ks), len(top) >= 2)] += 1
        row = (lib, r.query_name, r.reference_name, r.reference_start, rs, r.cigarstring, rl, m, len(ks), len(top), mx, fmt(ks, ov))
        if len(top) >= 2:
            acc["tie_unresolved_2bin" if len(top) == 2 else "tie_unresolved_3plus"] += 1; tie[m] += 1; ties.append(row); continue
        acc["assigned"] += 1; phys[m] += 1; summ[m] += 1; bc[binid[top[0]]] += 1
        if len(ks) >= 2:
            acc["assigned_touching_2plus_bins"] += 1
            if len(ex) < NEX: ex.append(row + (binid[top[0]],))
H = "library\tread\tchrom\tstart0\trna_strand\tcigar\tread_length\tmodel_id\tn_bins_touched\tn_bins_at_max\tmax_overlap\tbin_overlaps_desc"
with gzip.open(f"{outp}/{lib}.bin_counts.tsv.gz", "wt") as o:
    o.write("bin_id\tassigned_reads\n")
    for b in sorted(bc): o.write(f"{b}\t{bc[b]}\n")
with gzip.open(f"{outp}/{lib}.transcript_totals.tsv.gz", "wt") as o:
    o.write("model_id\tassigned_reads\tsummed_bin_counts\tunresolved_tie_reads\texcluded_3prime_overhang_reads\n")
    for m in sorted(set(phys) | set(tie) | set(x3)): o.write(f"{m}\t{phys[m]}\t{summ[m]}\t{tie[m]}\t{x3[m]}\n")
with gzip.open(f"{outp}/{lib}.excluded_3prime.tsv.gz", "wt") as o:
    o.write("library\tread\tchrom\tstart0\trna_strand\tcigar\tread_length\tmodel_id\taligned_bases\tbases_beyond_3prime_end\n")
    for t in x3rows: o.write("\t".join(map(str, t)) + "\n")
with gzip.open(f"{outp}/{lib}.ties.tsv.gz", "wt") as o:
    o.write(H + "\n")
    for t in ties: o.write("\t".join(map(str, t)) + "\n")
with gzip.open(f"{outp}/{lib}.multibin_examples.tsv.gz", "wt") as o:
    o.write(H + "\tassigned_bin\n")
    for t in ex: o.write("\t".join(map(str, t)) + "\n")
with gzip.open(f"{outp}/{lib}.overlap_patterns.tsv.gz", "wt") as o:
    o.write("n_bins_touched\toverlaps_desc\tis_tie\treads\n")
    for (n, v, t), cnt in sorted(pat.items()): o.write(f"{n}\t{','.join(map(str, v))}\t{t}\t{cnt}\n")
acc = dict(acc); acc.update(library=lib, kind=kind, bam=bam, overhang_rule=RULE, bins_with_counts=len(bc), transcripts_with_assigned_reads=len(phys),
                            transcripts_with_unresolved_ties=len(tie),
                            closes=acc.get("eligible_NH1_ge18", 0) == sum(acc.get(k, 0) for k in (
                                "no_binned_model", "ambiguous_model", "competitor_only", "excluded_3prime_overhang", "tie_unresolved_2bin", "tie_unresolved_3plus", "assigned")))
json.dump(acc, open(f"{outp}/{lib}.accounting.json", "w"), indent=1, sort_keys=True); print(json.dumps(acc, sort_keys=True))
