#!/usr/bin/env python3
"""
build_bin100.py -- Option-1 100-nt branch, 2026-10-08. Derived from frozen build_bin150.py (md5 846487725ad7...) by changing ONLY
  the bin size (150 -> 100) and output names; eligibility screen and model construction unchanged. Annotation only.

Supervisor-requested 150-nt bin analysis, part 1 + 2:
  (1) lncRNA eligibility screen under Caro's clarified rules (catalogue strand authoritative; exclude only
      SAME-strand overlap with PCG CDS/UTR, high-confidence tRNA (comprehensive genome-wide set), rRNA, OTHER_NC;
      opposite-strand overlap is not an exclusion; no catalogue_antisense units);
  (2) fixed 100-nt bins in spliced transcript coordinates for PCG (frozen primary CDS model) + eligible lncRNAs.

Does NOT modify any frozen file; reads read-only copies in ../inputs/ (md5-checked).
Coordinates: everything internal and in outputs is 0-based half-open [start, end).
Overlap: >= 1 shared nucleotide (identical to the frozen eligibility rule).
"""
import collections
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
IN = os.path.join(HERE, "..", "inputs")
OUT = os.path.join(HERE, "..", "outputs")
BIN = 100

EXPECTED_MD5 = {
    "intervals.tsv": "d83a0d1ddbc69524e729a365b9ab3906",
    "candidate_trna_reference_168.bed": "b2051a075b9bc1b8225c47d3d29b916e",
    "lncRNA_rDNA_provenance.tsv": "35a21820a7433bee9f50d5b3880c143f",
    "StringTie_rDNA_provenance.tsv": "cbf9d342961635ac73102734c04b2735",
    "flag_sets.json": "dfcce2d2d24d18762cac7fafb31f9d20",
    "v2_lncRNA_FROZEN_units.tsv": "5f7f596935dabc26a671f025a1837203",
}
REPEAT_REG = {"LNCRNA.1946", "LNCRNA.1057"}
PARTIAL_REG = {"LNCRNA.4454", "LNCRNA.4455", "LNCRNA.4897", "LNCRNA.4899", "LNCRNA.4901", "LNCRNA.4902", "LNCRNA.5575"}
GANNOT_PCG_GENES = {"B9J08_002935", "B9J08_005532", "B9J08_005533"}   # frozen gene-level flag (j03_hotspots.py)
GANNOT_LNC = {"LNCRNA.1946": "missing 5S B9J08_001056 provisionally inside (reference note C2); rDNA handling prevails"}
# exclusion classes, in the priority used for the single "primary_reason" column
SCREEN = ["PCG_CDS", "PCG_UTR", "tRNA_HC", "rRNA", "OTHER_NC"]
DESCRIPTIVE = ["tRNA_questionable", "tRNA_oldGTFspan", "STRINGTIE"]
FLIP = {"+": "-", "-": "+"}

gates = []


def gate(name, ok, detail=""):
    gates.append((name, "PASS" if ok else "FAIL", detail))


def md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()


# ------------------------------------------------------------------ inputs
for f, m in EXPECTED_MD5.items():
    got = md5(os.path.join(IN, f))
    gate(f"G0_md5_{f}", got == m, got)

feat = collections.defaultdict(list)   # class -> list of dict(chrom,s,e,strand,id,gene)
lnc = {}
pcg_tx = collections.defaultdict(list)
pcg_gene = {}
with open(os.path.join(IN, "intervals.tsv")) as fh:
    next(fh)
    for line in fh:
        c, s, e, cls, g, tx, st, _ = line.rstrip("\n").split("\t")
        s0, e0 = int(s) - 1, int(e)
        rec = dict(chrom=c, s=s0, e=e0, strand=st, id=g, gene=g)
        if cls == "LNC":
            assert g not in lnc, f"multi-interval LNC {g}"
            lnc[g] = rec
        elif cls == "PCG_CDS":
            pcg_tx[tx].append((c, s0, e0, st)); pcg_gene[tx] = g
            feat["PCG_CDS"].append(dict(rec, id=tx))
        elif cls == "PCG_UTR":
            feat["PCG_UTR"].append(dict(rec, id=tx))
        elif cls == "tRNA":
            feat["tRNA_oldGTFspan"].append(rec)
        elif cls in ("rRNA", "OTHER_NC", "STRINGTIE"):
            feat[cls].append(dict(rec, id=g))
with open(os.path.join(IN, "candidate_trna_reference_168.bed")) as fh:
    next(fh)
    for line in fh:
        c, s, e, name, _sc, st, fam, fid, intron, tier = line.rstrip("\n").split("\t")
        cls = "tRNA_HC" if tier == "A_core" else "tRNA_questionable"
        feat[cls].append(dict(chrom=c, s=int(s), e=int(e), strand=st, id=f"{name.split('|')[0]}|{fam}|{fid}", gene=fid))
gate("G1_counts_input", len(lnc) == 4588 and len(pcg_tx) == 5304 and len(feat["tRNA_HC"]) == 167
     and len(feat["tRNA_questionable"]) == 1,
     f"LNC={len(lnc)} PCGtx={len(pcg_tx)} tRNA_HC={len(feat['tRNA_HC'])} questionable={len(feat['tRNA_questionable'])}")

st_flag = {}
with open(os.path.join(IN, "StringTie_rDNA_provenance.tsv")) as fh:
    h = next(fh).rstrip("\n").split("\t")
    for line in fh:
        r = dict(zip(h, line.rstrip("\n").split("\t")))
        st_flag[r["StringTie_gene"]] = r["flag"]
with open(os.path.join(IN, "lncRNA_rDNA_provenance.tsv")) as fh:
    h = next(fh).rstrip("\n").split("\t")
    prov = [dict(zip(h, l.rstrip("\n").split("\t"))) for l in fh]
rep = {r["lncRNA_model"] for r in prov if r["flag"] == "rDNA_repeat_associated"}
par = {r["lncRNA_model"] for r in prov if r["flag"] == "partial_rDNA_homology_associated"}
gate("G2_rDNA_flags_eq_register", rep == REPEAT_REG and par == PARTIAL_REG)
fs = json.load(open(os.path.join(IN, "flag_sets.json")))
F1, F2 = set(fs["F1"]), set(fs["F2"])


# ------------------------------------------------------------------ overlap index
class Idx:
    def __init__(self, recs):
        self.by = collections.defaultdict(list)
        for r in recs:
            self.by[r["chrom"]].append(r)

    def hits(self, chrom, s, e):
        return [r for r in self.by.get(chrom, ()) if r["s"] < e and s < r["e"]]


idx = {k: Idx(v) for k, v in feat.items()}


def ovl_bp(a_s, a_e, b_s, b_e):
    return max(0, min(a_e, b_e) - max(a_s, b_s))


# ------------------------------------------------------------------ (1) lncRNA screen
audit = []
eligible = []
for g in sorted(lnc, key=lambda x: int(x.split(".")[1])):
    m = lnc[g]
    row = dict(model=g, chrom=m["chrom"], start0=m["s"], end=m["e"], length=m["e"] - m["s"], catalogue_strand=m["strand"])
    same_classes = []
    for cls in SCREEN + DESCRIPTIVE:
        hs = idx[cls].hits(m["chrom"], m["s"], m["e"]) if cls in idx else []
        for side, want in (("same", m["strand"]), ("opp", FLIP[m["strand"]])):
            sel = [h for h in hs if h["strand"] == want]
            ids = sorted({h["id"] for h in sel})
            bp = sum(ovl_bp(m["s"], m["e"], h["s"], h["e"]) for h in sel)
            row[f"{cls}_{side}_n"] = len(ids)
            row[f"{cls}_{side}_ids"] = ",".join(ids)
            row[f"{cls}_{side}_bp_sum"] = bp
            if side == "same" and cls in SCREEN and ids:
                same_classes.append(cls)
    opp_any = [c for c in SCREEN if row[f"{c}_opp_n"]]
    row["same_strand_exclusion_classes"] = ",".join(same_classes)
    row["opposite_strand_overlap_classes"] = ",".join(opp_any)
    row["eligible"] = "yes" if not same_classes else "no"
    row["primary_exclusion_reason"] = (f"same_strand_{same_classes[0]}" if same_classes else "")
    row["rdna_flag"] = ("rDNA_repeat_associated" if g in REPEAT_REG else
                        "partial_rDNA_homology_associated" if g in PARTIAL_REG else "")
    row["role_if_eligible"] = ("competitor_only" if g in REPEAT_REG else "ordinary") if not same_classes else ""
    row["GANNOT_note"] = GANNOT_LNC.get(g, "")
    st_ids = (row["STRINGTIE_same_ids"] + "," + row["STRINGTIE_opp_ids"]).strip(",")
    row["StringTie_rDNA_flags_overlapping"] = ",".join(sorted({f"{x}:{st_flag[x]}" for x in st_ids.split(",") if x in st_flag}))
    row["D4_mappability_F1"] = "yes" if f"{g}-T" in F1 else ""
    row["D4_mappability_F2"] = "yes" if f"{g}-T" in F2 else ""
    audit.append(row)
    if not same_classes:
        eligible.append(g)
elig = set(eligible)

# independent brute-force re-derivation of eligibility (all features, pairwise)
allscreen = [r for c in SCREEN for r in feat[c]]
brute = {g for g, m in lnc.items()
         if not any(r["chrom"] == m["chrom"] and r["strand"] == m["strand"] and r["s"] < m["e"] and m["s"] < r["e"]
                    for r in allscreen)}
gate("G3_eligible_eq_bruteforce", brute == elig, f"{len(elig)} vs {len(brute)}")

# same-strand lncRNA-lncRNA overlaps among eligible models
ll = Idx([dict(lnc[g], id=g) for g in eligible])
pairs = []
for g in eligible:
    m = lnc[g]
    for h in ll.hits(m["chrom"], m["s"], m["e"]):
        if h["id"] != g and h["strand"] == m["strand"] and g < h["id"]:
            pairs.append((g, h["id"], ovl_bp(m["s"], m["e"], h["s"], h["e"])))
in_pair = collections.defaultdict(list)
for a, b, bp in pairs:
    in_pair[a].append(b); in_pair[b].append(a)
for row in audit:
    row["same_strand_eligible_lncRNA_partners"] = ",".join(sorted(in_pair.get(row["model"], []))) if row["eligible"] == "yes" else ""

# old rule (frozen 674) re-derived from the same intervals, and compared with the frozen units file
OLD_EXCL = ("PCG_CDS", "PCG_UTR", "tRNA_oldGTFspan", "rRNA", "OTHER_NC")
old_reason = {}
for g, m in lnc.items():
    r = set()
    for c in OLD_EXCL:
        if idx[c].hits(m["chrom"], m["s"], m["e"]):
            r.add(c)
    old_reason[g] = r
old_elig = {g for g, r in old_reason.items() if not r}
frozen_models = set()
with open(os.path.join(IN, "v2_lncRNA_FROZEN_units.tsv")) as fh:
    next(fh)
    for line in fh:
        frozen_models.add(line.split("\t")[1])
gate("G4_old_rule_reproduces_frozen_674", old_elig == frozen_models and len(frozen_models) == 674,
     f"rederived={len(old_elig)} frozen={len(frozen_models)}")

for row in audit:
    g = row["model"]
    o, n = g in old_elig, g in elig
    row["old_674_eligible"] = "yes" if o else "no"
    row["old_exclusion_classes_either_strand"] = ",".join(sorted(old_reason[g]))
    if o and n:
        cmp_ = "shared"
    elif n and not o:
        cmp_ = "newly_recovered"
    elif o and not n:
        cmp_ = "newly_excluded"
    else:
        cmp_ = "excluded_both"
    row["old_vs_new"] = cmp_

# ------------------------------------------------------------------ (2) bins
def make_model(chrom, strand, exons):
    gs = sorted(exons)
    tord = gs if strand == "+" else gs[::-1]
    segs, off = [], 0
    for (s, e) in tord:
        segs.append((s, e, off)); off += e - s
    return dict(chrom=chrom, strand=strand, exons=gs, segs=segs, length=off)


def tx_to_genome(model, a, b):
    """transcript [a,b) -> genomic blocks (0-based half-open) listed in transcript 5'->3' order."""
    out = []
    for s, e, o in model["segs"]:
        lo, hi = max(a, o), min(b, o + e - s)
        if lo >= hi:
            continue
        if model["strand"] == "+":
            out.append((s + lo - o, s + hi - o))
        else:
            out.append((e - (hi - o), e - (lo - o)))
    return out


models = []
for tx in sorted(pcg_tx):
    segs = pcg_tx[tx]
    assert len({x[0] for x in segs}) == 1 and len({x[3] for x in segs}) == 1, tx
    m = make_model(segs[0][0], segs[0][3], [(x[1], x[2]) for x in segs])
    gid = pcg_gene[tx]
    m.update(rna_class="PCG_CDS", model_id=tx, gene_id=gid, role="ordinary", rdna_flag="",
             gannot=("gene-level ambiguity: frozen reference note C2 (not direct overlap)" if gid in GANNOT_PCG_GENES else ""),
             F1="yes" if tx in F1 else "", F2="yes" if tx in F2 else "")
    models.append(m)
for g in eligible:
    x = lnc[g]
    m = make_model(x["chrom"], x["strand"], [(x["s"], x["e"])])
    m.update(rna_class="lncRNA", model_id=g, gene_id=g,
             role="competitor_only" if g in REPEAT_REG else "ordinary",
             rdna_flag="rDNA_repeat_associated" if g in REPEAT_REG else ("partial_rDNA_homology_associated" if g in PARTIAL_REG else ""),
             gannot=GANNOT_LNC.get(g, ""), F1="yes" if f"{g}-T" in F1 else "", F2="yes" if f"{g}-T" in F2 else "")
    models.append(m)
gate("G5_PCG_CDS_total_nt", sum(m["length"] for m in models if m["rna_class"] == "PCG_CDS") == 7953370)

# index of all binned models' exons for same-strand cross-model overlap flags
mx = collections.defaultdict(list)
for k, m in enumerate(models):
    for (s, e) in m["exons"]:
        mx[m["chrom"]].append((s, e, m["strand"], k))
st_idx = idx["STRINGTIE"]

bins = []
v_ok = True
for k, m in enumerate(models):
    L = m["length"]
    nb = (L + BIN - 1) // BIN
    covered = []
    for i in range(nb):
        a, b = i * BIN, min((i + 1) * BIN, L)
        blocks = tx_to_genome(m, a, b)
        covered.extend(blocks)
        # bin-level flags
        same_other, trna_same, trna_opp, st_hits = set(), set(), set(), set()
        for (s, e) in blocks:
            for (s2, e2, st2, k2) in mx[m["chrom"]]:
                if k2 != k and st2 == m["strand"] and s2 < e and s < e2:
                    same_other.add(models[k2]["model_id"])
            for h in idx["tRNA_HC"].hits(m["chrom"], s, e):
                (trna_same if h["strand"] == m["strand"] else trna_opp).add(h["id"])
            for h in st_idx.hits(m["chrom"], s, e):
                st_hits.add(h["id"] + ("(" + st_flag[h["id"]] + ")" if h["id"] in st_flag else ""))
        bins.append(dict(
            bin_id=f"{m['model_id']}|bin{i + 1:03d}", rna_class=m["rna_class"], model_id=m["model_id"],
            gene_id=m["gene_id"], chrom=m["chrom"], strand=m["strand"], role=m["role"],
            transcript_length=L, n_bins=nb, bin_number=i + 1, tx_start0=a, tx_end=b, bin_length=b - a,
            final_bin_shorter=("yes" if (i == nb - 1 and b - a < BIN) else ""),
            single_bin_no_positional_contrast=("yes" if L <= BIN else ""),
            genomic_blocks_5to3=";".join(f"{s}-{e}" for s, e in blocks),
            rdna_flag=m["rdna_flag"], GANNOT=m["gannot"],
            D4_mappability_F1_transcript=m["F1"], D4_mappability_F2_transcript=m["F2"],
            poorly_mappable_v2="NOT_DEFINED(register open item)",
            same_strand_binned_model_overlap=",".join(sorted(same_other)),
            HC_tRNA_same_strand=",".join(sorted(trna_same)), HC_tRNA_opposite_strand=",".join(sorted(trna_opp)),
            StringTie_overlap_either_strand=",".join(sorted(st_hits))))
    # V1: bin lengths sum to transcript length; all bins full except possibly the last
    lens = [bb["bin_length"] for bb in bins[-nb:]]
    if sum(lens) != L or any(x != BIN for x in lens[:-1]) or not (0 < lens[-1] <= BIN):
        v_ok = False
    # V2 (independent per-base): enumerate model bases in transcript order directly from exons;
    # base j must fall in bin j//BIN and be listed exactly once across that model's bin blocks
    bases = []
    for (s, e) in m["exons"]:
        bases.extend(range(s, e))
    if m["strand"] == "-":
        bases.reverse()
    perbin = []
    for bb in bins[-nb:]:
        bl = []
        for blk in bb["genomic_blocks_5to3"].split(";"):
            s, e = map(int, blk.split("-"))
            rng = range(s, e) if m["strand"] == "+" else range(e - 1, s - 1, -1)
            bl.extend(rng)
        perbin.append(bl)
    flat = [g for bl in perbin for g in bl]
    if flat != bases or len(set(flat)) != len(flat) or any(len(bl) != lens[i] for i, bl in enumerate(perbin)):
        v_ok = False
gate("G6_bins_partition_every_position_exactly_once_and_sum_to_length", v_ok)
gate("G7_lncRNA_models_binned_eq_eligible",
     {m["model_id"] for m in models if m["rna_class"] == "lncRNA"} == elig)

# ------------------------------------------------------------------ write
os.makedirs(OUT, exist_ok=True)


def write(path, rows, cols=None):
    cols = cols or list(rows[0].keys())
    with open(path, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows:
            o.write("\t".join(str(r[c]) for c in cols) + "\n")


write(os.path.join(OUT, "lncRNA_eligibility_audit.tsv"), audit)
write(os.path.join(OUT, "lncRNA_same_strand_eligible_pairs.tsv"),
      [dict(model_a=a, model_b=b, overlap_bp=bp,
            roles=f"{'competitor_only' if a in REPEAT_REG else 'ordinary'}+{'competitor_only' if b in REPEAT_REG else 'ordinary'}")
       for a, b, bp in pairs] or [dict(model_a="", model_b="", overlap_bp="", roles="")])
write(os.path.join(OUT, "bin100_master_annotation.tsv"), bins)

# summary
A = audit
n = lambda f: sum(1 for r in A if f(r))
summ = collections.OrderedDict()
summ["catalogue_lncRNA_models"] = len(A)
summ["retained_eligible"] = len(elig)
summ["  of which ordinary"] = len(elig - REPEAT_REG)
summ["  of which competitor_only (rDNA_repeat_associated)"] = len(elig & REPEAT_REG)
summ["  of which partial_rDNA flagged"] = len(elig & PARTIAL_REG)
summ["excluded_total"] = len(A) - len(elig)
for c in SCREEN:
    summ[f"excluded: same-strand {c} (any)"] = n(lambda r, c=c: r[f"{c}_same_n"] > 0)
    summ[f"excluded: primary reason same-strand {c}"] = n(lambda r, c=c: r["primary_exclusion_reason"] == f"same_strand_{c}")
summ["excluded: same-strand PCG (CDS or UTR, any)"] = n(lambda r: r["PCG_CDS_same_n"] + r["PCG_UTR_same_n"] > 0)
summ["excluded: same-strand tRNA_HC only (no same-strand PCG)"] = n(
    lambda r: r["tRNA_HC_same_n"] > 0 and r["PCG_CDS_same_n"] + r["PCG_UTR_same_n"] == 0)
summ["excluded: same-strand rRNA/OTHER_NC only (no same-strand PCG or tRNA_HC)"] = n(
    lambda r: (r["rRNA_same_n"] + r["OTHER_NC_same_n"]) > 0 and r["PCG_CDS_same_n"] + r["PCG_UTR_same_n"] + r["tRNA_HC_same_n"] == 0)
summ["retained despite opposite-strand overlap (any screen class)"] = n(lambda r: r["eligible"] == "yes" and r["opposite_strand_overlap_classes"])
for c in SCREEN:
    summ[f"  retained with opposite-strand {c}"] = n(lambda r, c=c: r["eligible"] == "yes" and r[f"{c}_opp_n"] > 0)
summ["retained with questionable tRNA (B9J08_001629) overlap, either strand"] = n(
    lambda r: r["eligible"] == "yes" and (r["tRNA_questionable_same_n"] + r["tRNA_questionable_opp_n"]) > 0)
summ["eligible models in >=1 same-strand eligible lncRNA pair"] = len(in_pair)
summ["same-strand eligible lncRNA pairs"] = len(pairs)
summ["old_vs_new: shared"] = n(lambda r: r["old_vs_new"] == "shared")
summ["old_vs_new: newly_recovered"] = n(lambda r: r["old_vs_new"] == "newly_recovered")
summ["old_vs_new: newly_excluded"] = n(lambda r: r["old_vs_new"] == "newly_excluded")
summ["old_vs_new: excluded_both"] = n(lambda r: r["old_vs_new"] == "excluded_both")
for cls in ("PCG_CDS", "lncRNA"):
    bs = [b for b in bins if b["rna_class"] == cls]
    ms = {b["model_id"] for b in bs}
    summ[f"bins {cls}: models"] = len(ms)
    summ[f"bins {cls}: bins"] = len(bs)
    summ[f"bins {cls}: total nt"] = sum(b["bin_length"] for b in bs)
    summ[f"bins {cls}: single-bin models (<=100 nt)"] = len({b["model_id"] for b in bs if b["single_bin_no_positional_contrast"]})
    summ[f"bins {cls}: shorter final bins"] = sum(1 for b in bs if b["final_bin_shorter"])
    summ[f"bins {cls}: competitor_only bins"] = sum(1 for b in bs if b["role"] == "competitor_only")
    summ[f"bins {cls}: bins with same-strand overlap to another binned model"] = sum(1 for b in bs if b["same_strand_binned_model_overlap"])
    summ[f"bins {cls}: bins overlapping HC tRNA same strand"] = sum(1 for b in bs if b["HC_tRNA_same_strand"])
    summ[f"bins {cls}: bins overlapping HC tRNA opposite strand"] = sum(1 for b in bs if b["HC_tRNA_opposite_strand"])
    summ[f"bins {cls}: bins with GANNOT"] = sum(1 for b in bs if b["GANNOT"])
    summ[f"bins {cls}: bins with rDNA flag"] = sum(1 for b in bs if b["rdna_flag"])
with open(os.path.join(OUT, "summary.tsv"), "w") as o:
    for k, v in summ.items():
        o.write(f"{k}\t{v}\n")
    for g_ in gates:
        o.write("\t".join(g_) + "\n")
with open(os.path.join(OUT, "outputs.md5"), "w") as o:
    for f in sorted(os.listdir(OUT)):
        if f != "outputs.md5":
            o.write(f"{md5(os.path.join(OUT, f))}  {f}\n")
print(open(os.path.join(OUT, "summary.tsv")).read())
sys.exit(0 if all(g_[1] == "PASS" for g_ in gates) else 1)
