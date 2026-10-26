#!/usr/bin/env python3
"""
codon_scorer.py  (v2 Step 6C (3), 2026-10-12; r2 2026-10-12)

r2 (Justine 2026-10-12, Q2 = a): optional --primary_only filter. When ON, a protein_coding CDS
record is kept only if its Ensembl header location is chromosome:GRCh38:<1-22, X, Y, MT>:...;
records on alternative haplotypes / patches (CHR_*) and scaffolds are dropped BEFORE the
longest-CDS-per-gene choice. When OFF (default) behaviour and every output file are identical to r1.

Human gene codon-bias composite score, copied from v1
EV_RNA_analysis/10_Codon_bias_human_pathway_link/scripts/20_human_gene_codon_scoring.py
with ONE change: the codon set is a parameter (--codons) instead of hard-coded.
Everything else (protein_coding filter, longest CDS per gene, length filters,
last codon dropped, >=5 codons per aa, >=3 aa per gene, composite = mean,
rounding, ranking, top 10%) is unchanged, so the v1 codon set must reproduce
v1 human_gene_codon_scores.tsv byte for byte.

Usage
    codon_scorer.py --fasta Homo_sapiens.GRCh38.cds.all.fa[.gz] \
                    --codons "Phe:UUC;Lys:AAG;His:CAC;Gln:CAG;Ser:UCU,UCA;Ala:GCA" \
                    --outdir DIR [--primary_only]

The functions load_gene_codon_counts() and score_genes() are imported by
make_null_foregrounds.py so the null sets use exactly the same code path.
"""
import argparse
import gzip
import os
from collections import defaultdict

CODON_TABLE = {
    "UUU": "Phe", "UUC": "Phe", "UUA": "Leu", "UUG": "Leu",
    "CUU": "Leu", "CUC": "Leu", "CUA": "Leu", "CUG": "Leu",  # human standard code
    "AUU": "Ile", "AUC": "Ile", "AUA": "Ile", "AUG": "Met",
    "GUU": "Val", "GUC": "Val", "GUA": "Val", "GUG": "Val",
    "UCU": "Ser", "UCC": "Ser", "UCA": "Ser", "UCG": "Ser",
    "CCU": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACU": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCU": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "UAU": "Tyr", "UAC": "Tyr", "UAA": "STOP", "UAG": "STOP",
    "CAU": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "AAU": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "GAU": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "UGU": "Cys", "UGC": "Cys", "UGA": "STOP", "UGG": "Trp",
    "CGU": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGU": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GGU": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}

PRIMARY_SEQS = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}

MIN_CODONS_PER_AA = 5
MIN_AAS_REQUIRED = 3


def parse_codon_set(spec):
    """'Phe:UUC;Ser:UCU,UCA' -> {'Phe': {'UUC'}, 'Ser': {'UCU','UCA'}} (order kept)."""
    out = {}
    for item in spec.split(";"):
        item = item.strip()
        if not item:
            continue
        aa, codons = item.split(":")
        cset = {c.strip().upper().replace("T", "U") for c in codons.split(",")}
        for c in cset:
            if CODON_TABLE.get(c) != aa:
                raise ValueError(f"codon {c} is not {aa} in the standard code")
        if aa in out:
            raise ValueError(f"amino acid {aa} given twice")
        out[aa] = cset
    return out


def header_location(line):
    """Ensembl cds header: '>ENST... cds chromosome:GRCh38:6:start:end:strand gene:...' -> 3rd token."""
    parts = line[1:].split()
    return parts[2] if len(parts) > 2 else ""


def is_primary(location):
    f = location.split(":")
    return len(f) >= 3 and f[0] == "chromosome" and f[1] == "GRCh38" and f[2] in PRIMARY_SEQS


def parse_header(line):
    parts = line[1:].split()
    transcript_id = parts[0]
    gene_id = gene_biotype = gene_symbol = None
    for p in parts[1:]:
        if p.startswith("gene:"):
            gene_id = p.split(":", 1)[1].split(".")[0]
        elif p.startswith("gene_biotype:"):
            gene_biotype = p.split(":", 1)[1]
        elif p.startswith("gene_symbol:"):
            gene_symbol = p.split(":", 1)[1]
    return transcript_id, gene_id, gene_biotype, gene_symbol


def _open(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def load_fasta_filtered(path, with_location=False):
    header = None
    loc = None
    parts = []
    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    tid, gid, biotype, sym = header
                    if biotype == "protein_coding" and gid:
                        yield (tid, gid, sym, "".join(parts), loc) if with_location else (tid, gid, sym, "".join(parts))
                header = parse_header(line)
                loc = header_location(line)
                parts = []
            else:
                parts.append(line)
        if header is not None:
            tid, gid, biotype, sym = header
            if biotype == "protein_coding" and gid:
                yield (tid, gid, sym, "".join(parts), loc) if with_location else (tid, gid, sym, "".join(parts))


def load_gene_codon_counts(fasta, primary_only=False):
    """Return (genes, stats). genes = list of (gene_id, symbol, aa_counts) in v1 order
    (dict insertion order of first appearance); aa_counts[aa][codon] = n.
    r2: primary_only=True drops protein_coding records whose location is not a primary chromosome;
    stats["dropped"] lists them (tid, gid, symbol, location)."""
    by_gene = {}
    n_total = 0
    dropped = []
    primary_genes = set()
    for tid, gid, sym, seq, loc in load_fasta_filtered(fasta, with_location=True):
        n_total += 1
        if primary_only:
            if not is_primary(loc):
                dropped.append((tid, gid, sym or "", loc))
                continue
            primary_genes.add(gid)
        if gid not in by_gene or len(seq) > len(by_gene[gid][1]):
            by_gene[gid] = (sym, seq)
    genes = []
    n_dropped_short = 0
    for gid, (sym, seq) in by_gene.items():
        seq = seq.upper().replace("T", "U")
        if len(seq) < 30 or len(seq) % 3 != 0:
            n_dropped_short += 1
            continue
        codons = [seq[i:i+3] for i in range(0, len(seq) - 3, 3)]  # drop last codon (stop)
        aa_counts = defaultdict(lambda: defaultdict(int))
        for c in codons:
            aa = CODON_TABLE.get(c)
            if aa and aa != "STOP":
                aa_counts[aa][c] += 1
        genes.append((gid, sym, aa_counts))
    stats = {"n_total": n_total, "n_genes": len(by_gene), "n_dropped_short": n_dropped_short,
             "primary_only": primary_only, "dropped": dropped,
             "n_dropped_genes_entirely": len({d[1] for d in dropped} - primary_genes),
             "n_genes_mixed_primary_nonprimary": len({d[1] for d in dropped} & primary_genes)}
    return genes, stats


def score_genes(genes, target_codons):
    """v1 scoring for one codon set. Returns rows sorted as v1 (stable sort on rounded score)."""
    rows = []
    for gid, sym, aa_counts in genes:
        per_aa_pct = {}
        for aa, target_set in target_codons.items():
            counts = aa_counts.get(aa, {})
            total_aa = sum(counts.values())
            if total_aa < MIN_CODONS_PER_AA:
                continue
            target_count = sum(counts.get(c, 0) for c in target_set)
            per_aa_pct[aa] = target_count / total_aa * 100
        if len(per_aa_pct) < MIN_AAS_REQUIRED:
            continue
        composite = sum(per_aa_pct.values()) / len(per_aa_pct)
        row = {"gene_id": gid, "symbol": sym or "", "composite_score": round(composite, 2),
               "n_aas_scored": len(per_aa_pct)}
        for aa in target_codons:
            row[f"{aa}_pct"] = round(per_aa_pct.get(aa, float("nan")), 1) if aa in per_aa_pct else ""
        rows.append(row)
    rows.sort(key=lambda r: -r["composite_score"])
    return rows


def n_foreground(rows):
    return max(1, round(len(rows) * 0.10))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--codons", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--primary_only", action="store_true",
                    help="r2: keep only CDS on chromosome:GRCh38:1-22/X/Y/MT (default OFF = r1 behaviour)")
    a = ap.parse_args()
    target = parse_codon_set(a.codons)
    os.makedirs(a.outdir, exist_ok=True)
    genes, st = load_gene_codon_counts(a.fasta, primary_only=a.primary_only)
    rows = score_genes(genes, target)

    cols = ["gene_id", "symbol", "composite_score", "n_aas_scored"] + [f"{aa}_pct" for aa in target]
    with open(os.path.join(a.outdir, "human_gene_codon_scores.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    n_fg = n_foreground(rows)
    with open(os.path.join(a.outdir, "foreground_top10pct.txt"), "w", encoding="utf-8") as fh:
        for r in rows[:n_fg]:
            fh.write(r["gene_id"] + "\n")
    with open(os.path.join(a.outdir, "background_all_scored.txt"), "w", encoding="utf-8") as fh:
        for r in rows:
            fh.write(r["gene_id"] + "\n")
    with open(os.path.join(a.outdir, "scoring_QC.txt"), "w", encoding="utf-8") as fh:
        fh.write("Human gene codon-bias composite scoring QC\n")
        fh.write("=" * 60 + "\n")
        fh.write(f"codon set: {a.codons}\n")
        fh.write(f"protein_coding transcripts scanned: {st['n_total']}\n")
        if st["primary_only"]:  # r2 lines; absent when the filter is OFF so the QC file stays r1-identical
            fh.write("primary-assembly filter: ON (chromosome:GRCh38:1-22,X,Y,MT)\n")
            fh.write(f"protein_coding transcripts dropped (non-primary location): {len(st['dropped'])}\n")
            fh.write(f"genes with no primary-chromosome CDS (dropped entirely): {st['n_dropped_genes_entirely']}\n")
            fh.write(f"genes with CDS on both primary and non-primary sequences: {st['n_genes_mixed_primary_nonprimary']}\n")
        fh.write(f"unique protein_coding genes (longest CDS kept): {st['n_genes']}\n")
        fh.write(f"genes dropped (too short/bad length): {st['n_dropped_short']}\n")
        fh.write(f"genes scored: {len(rows)}\n")
        fh.write(f"foreground (top 10%): {n_fg}\n")
        fh.write(f"background (all scored): {len(rows)}\n\n")
        fh.write("Top 15 genes by composite score:\n")
        for r in rows[:15]:
            fh.write(f"  {r['symbol']:<15}{r['gene_id']:<20}composite={r['composite_score']}\n")
        if rows:
            fh.write(f"\nMedian composite score: {sorted(r['composite_score'] for r in rows)[len(rows)//2]:.1f}\n")
    if st["primary_only"]:
        with open(os.path.join(a.outdir, "nonprimary_dropped_transcripts.tsv"), "w", encoding="utf-8") as fh:
            fh.write("transcript_id\tgene_id\tsymbol\tlocation\n")
            for d in st["dropped"]:
                fh.write("\t".join(d) + "\n")
    print(f"genes scored {len(rows)}; foreground {n_fg}")


if __name__ == "__main__":
    main()
