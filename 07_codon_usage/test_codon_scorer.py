#!/usr/bin/env python3
"""Known-answer test for codon_scorer.py on a tiny synthetic CDS set (hand-computed answers).

Codon set: Phe:UUC; Lys:AAG; His:CAC; Gln:CAG; Ser:UCU,UCA
  G1  ATG UUCx5 UUUx5 AAGx6 CACx1 CATx4 TAA
      Phe 5/10=50.0  Lys 6/6=100.0  His 1/5=20.0  (Gln, Ser absent) -> 3 aa, composite 56.67
  G2  two transcripts; the longer one is kept:
      ATG TTTx5 AAAx5 CAAx5 TCTx3 TCAx2 TCCx5 TAA
      Phe 0.0 Lys 0.0 Gln 0.0 Ser 5/10=50.0 -> 4 aa, composite 12.5
  G3  biotype lncRNA -> excluded
  G4  ATG TTCx5 AAGx5 CACx4 TAA -> His has 4 (<5) -> 2 aa -> dropped
  G5  length not divisible by 3 -> dropped (counted as short/bad length)
  G6  24 nt (<30) -> dropped (counted as short/bad length)
  G7  ATG TTCx6 AAGx5 CACx5 TAA -> Phe 6/6, Lys 5/5, His 5/5 -> composite 100.0
Expected order: G7 (100.0), G1 (56.67), G2 (12.5); foreground = max(1, round(0.3)) = 1 -> G7.

r2 part B (--primary_only), codon set Phe:UUC; Lys:AAG; His:CAC; full Ensembl-style locations:
  P1  chromosome:GRCh38:6         ATG TTCx6 AAGx5 CACx5 TAA                   -> 100.0
  P2  chromosome:GRCh38:CHR_HSCHR6_MHC_COX_CTG1 (haplotype)
                                  ATG TTCx5 TTTx5 AAGx6 CAC CATx4 TAA         -> 56.67 (dropped when ON)
  P3  scaffold:GRCh38:KI270728.1  ATG TTCx5 AAGx5 CACx5 TAA                   -> 100.0 (dropped when ON)
  P4  T4a chromosome:GRCh38:X     ATG TTTx5 AAGx5 CACx5 TAA (shorter)         -> 66.67
      T4b chromosome:GRCh38:CHR_HG2_PATCH (longer) ATG TTCx10 AAGx5 CACx5 TAA -> 100.0
      OFF: longest = T4b -> 100.0;  ON: T4b dropped -> T4a -> 66.67
  P5  chromosome:GRCh38:MT        ATG TTCx5 TTTx5 AAGx5 CACx5 TAA             -> 83.33
  P6  chromosome:GRCh38:Y         ATG TTTx5 AAAx5 CATx5 TAA                   -> 0.0
  OFF order (stable sort, first-appearance order on ties): P1, P3, P4 (100), P5, P2, P6
  ON  order: P1 (100), P5 (83.33), P4 (66.67), P6 (0.0); dropped transcripts 3 (P2, P3, T4b);
      genes dropped entirely 2 (P2, P3); genes on both primary and non-primary 1 (P4).
"""
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))


def rec(tid, gid, biotype, sym, seq):
    return f">{tid} cds chromosome:X gene:{gid}.1 gene_biotype:{biotype} transcript_biotype:{biotype} gene_symbol:{sym}\n{seq}\n"


def main():
    fa = ""
    fa += rec("T1", "G1", "protein_coding", "A1", "ATG" + "TTC"*5 + "TTT"*5 + "AAG"*6 + "CAC" + "CAT"*4 + "TAA")
    fa += rec("T2s", "G2", "protein_coding", "A2", "ATG" + "TTC"*5 + "AAG"*5 + "CAG"*5 + "TAA")
    fa += rec("T2l", "G2", "protein_coding", "A2", "ATG" + "TTT"*5 + "AAA"*5 + "CAA"*5 + "TCT"*3 + "TCA"*2 + "TCC"*5 + "TAA")
    fa += rec("T3", "G3", "lncRNA", "A3", "ATG" + "TTC"*10 + "AAG"*10 + "CAC"*10 + "TAA")
    fa += rec("T4", "G4", "protein_coding", "A4", "ATG" + "TTC"*5 + "AAG"*5 + "CAC"*4 + "TAA")
    fa += rec("T5", "G5", "protein_coding", "A5", "ATG" + "TTC"*10 + "AAG"*10 + "CAC"*10 + "TAAG")
    fa += rec("T6", "G6", "protein_coding", "A6", "ATG" + "TTC"*6 + "TAA")
    fa += rec("T7", "G7", "protein_coding", "A7", "ATG" + "TTC"*6 + "AAG"*5 + "CAC"*5 + "TAA")
    expected = (
        "gene_id\tsymbol\tcomposite_score\tn_aas_scored\tPhe_pct\tLys_pct\tHis_pct\tGln_pct\tSer_pct\n"
        "G7\tA7\t100.0\t3\t100.0\t100.0\t100.0\t\t\n"
        "G1\tA1\t56.67\t3\t50.0\t100.0\t20.0\t\t\n"
        "G2\tA2\t12.5\t4\t0.0\t0.0\t\t0.0\t50.0\n"
    )
    with tempfile.TemporaryDirectory() as d:
        fp = os.path.join(d, "syn.fa")
        open(fp, "w").write(fa)
        out = os.path.join(d, "out")
        subprocess.run([sys.executable, os.path.join(HERE, "codon_scorer.py"), "--fasta", fp,
                        "--codons", "Phe:UUC;Lys:AAG;His:CAC;Gln:CAG;Ser:UCU,UCA", "--outdir", out], check=True)
        got = open(os.path.join(out, "human_gene_codon_scores.tsv")).read()
        fg = open(os.path.join(out, "foreground_top10pct.txt")).read()
        bg = open(os.path.join(out, "background_all_scored.txt")).read()
        qc = open(os.path.join(out, "scoring_QC.txt")).read()
    checks = {
        "scores table": got == expected,
        "foreground": fg == "G7\n",
        "background": bg == "G7\nG1\nG2\n",
        "QC transcripts scanned 7": "protein_coding transcripts scanned: 7" in qc,
        "QC unique genes 6": "unique protein_coding genes (longest CDS kept): 6" in qc,
        "QC dropped short/bad 2": "genes dropped (too short/bad length): 2" in qc,
    }
    # codon-set parser rejects a wrong amino-acid/codon pair
    sys.path.insert(0, HERE)
    import codon_scorer
    try:
        codon_scorer.parse_codon_set("Phe:UUA")
        checks["parser rejects Phe:UUA"] = False
    except ValueError:
        checks["parser rejects Phe:UUA"] = True
    # ---- r2 part B: primary-assembly filter ----
    def recl(tid, gid, loc, seq):
        return (f">{tid} cds {loc}:1:100:1 gene:{gid}.1 gene_biotype:protein_coding "
                f"transcript_biotype:protein_coding gene_symbol:S{gid}\n{seq}\n")
    fb = ""
    fb += recl("U1", "P1", "chromosome:GRCh38:6", "ATG" + "TTC"*6 + "AAG"*5 + "CAC"*5 + "TAA")
    fb += recl("U2", "P2", "chromosome:GRCh38:CHR_HSCHR6_MHC_COX_CTG1", "ATG" + "TTC"*5 + "TTT"*5 + "AAG"*6 + "CAC" + "CAT"*4 + "TAA")
    fb += recl("U3", "P3", "scaffold:GRCh38:KI270728.1", "ATG" + "TTC"*5 + "AAG"*5 + "CAC"*5 + "TAA")
    fb += recl("U4a", "P4", "chromosome:GRCh38:X", "ATG" + "TTT"*5 + "AAG"*5 + "CAC"*5 + "TAA")
    fb += recl("U4b", "P4", "chromosome:GRCh38:CHR_HG2_PATCH", "ATG" + "TTC"*10 + "AAG"*5 + "CAC"*5 + "TAA")
    fb += recl("U5", "P5", "chromosome:GRCh38:MT", "ATG" + "TTC"*5 + "TTT"*5 + "AAG"*5 + "CAC"*5 + "TAA")
    fb += recl("U6", "P6", "chromosome:GRCh38:Y", "ATG" + "TTT"*5 + "AAA"*5 + "CAT"*5 + "TAA")
    hdr = "gene_id\tsymbol\tcomposite_score\tn_aas_scored\tPhe_pct\tLys_pct\tHis_pct\n"
    exp_off = hdr + ("P1\tSP1\t100.0\t3\t100.0\t100.0\t100.0\n"
                     "P3\tSP3\t100.0\t3\t100.0\t100.0\t100.0\n"
                     "P4\tSP4\t100.0\t3\t100.0\t100.0\t100.0\n"
                     "P5\tSP5\t83.33\t3\t50.0\t100.0\t100.0\n"
                     "P2\tSP2\t56.67\t3\t50.0\t100.0\t20.0\n"
                     "P6\tSP6\t0.0\t3\t0.0\t0.0\t0.0\n")
    exp_on = hdr + ("P1\tSP1\t100.0\t3\t100.0\t100.0\t100.0\n"
                    "P5\tSP5\t83.33\t3\t50.0\t100.0\t100.0\n"
                    "P4\tSP4\t66.67\t3\t0.0\t100.0\t100.0\n"
                    "P6\tSP6\t0.0\t3\t0.0\t0.0\t0.0\n")
    with tempfile.TemporaryDirectory() as d:
        fp = os.path.join(d, "synB.fa")
        open(fp, "w").write(fb)
        res = {}
        for mode, extra in (("off", []), ("on", ["--primary_only"])):
            out = os.path.join(d, mode)
            subprocess.run([sys.executable, os.path.join(HERE, "codon_scorer.py"), "--fasta", fp,
                            "--codons", "Phe:UUC;Lys:AAG;His:CAC", "--outdir", out] + extra, check=True)
            res[mode] = {f: (open(os.path.join(out, f)).read() if os.path.exists(os.path.join(out, f)) else None)
                         for f in ("human_gene_codon_scores.tsv", "scoring_QC.txt", "nonprimary_dropped_transcripts.tsv")}
    checks["B filter OFF scores (non-primary kept, longest T4b)"] = res["off"]["human_gene_codon_scores.tsv"] == exp_off
    checks["B filter OFF writes no dropped list / no filter QC lines"] = (
        res["off"]["nonprimary_dropped_transcripts.tsv"] is None and "primary-assembly" not in res["off"]["scoring_QC.txt"])
    checks["B filter ON scores (P2, P3 gone; P4 from T4a)"] = res["on"]["human_gene_codon_scores.tsv"] == exp_on
    qon = res["on"]["scoring_QC.txt"] or ""
    checks["B QC dropped transcripts 3"] = "protein_coding transcripts dropped (non-primary location): 3" in qon
    checks["B QC genes dropped entirely 2"] = "genes with no primary-chromosome CDS (dropped entirely): 2" in qon
    checks["B QC mixed genes 1"] = "genes with CDS on both primary and non-primary sequences: 1" in qon
    dl = (res["on"]["nonprimary_dropped_transcripts.tsv"] or "").splitlines()
    checks["B dropped list = U2, U3, U4b"] = [x.split("\t")[0] for x in dl[1:]] == ["U2", "U3", "U4b"]
    checks["B is_primary edge cases"] = (
        codon_scorer.is_primary("chromosome:GRCh38:22:1:2:1") and codon_scorer.is_primary("chromosome:GRCh38:MT:1:2:1")
        and not codon_scorer.is_primary("chromosome:GRCh38:23:1:2:1")
        and not codon_scorer.is_primary("chromosome:GRCh38:CHR_HSCHR6_MHC_COX_CTG1:1:2:1")
        and not codon_scorer.is_primary("scaffold:GRCh38:KI270728.1:1:2:1")
        and not codon_scorer.is_primary("chromosome:GRCh37:6:1:2:1") and not codon_scorer.is_primary(""))
    ok = True
    for k, v in checks.items():
        print(f"{'PASS' if v else 'FAIL'}\t{k}")
        ok &= v
    if not ok:
        print("GOT:\n" + got)
    print("KNOWN-ANSWER TEST:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
