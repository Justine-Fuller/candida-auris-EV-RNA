#!/usr/bin/env python3
"""
compare_runs.py  (copy of compare_r1_r2.py, used 2026-10-13 for r2 vs r3; labels "r1"/"r2" in outputs = first/second run)
Side-by-side r1 (all Ensembl cds.all CDS) vs r2 (primary assembly only). Read-only on both runs.
Outputs: counts_r1_vs_r2.tsv, terms_sig_either_run.tsv, focus_MHC_antigen_terms.tsv, top_specific_r2.tsv
Args: <r1_run_dir> <r2_run_dir> <outdir>
"""
import os
import re
import sys

ONTS = [("BP", "GO_BP"), ("MF", "GO_MF"), ("CC", "GO_CC"), ("KEGG", "KEGG")]
FOCUS = re.compile(r"MHC|HLA|antigen|allograft|graft-versus-host|autoimmune thyroid|type I diabetes", re.I)


def read_tsv(p):
    with open(p) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, l.rstrip("\n").split("\t"))) for l in fh]


def load(run):
    allt, spec = {}, {}
    for o, tag in ONTS:
        for r in read_tsv(os.path.join(run, "S3_enrichment_real", f"{tag}_all.tsv")):
            allt[(o, r["ID"])] = r
    for r in read_tsv(os.path.join(run, "S6_summary", "specificity_tracked_terms.tsv")):
        spec[(r["ontology"], r["ID"])] = r
    return allt, spec


def main():
    r1d, r2d, od = sys.argv[1:4]
    os.makedirs(od, exist_ok=True)
    a1, s1 = load(r1d)
    a2, s2 = load(r2d)
    sig = lambda a, k: k in a and float(a[k]["p.adjust"]) < 0.05
    spc = lambda s, a, k: sig(a, k) and k in s and s[k]["specific"] == "yes"
    keys = sorted(set(a1) | set(a2))
    with open(os.path.join(od, "counts_r1_vs_r2.tsv"), "w") as fh:
        cols = ["ontology", "r1_sig", "r1_specific", "r2_sig", "r2_specific", "sig_both", "sig_r1_only", "sig_r2_only",
                "specific_both", "specific_r1_only", "specific_r2_only"]
        fh.write("\t".join(cols) + "\n")
        tot = [0] * 10
        for o, _ in ONTS + [("Total", None)]:
            ks = [k for k in keys if k[0] == o] if o != "Total" else keys
            v = [sum(sig(a1, k) for k in ks), sum(spc(s1, a1, k) for k in ks), sum(sig(a2, k) for k in ks),
                 sum(spc(s2, a2, k) for k in ks), sum(sig(a1, k) and sig(a2, k) for k in ks),
                 sum(sig(a1, k) and not sig(a2, k) for k in ks), sum(sig(a2, k) and not sig(a1, k) for k in ks),
                 sum(spc(s1, a1, k) and spc(s2, a2, k) for k in ks), sum(spc(s1, a1, k) and not spc(s2, a2, k) for k in ks),
                 sum(spc(s2, a2, k) and not spc(s1, a1, k) for k in ks)]
            fh.write(o + "\t" + "\t".join(map(str, v)) + "\n")

    def row(k):
        d = (a2.get(k) or a1.get(k))["Description"]
        out = [k[0], k[1], d]
        for a, s in ((a1, s1), (a2, s2)):
            r = a.get(k)
            out += [r["p.adjust"] if r else "NA(no fg gene)", r["Count"] if r else "0",
                    s[k]["empirical_fraction"] if k in s else "NA(not tracked)",
                    "yes" if sig(a, k) else "no", "yes" if spc(s, a, k) else "no"]
        return out
    head = ["ontology", "ID", "Description", "r1_q", "r1_count", "r1_fraction", "r1_sig", "r1_specific",
            "r2_q", "r2_count", "r2_fraction", "r2_sig", "r2_specific"]
    with open(os.path.join(od, "terms_sig_either_run.tsv"), "w") as fh:
        fh.write("\t".join(head) + "\n")
        for k in keys:
            if sig(a1, k) or sig(a2, k):
                fh.write("\t".join(row(k)) + "\n")
    with open(os.path.join(od, "focus_MHC_antigen_terms.tsv"), "w") as fh:
        fh.write("\t".join(head) + "\n")
        for k in keys:
            d = (a2.get(k) or a1.get(k))["Description"]
            if FOCUS.search(d) and (sig(a1, k) or sig(a2, k) or k in s1 or k in s2):
                fh.write("\t".join(row(k)) + "\n")
    with open(os.path.join(od, "top_specific_r2.tsv"), "w") as fh:
        fh.write("\t".join(head) + "\n")
        for o, _ in ONTS:
            ks = sorted([k for k in keys if k[0] == o and spc(s2, a2, k)], key=lambda k: float(a2[k]["pvalue"]))
            for k in ks:
                fh.write("\t".join(row(k)) + "\n")
    print(open(os.path.join(od, "counts_r1_vs_r2.tsv")).read())


if __name__ == "__main__":
    main()
