#!/usr/bin/env python3
"""
================================================================================
 build_pcg_annotation_table.py
================================================================================

SCRIPT TITLE
    Master C. auris Protein-Coding Gene (PCG) Functional Annotation Builder

PURPOSE
    Build ONE master functional-annotation table -- one row per Candida
    auris B8441 protein-coding gene -- by combining the outputs of six
    independently run annotation tools: eggNOG-mapper, InterProScan,
    OrthoFinder, SignalP 6.0, DeepTMHMM and PredGPI.

PIPELINE CONTEXT (C. auris caspofungin EV/cell RNA-seq analysis pipeline)
    Step 1-6  : read processing, alignment, DESeq2 differential expression
                (cells and EVs), strain comparison
    Step 7    : join genomic lncRNA-nearest-PCG neighbour pairs to
                strain-specific DESeq2 results
    Step 8    : classify each lncRNA-PCG pair's co-regulation relationship
                (co-up/co-down/opposite/etc.) per strain
    Step 9    : hierarchically rank the strongest concordantly regulated
                lncRNA-PCG neighbour pairs
  > Step 10   : (THIS SCRIPT) build a master functional-annotation table
                for every C. auris PCG -- run independently of the
                lncRNA-PCG pairing work above; no join to those tables
                happens here
    Step 11+  : (not part of this script) join this annotation table onto
                the ranked lncRNA-PCG pairs from Step 9, to add
                functional/orthology context to candidate cis-associated
                pairs

BIOLOGICAL OBJECTIVE
    The DESeq2 + lncRNA-PCG pairing steps identify WHICH genes respond to
    caspofungin and WHICH neighbouring lncRNA-PCG pairs move together, but
    say nothing about what those genes actually do. This step attaches
    orthogonal functional evidence -- known function/domains (eggNOG,
    InterProScan), cross-species orthology to well-studied yeasts
    (OrthoFinder), and membrane/secretion topology (SignalP, DeepTMHMM,
    PredGPI) -- to every PCG in the genome, so downstream candidate
    lncRNA-PCG pairs can eventually be interpreted biologically (e.g.
    "this candidate pair's PCG is a predicted GPI-anchored cell-wall
    protein orthologous to a characterised S. cerevisiae gene").

INPUT FILES / DIRECTORIES (paths relative to --project-root, default
/scratch/flljus002; see build_config() for the exact join logic)
    Backbone (defines the full PCG gene set, one row per gene):
        orthofinder/Proteomes/FungiDB-68_CaurisB8441_AnnotatedProteins.fasta

    eggNOG-mapper (full proteome; ran on ALL PCGs):
        eggnog_full_results/full_proteome_annotation.emapper.annotations

    InterProScan (subset: significant up/down PCGs only):
        interproscan_results/up/up_interproscan.tsv
        interproscan_results/down/down_interproscan.tsv

    OrthoFinder (full proteome; 4-species comparison):
        orthofinder/Proteomes/OrthoFinder/Results_Four_fungi_comparison/
            Orthogroups/Orthogroups.tsv
            Orthogroups/Orthogroups_UnassignedGenes.tsv
            Orthologues/Orthologues_*Cauris*/*__v__*.tsv   (3 pairwise
                files, discovered by species-name pattern -- not
                hard-coded, see SPECIES_KEYWORDS)

    SignalP 6.0 (subset: significant up/down PCGs only):
        signalp_up_results/up_prediction_results.txt
        signalp_down_results/down_prediction_results.txt

    DeepTMHMM (subset: significant up/down PCGs only):
        deeptmhmm_results/up/biolib_results/TMRs.gff3
        deeptmhmm_results/up/biolib_results/predicted_topologies.3line
        deeptmhmm_results/down/biolib_results/TMRs.gff3
        deeptmhmm_results/down/biolib_results/predicted_topologies.3line

    PredGPI (subset: significant up/down PCGs only):
        predgpi_results/up_predgpi.json
        predgpi_results/down_predgpi.json

OUTPUT FILES (written to --output-dir, default
<project-root>/PCG_functional_annotation/)
    PCG_functional_annotation.tsv       master table, tab-delimited
    PCG_functional_annotation.csv       master table, comma-delimited
    eggnog_collapsed.tsv                intermediate: one row per gene
    interpro_collapsed.tsv              intermediate: one row per gene
    orthofinder_collapsed.tsv           intermediate: one row per gene
    signalp_collapsed.tsv               intermediate: one row per gene
    deeptmhmm_collapsed.tsv             intermediate: one row per gene
    predgpi_collapsed.tsv               intermediate: one row per gene
    PCG_functional_annotation_QC.txt    full run log = QC report (every
                                         INFO/WARNING emitted during the
                                         run is captured here automatically
                                         via setup_logging())

EXTERNAL SOFTWARE / DATABASES (already run upstream; NOT invoked by this
script -- this script only parses their pre-computed output files)
    eggNOG-mapper v2.1.15   DIAMOND search against the eggNOG 5 database
    InterProScan            CDD, Coils, FunFam, Gene3D, Hamap, MobiDBLite,
                             NCBIfam, PANTHER, Pfam, PIRSF, PRINTS,
                             ProSitePatterns, ProSiteProfiles, SFLD, SMART,
                             SUPERFAMILY member databases, plus GO/pathway
                             mapping
    OrthoFinder v3.1.5       diamond search, famsa alignment, fasttree;
                             4-way comparison: C. auris B8441,
                             C. albicans SC5314, S. cerevisiae S288C,
                             S. pombe 972h
    SignalP 6.0              Eukarya mode
    DeepTMHMM                BioLib-hosted
    PredGPI                  EBI-style JSON export

PYTHON PACKAGE DEPENDENCIES
    Developed and validated on Python 3.12.1 (standard-library features
    used -- pathlib, typing, `from __future__ import annotations` -- are
    compatible with Python 3.8+).
    pandas    the only third-party dependency (validated with pandas 3.0.5)
    Standard library only otherwise: argparse, csv, json, logging, os, re,
    sys, pathlib, typing

HIGH-LEVEL WORKFLOW (see build_master_table() for the orchestration code;
every step follows the same PARSE -> QC -> LEFT-MERGE onto the running
master -> QC pattern)
    0. BACKBONE     Read the C. auris proteome FASTA; every entry's own
                     'gene=' / 'transcript=' header fields (not a guessed
                     suffix) become gene_id / transcript_id / protein_id.
                     This -- not eggNOG's output -- is the master gene
                     list, because eggNOG silently OMITS genes with no
                     DIAMOND hit rather than emitting an empty row for
                     them.
    1. eggNOG       Parse the full-proteome annotation table; left-merge
                     onto the backbone on gene_id.
    2. InterProScan Parse+collapse the up/down subset TSVs (many
                     domain-hit rows per protein -> one row per gene; see
                     "DUPLICATE HANDLING" below); left-merge.
    3. OrthoFinder  Parse the full-proteome Orthogroup assignment AND the
                     three pairwise ortholog tables (C. albicans,
                     S. cerevisiae, S. pombe); left-merge.
    4. SignalP      Parse+collapse the up/down subset predictions;
                     left-merge.
    5. DeepTMHMM    Parse+collapse the up/down subset TM predictions
                     (TMRs.gff3 for counts/regions, .3line for the overall
                     topology class); left-merge.
    6. PredGPI      Parse+collapse the up/down subset GPI-anchor JSON
                     predictions; left-merge.
    7. FINAL QC, fixed column ordering, and save every TSV/CSV output.

QUALITY CONTROL PERFORMED (every check below either logs a PASS to the QC
report or raises an exception immediately -- nothing fails silently)
    - Every raw ID from every source is normalised to a B9J08_###### gene
      ID and verified to actually match that pattern (guards against
      silently treating an unrelated/malformed ID as a C. auris gene) --
      see normalize_id_series().
    - Duplicate raw IDs are reported for every source before collapsing --
      flagged as EXPECTED for InterProScan (many domain hits per protein)
      and as a WARNING for every other source (which should already be
      one-row-per-protein) -- see report_duplicates().
    - After collapsing, every per-source table is asserted to have exactly
      one row per gene_id -- see assert_one_row_per_gene().
    - Before and after every merge, the master table's row count is
      logged and compared; a merge that changes the row count, or that
      produces a duplicated gene_id, raises immediately -- see
      qc_before_merge() / qc_after_merge().
    - After every merge, coverage is reported: how many / how few PCGs
      received a non-null value from that source -- see
      report_source_coverage().
    - A DeepTMHMM-specific cross-check: the parsed count of membrane-
      spanning segments (TMhelix + Beta sheet) is verified against that
      protein's own self-reported 'Number of predicted TMRs' value -- see
      _parse_tmrs_gff3().
    - Final QC re-confirms the master table's row count still equals the
      backbone's gene count, and that gene_id has no duplicates or nulls,
      before any file is written.

EXPECTED OUTPUT DIMENSIONS
    One row per C. auris protein-coding gene -- 5,417 rows in the dataset
    this script was validated against (the full FungiDB-68
    CaurisB8441_AnnotatedProteins.fasta gene count), 36 columns. This row
    count must be IDENTICAL to the backbone gene count at every stage; the
    QC checks above exist specifically to guarantee that.

WHY LEFT JOINS, AND HOW DUPLICATE ANNOTATIONS ARE COLLAPSED
    Every per-source merge in build_master_table() is a LEFT join of the
    running master (backbone-derived, always one row per gene) onto each
    newly parsed source table. This is deliberate: eggNOG and OrthoFinder
    ran on the FULL proteome but eggNOG still silently omits no-hit genes,
    while InterProScan/SignalP/DeepTMHMM/PredGPI ran on SUBSETS
    (significant up/down PCGs only). A left join guarantees every PCG in
    the backbone survives every merge regardless of which source(s) have
    (or lack) data for it -- an inner join would silently shrink the table
    every time a source doesn't cover a gene, and this pipeline's QC
    explicitly forbids the master row count ever changing.

    Because a left join can only ever ADD columns (never rows) when the
    right-hand table is guaranteed one-row-per-key, every source table is
    collapsed to exactly one row per gene_id BEFORE it is merged.
    InterProScan/SignalP/DeepTMHMM/PredGPI are each reduced to
    one-row-per-protein via a groupby(gene_id) + per-column aggregation:
      - single-value fields use unique_join() -- sorted, deduplicated,
        ';'-joined;
      - fields that are themselves multi-value per row (e.g.
        InterProScan's '|'-joined GO terms) use split_and_collect(), which
        flattens every row's tokens before deduplicating, so no domain
        hit, GO term or pathway annotation is ever discarded -- only
        de-duplicated.
    This also transparently handles the rare genes (3 in the validated
    dataset) that appear in BOTH the up and down subset files: their rows
    are concatenated before grouping, so they still collapse to a single,
    complete row rather than being overwritten or duplicated.

This script does NOT:
    - rerun any annotation tool
    - rerun DESeq2
    - redo lncRNA-PCG genomic pairing
    - join anything to the lncRNA tables (that is a later step)

USAGE
    python3 build_pcg_annotation_table.py
    python3 build_pcg_annotation_table.py --project-root /scratch/flljus002 \
        --output-dir /scratch/flljus002/PCG_functional_annotation
================================================================================
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd

# =============================================================================
# CONFIGURATION -- all paths live here, and ONLY here. No function below
# should ever contain a literal file name; every path is passed in as an
# argument, ultimately sourced from this block (or its CLI overrides).
# =============================================================================


def build_config(project_root: Path, output_dir: Path) -> Dict[str, Path]:
    """
    Build the full set of input/output paths from a single project root and
    output directory. Centralising this in one function means the entire
    script can be repointed at a different run (e.g. a re-analysis with new
    tool outputs) by changing two paths, without hunting through the code.
    """
    orthofinder_results = (
        project_root
        / "orthofinder"
        / "Proteomes"
        / "OrthoFinder"
        / "Results_Four_fungi_comparison"
    )

    return {
        # --- backbone ---
        "fasta_backbone": project_root
        / "orthofinder"
        / "Proteomes"
        / "FungiDB-68_CaurisB8441_AnnotatedProteins.fasta",
        # --- eggNOG (full proteome) ---
        "eggnog_annotations": project_root
        / "eggnog_full_results"
        / "full_proteome_annotation.emapper.annotations",
        # --- InterProScan (up/down subsets) ---
        "interpro_up": project_root / "interproscan_results" / "up" / "up_interproscan.tsv",
        "interpro_down": project_root / "interproscan_results" / "down" / "down_interproscan.tsv",
        # --- OrthoFinder (full proteome) ---
        "orthofinder_results_dir": orthofinder_results,
        "orthogroups_tsv": orthofinder_results / "Orthogroups" / "Orthogroups.tsv",
        "orthogroups_unassigned_tsv": orthofinder_results
        / "Orthogroups"
        / "Orthogroups_UnassignedGenes.tsv",
        "orthologues_dir": orthofinder_results / "Orthologues",
        # --- SignalP (up/down subsets) ---
        "signalp_up": project_root / "signalp_up_results" / "up_prediction_results.txt",
        "signalp_down": project_root / "signalp_down_results" / "down_prediction_results.txt",
        # --- DeepTMHMM (up/down subsets) ---
        "deeptmhmm_up_gff3": project_root
        / "deeptmhmm_results" / "up" / "biolib_results" / "TMRs.gff3",
        "deeptmhmm_down_gff3": project_root
        / "deeptmhmm_results" / "down" / "biolib_results" / "TMRs.gff3",
        "deeptmhmm_up_3line": project_root
        / "deeptmhmm_results" / "up" / "biolib_results" / "predicted_topologies.3line",
        "deeptmhmm_down_3line": project_root
        / "deeptmhmm_results" / "down" / "biolib_results" / "predicted_topologies.3line",
        # --- PredGPI (up/down subsets) ---
        "predgpi_up": project_root / "predgpi_results" / "up_predgpi.json",
        "predgpi_down": project_root / "predgpi_results" / "down_predgpi.json",
        # --- output ---
        "output_dir": output_dir,
    }


# Species-name keywords used to *discover* OrthoFinder's files by pattern
# rather than hard-coding exact FungiDB build strings. If FungiDB is
# updated (e.g. FungiDB-69 instead of FungiDB-68), OrthoFinder's directory
# and file names will change, but these keywords will still match.
SPECIES_KEYWORDS = {
    "focal": "Cauris",
    "candida_albicans": "Calbicans",
    "saccharomyces_cerevisiae": "Scerevisiae",
    "schizosaccharomyces_pombe": "Spombe",
}

# Delimiter used for every multi-value field in the final table.
DELIM = ";"

# Tokens that mean "no value" across the various raw source formats.
MISSING_TOKENS = {"-", "", "nan", "NaN", "None", "N/A"}

GENE_ID_PATTERN = re.compile(r"B9J08_[0-9]+")

EGGNOG_COLUMNS = [
    "query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl",
    "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko",
    "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE",
    "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs",
]

INTERPRO_COLUMNS = [
    "protein_accession", "md5", "seq_length", "analysis", "signature_accession",
    "signature_description", "start", "stop", "score", "status", "date",
    "interpro_accession", "interpro_description", "go_terms", "pathways",
]

SIGNALP_COLUMNS = ["protein_accession", "prediction", "prob_other", "prob_sp", "cs_position"]


# =============================================================================
# LOGGING
# =============================================================================


def setup_logging(log_path: Path) -> logging.Logger:
    """
    Configure a logger that writes INFO-and-above to both the console and a
    QC log file. This log file IS the QC report deliverable
    (PCG_functional_annotation_QC.txt) -- every QC statement made while the
    script runs is captured here automatically, so there is no separate
    bookkeeping step that could drift out of sync with what actually ran.
    """
    logger = logging.getLogger("pcg_annotation")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    formatter = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%Y-%m-%d %H:%M:%S")

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    log_path.parent.mkdir(parents=True, exist_ok=True)
    file_handler = logging.FileHandler(log_path, mode="w")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return logger


def section(logger: logging.Logger, title: str) -> None:
    """Print a visually distinct section header in the log/QC report."""
    bar = "=" * 78
    logger.info(bar)
    logger.info(title)
    logger.info(bar)


# =============================================================================
# GENERIC HELPERS (error checking, ID normalisation, value collapsing)
# =============================================================================


def require_file(path: Path, description: str) -> None:
    """Raise a clear, actionable error if a required input file is missing."""
    if not path.is_file():
        raise FileNotFoundError(
            f"Required input ({description}) not found at: {path}\n"
            f"Check the CONFIG paths at the top of this script / your --project-root."
        )


def require_dir(path: Path, description: str) -> None:
    """Raise a clear, actionable error if a required input directory is missing."""
    if not path.is_dir():
        raise FileNotFoundError(
            f"Required input directory ({description}) not found at: {path}\n"
            f"Check the CONFIG paths at the top of this script / your --project-root."
        )


def clean_missing(value) -> Optional[str]:
    """
    Normalise a raw cell value to either a clean stripped string, or None if
    it represents "no value" in any of the conventions used across these
    six tools ('-' for eggNOG/InterPro, '' for SignalP CS position, NaN from
    pandas, etc.).
    """
    if value is None:
        return None
    if isinstance(value, float) and pd.isna(value):
        return None
    text = str(value).strip()
    if text in MISSING_TOKENS:
        return None
    return text


def normalize_gene_id(raw_id) -> Optional[str]:
    """Extract the canonical B9J08_###### gene ID from any raw identifier."""
    if raw_id is None:
        return None
    match = GENE_ID_PATTERN.search(str(raw_id))
    return match.group(0) if match else None


def normalize_id_series(
    raw_ids: Sequence[str], source_name: str, logger: logging.Logger
) -> List[str]:
    """
    Normalise a whole list of raw IDs to gene IDs, and HARD FAIL (rather
    than silently dropping rows) if any raw ID does not contain the
    expected B9J08_###### pattern. This is the "confirm no unexpected
    non-B9J08 IDs are being treated as C. auris genes" check, applied
    up front to every source before it touches the merge.
    """
    normalized = [normalize_gene_id(r) for r in raw_ids]
    bad = [raw for raw, norm in zip(raw_ids, normalized) if norm is None]
    if bad:
        examples = bad[:10]
        raise ValueError(
            f"{source_name}: {len(bad)} raw ID(s) did not match the expected "
            f"'B9J08_<digits>' pattern and cannot be treated as C. auris gene "
            f"IDs. Examples: {examples}"
        )
    logger.info(f"{source_name}: normalised {len(raw_ids)} raw ID(s), all matched B9J08_<digits>.")
    return normalized


def unique_join(values: Iterable, delimiter: str = DELIM) -> Optional[str]:
    """Collapse an iterable of scalar values to a sorted, unique, delimited string."""
    cleaned = {clean_missing(v) for v in values}
    cleaned.discard(None)
    if not cleaned:
        return None
    return delimiter.join(sorted(cleaned))


def split_and_collect(values: Iterable, split_on: str = "|", delimiter: str = DELIM) -> Optional[str]:
    """
    Like unique_join, but each input value may itself be a delimited list
    (e.g. InterProScan's GO/pathway cells are '|'-joined). Flattens first,
    then deduplicates. Used so multi-value cells from different domain-hit
    rows are merged token-by-token rather than kept as separate opaque
    strings.
    """
    tokens: List[str] = []
    for v in values:
        cv = clean_missing(v)
        if cv is None:
            continue
        tokens.extend(part.strip() for part in cv.split(split_on) if part.strip())
    return unique_join(tokens, delimiter=delimiter)


def report_duplicates(
    df: pd.DataFrame, id_col: str, source_name: str, logger: logging.Logger, expected: bool
) -> pd.Series:
    """
    QC requirement: "print duplicated PCG IDs before collapsing". Returns
    the duplicate-count-per-ID series for the caller to use downstream.

    expected=True  -> duplicates are normal (e.g. InterProScan has many
                       domain-hit rows per protein); just report the count.
    expected=False -> duplicates would be a red flag for a source that
                       should already be one-row-per-protein (eggNOG,
                       SignalP, DeepTMHMM, PredGPI); report loudly, but do
                       not raise here -- the caller decides whether to
                       collapse defensively or fail, since a handful of
                       genuine re-runs/appends could plausibly happen.
    """
    counts = df[id_col].value_counts()
    dup_counts = counts[counts > 1]
    if dup_counts.empty:
        logger.info(f"{source_name}: no duplicate {id_col} values before collapsing.")
    else:
        level = logging.INFO if expected else logging.WARNING
        logger.log(
            level,
            f"{source_name}: {len(dup_counts)} {id_col} value(s) have multiple rows "
            f"before collapsing (expected={expected}). "
            f"Max rows for one ID: {dup_counts.max()}. "
            f"Examples: {dup_counts.head(5).to_dict()}",
        )
    return dup_counts


def assert_one_row_per_gene(df: pd.DataFrame, id_col: str, source_name: str) -> None:
    """Hard-fail if a table that is supposed to be one-row-per-gene is not."""
    n_rows = len(df)
    n_unique = df[id_col].nunique()
    if n_rows != n_unique:
        dup_ids = df[df.duplicated(id_col, keep=False)][id_col].unique()
        raise ValueError(
            f"{source_name}: expected exactly one row per {id_col} after "
            f"collapsing, but got {n_rows} rows for {n_unique} unique IDs. "
            f"Duplicated IDs: {list(dup_ids)[:20]}"
        )


# =============================================================================
# BACKBONE: the full list of C. auris PCGs, from the proteome FASTA that
# was actually fed into eggNOG and OrthoFinder (NOT from eggNOG's output,
# which silently omits genes with no DIAMOND hit).
# =============================================================================


def parse_fasta_header_fields(header: str) -> Dict[str, str]:
    """
    Parse a FungiDB-style pipe-delimited FASTA header into a dict of its
    'key=value' fields, plus the leading accession token under key 'accession'.

    Example header:
        B9J08_000001-t37_1-p1 | transcript=B9J08_000001-t37_1 | gene=B9J08_000001 | ...
    """
    parts = [p.strip() for p in header.split("|")]
    fields: Dict[str, str] = {"accession": parts[0].split()[0]}
    for part in parts[1:]:
        if "=" in part:
            key, _, value = part.partition("=")
            fields[key.strip()] = value.strip()
    return fields


def build_gene_backbone(fasta_path: Path, logger: logging.Logger) -> pd.DataFrame:
    """
    Build the master backbone: one row per C. auris PCG, with gene_id,
    transcript_id and protein_id taken directly from the FASTA header's own
    structured fields (not guessed via suffix-stripping).

    Also asserts that every entry is a true, non-pseudo protein-coding gene
    (is_pseudo=false, SO=protein_coding_gene), matching what was confirmed
    by direct inspection of this file -- if a future re-run of this
    pipeline ever includes a pseudogene, this assertion will catch it
    rather than silently folding it into the PCG table.
    """
    require_file(fasta_path, "C. auris proteome FASTA (backbone)")

    records = []
    with open(fasta_path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            fields = parse_fasta_header_fields(line[1:].strip())
            records.append(fields)

    logger.info(f"Backbone FASTA: read {len(records)} header(s) from {fasta_path}")

    missing_fields = [r for r in records if "gene" not in r or "transcript" not in r]
    if missing_fields:
        raise ValueError(
            f"{len(missing_fields)} FASTA header(s) are missing an expected "
            f"'gene=' or 'transcript=' field. Example: {missing_fields[0]}"
        )

    non_pcg = [
        r
        for r in records
        if r.get("is_pseudo", "false").lower() != "false"
        or r.get("SO", "protein_coding_gene") != "protein_coding_gene"
    ]
    if non_pcg:
        logger.warning(
            f"Backbone FASTA contains {len(non_pcg)} entr(y/ies) NOT flagged as a "
            f"non-pseudo protein_coding_gene. These are being KEPT (not silently "
            f"dropped) -- review manually. Examples: {non_pcg[:3]}"
        )
    else:
        logger.info(
            "Backbone FASTA: confirmed all entries are is_pseudo=false, "
            "SO=protein_coding_gene."
        )

    protein_ids = [r["accession"] for r in records]
    transcript_ids = [r["transcript"] for r in records]
    raw_gene_ids = [r["gene"] for r in records]

    gene_ids = normalize_id_series(raw_gene_ids, "backbone FASTA gene= field", logger)

    df = pd.DataFrame(
        {"gene_id": gene_ids, "transcript_id": transcript_ids, "protein_id": protein_ids}
    )

    assert_one_row_per_gene(df, "gene_id", "backbone FASTA")
    logger.info(f"Backbone: {len(df)} unique C. auris PCGs established as the master gene set.")
    return df


# =============================================================================
# SOURCE 1: eggNOG-mapper (full proteome)
# =============================================================================


def parse_eggnog(path: Path, logger: logging.Logger) -> pd.DataFrame:
    """
    Parse eggNOG-mapper's *.emapper.annotations file.

    Format: emapper prefixes metadata/footer lines with '#' (both the
    '##...' run-info lines and the single '#query...' header line); all
    such lines are skipped, and the fixed 21-column schema (confirmed by
    direct inspection) is supplied explicitly rather than trusted from the
    file's own header text.

    Genes with no DIAMOND hit at all are simply ABSENT from this file
    (not present as an all-'-' row) -- this is why the master backbone
    is NOT built from eggNOG output.
    """
    require_file(path, "eggNOG annotations")

    data_lines = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            line = line.rstrip("\n")
            if not line:
                continue
            data_lines.append(line.split("\t"))

    logger.info(f"eggNOG: {len(data_lines)} raw data row(s) read from {path}")

    bad_rows = [row for row in data_lines if len(row) != len(EGGNOG_COLUMNS)]
    if bad_rows:
        raise ValueError(
            f"eggNOG: {len(bad_rows)} row(s) do not have the expected "
            f"{len(EGGNOG_COLUMNS)} columns. First bad row: {bad_rows[0]}"
        )

    df = pd.DataFrame(data_lines, columns=EGGNOG_COLUMNS)

    report_duplicates(df, "query", "eggNOG (raw, pre-normalisation)", logger, expected=False)

    df["gene_id"] = normalize_id_series(df["query"].tolist(), "eggNOG", logger)
    assert_one_row_per_gene(df, "gene_id", "eggNOG")

    scalar_cols = {
        "Description": "eggNOG_Description",
        "Preferred_name": "eggNOG_Preferred_name",
        "seed_ortholog": "eggNOG_seed_ortholog",
        "COG_category": "eggNOG_COG_category",
    }
    multi_cols = {
        "GOs": "eggNOG_GOs",
        "KEGG_ko": "eggNOG_KEGG_ko",
        "KEGG_Pathway": "eggNOG_KEGG_Pathway",
        "PFAMs": "eggNOG_PFAMs",
    }

    out = pd.DataFrame({"gene_id": df["gene_id"]})
    for src, dst in scalar_cols.items():
        out[dst] = df[src].apply(clean_missing)
    for src, dst in multi_cols.items():
        # eggNOG already comma-delimits multi-value fields; re-delimit to
        # ';' for consistency with every other source in the final table.
        out[dst] = df[src].apply(
            lambda v: clean_missing(v).replace(",", DELIM) if clean_missing(v) else None
        )

    n_annotated = out["eggNOG_Description"].notna().sum()
    logger.info(
        f"eggNOG: {len(out)} genes parsed; {n_annotated} have a non-empty Description "
        f"({len(out) - n_annotated} scanned but no functional hit)."
    )
    return out


# =============================================================================
# SOURCE 2: InterProScan (up/down significant-gene subsets)
# =============================================================================


def _read_interpro_tsv(path: Path, logger: logging.Logger) -> pd.DataFrame:
    """
    Read one raw InterProScan TSV (headerless, 15 fixed columns) into a
    DataFrame of strings, and confirm no completely empty row slipped in.
    """
    require_file(path, f"InterProScan TSV ({path.name})")
    # quoting=csv.QUOTE_NONE: InterProScan free-text descriptions can contain
    # literal double-quote characters (e.g. "Winged helix" DNA-binding
    # domain) that are NOT CSV quoting -- without this, pandas misparses
    # the file looking for a matching closing quote.
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=INTERPRO_COLUMNS,
        dtype=str,
        engine="python",
        quoting=csv.QUOTE_NONE,
    )
    bad_rows = df[df.isna().all(axis=1)]
    if len(bad_rows):
        raise ValueError(f"{path}: {len(bad_rows)} completely empty row(s) found.")
    logger.info(f"InterProScan: {len(df)} raw row(s) read from {path}")
    return df


def parse_interpro(
    up_path: Path, down_path: Path, logger: logging.Logger
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Parse and collapse InterProScan's up+down TSV outputs to one row per
    gene. Returns (collapsed_df, raw_concat_df) -- the raw frame is
    returned too so the caller can compute InterPro_tested from the exact
    set of genes actually submitted to InterProScan.

    Multiple domain/database hits per protein are NOT discarded: every
    field is collapsed to the sorted-unique set of values seen across all
    of that gene's hit rows, joined with ';'.
    """
    up_raw = _read_interpro_tsv(up_path, logger)
    down_raw = _read_interpro_tsv(down_path, logger)

    # A gene can legitimately appear in both the up and down subset files
    # (confirmed by inspection: 3 genes do here). Concatenate BEFORE
    # grouping so such genes get one unified collapsed row, not two
    # partial ones that would silently overwrite each other on merge.
    raw = pd.concat([up_raw, down_raw], ignore_index=True)
    overlap = set(up_raw["protein_accession"]) & set(down_raw["protein_accession"])
    if overlap:
        logger.warning(
            f"InterProScan: {len(overlap)} raw protein ID(s) appear in BOTH the up and "
            f"down files (e.g. {sorted(overlap)[:5]}); their hit rows have been "
            f"concatenated before collapsing so no information is lost."
        )

    raw["gene_id"] = normalize_id_series(
        raw["protein_accession"].tolist(), "InterProScan", logger
    )

    dup_counts = report_duplicates(raw, "gene_id", "InterProScan (raw)", logger, expected=True)
    logger.info(
        f"InterProScan: {raw['gene_id'].nunique()} unique gene(s) tested, "
        f"{len(raw)} total domain-hit rows "
        f"(mean {len(raw) / raw['gene_id'].nunique():.1f} hits/gene)."
    )

    def collapse_group(g: pd.DataFrame) -> pd.Series:
        return pd.Series(
            {
                "InterPro_accessions": unique_join(g["interpro_accession"]),
                "InterPro_descriptions": unique_join(g["interpro_description"]),
                "InterPro_databases": unique_join(g["analysis"]),
                "InterPro_signature_accessions": unique_join(g["signature_accession"]),
                "InterPro_signature_descriptions": unique_join(g["signature_description"]),
                "InterPro_GOs": split_and_collect(g["go_terms"], split_on="|"),
                "InterPro_pathways": split_and_collect(g["pathways"], split_on="|"),
            }
        )

    collapsed = raw.groupby("gene_id", sort=False).apply(collapse_group).reset_index()
    assert_one_row_per_gene(collapsed, "gene_id", "InterProScan (collapsed)")
    logger.info(f"InterProScan: collapsed to {len(collapsed)} row(s), one per gene.")

    return collapsed, raw


# =============================================================================
# SOURCE 3: OrthoFinder (full proteome)
# =============================================================================


def _find_orthofinder_species_column(header_cols: Sequence[str], keyword: str) -> str:
    """
    Locate the single OrthoFinder column whose header contains a species
    keyword (see SPECIES_KEYWORDS), rather than assuming a fixed column
    position -- OrthoFinder's own column order differs between
    Orthogroups.tsv (alphabetical) and the pairwise Orthologues/ files
    (focal-species-first), so position-based lookup would be wrong for one
    or the other. Raises if the match is not exactly one column.
    """
    matches = [c for c in header_cols if keyword.lower() in c.lower()]
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one OrthoFinder column matching species keyword "
            f"'{keyword}', found {len(matches)}: {matches}. Columns present: {header_cols}"
        )
    return matches[0]


def _explode_comma_column(
    df: pd.DataFrame, id_col: str, value_col: Optional[str] = None
) -> pd.DataFrame:
    """
    Explode a comma-separated gene-list column so each individual raw ID
    gets its own row, optionally carrying along one other column's value
    (e.g. the row's Orthogroup, or the paired other-species gene list) --
    used for both Orthogroups.tsv-style and Orthologues-style tables.
    """
    records = []
    for _, row in df.iterrows():
        cell = clean_missing(row[id_col])
        if cell is None:
            continue
        ids = [x.strip() for x in cell.split(",") if x.strip()]
        for raw_id in ids:
            rec = {"raw_id": raw_id}
            if value_col is not None:
                rec["value"] = row[value_col]
            records.append(rec)
    return pd.DataFrame(records)


def parse_orthofinder_orthogroups(
    orthogroups_path: Path, unassigned_path: Path, logger: logging.Logger
) -> pd.DataFrame:
    """
    Parse Orthogroups.tsv + Orthogroups_UnassignedGenes.tsv to assign every
    C. auris gene its Orthogroup ID. Together these two files give FULL
    proteome coverage (confirmed by inspection: 5093 + 324 = 5417), unlike
    the pairwise Orthologues/ files used below for actual ortholog calls.

    Column position for the C. auris species is NOT hard-coded -- these two
    files order species alphabetically, which differs from the pairwise
    Orthologues/ files (focal-species-first), so the column is located by
    name every time.
    """
    require_file(orthogroups_path, "OrthoFinder Orthogroups.tsv")
    require_file(unassigned_path, "OrthoFinder Orthogroups_UnassignedGenes.tsv")

    frames = []
    for path in (orthogroups_path, unassigned_path):
        df = pd.read_csv(path, sep="\t", dtype=str, quoting=csv.QUOTE_NONE)
        cauris_col = _find_orthofinder_species_column(df.columns, SPECIES_KEYWORDS["focal"])
        exploded = _explode_comma_column(df, cauris_col, value_col="Orthogroup")
        exploded = exploded.rename(columns={"value": "Orthogroup"})
        frames.append(exploded)
        logger.info(f"OrthoFinder: {path.name} -> {len(exploded)} gene-OG assignment row(s).")

    combined = pd.concat(frames, ignore_index=True)
    combined["gene_id"] = normalize_id_series(
        combined["raw_id"].tolist(), "OrthoFinder Orthogroup assignment", logger
    )

    report_duplicates(combined, "gene_id", "OrthoFinder Orthogroup assignment", logger, expected=False)
    assert_one_row_per_gene(combined, "gene_id", "OrthoFinder Orthogroup assignment")

    out = combined[["gene_id", "Orthogroup"]].copy()
    logger.info(f"OrthoFinder: {len(out)} genes assigned an Orthogroup ID.")
    return out


def _find_orthologues_file(orthologues_dir: Path, focal_keyword: str, target_keyword: str) -> Path:
    """
    Discover the pairwise Orthologues/*.tsv file for (focal, target)
    species by matching directory/file names against species-name
    keywords, rather than hard-coding the exact FungiDB build filenames
    (which change if FungiDB is updated, e.g. FungiDB-68 -> FungiDB-69).
    Raises if the focal subdirectory or the target file is not found
    uniquely.
    """
    focal_dirs = [
        d
        for d in orthologues_dir.iterdir()
        if d.is_dir() and focal_keyword.lower() in d.name.lower()
    ]
    if len(focal_dirs) != 1:
        raise FileNotFoundError(
            f"Expected exactly one Orthologues_* subdirectory matching "
            f"'{focal_keyword}' under {orthologues_dir}, found {len(focal_dirs)}."
        )
    focal_dir = focal_dirs[0]

    candidate_files = [
        f
        for f in focal_dir.glob("*.tsv")
        if focal_keyword.lower() in f.name.lower() and target_keyword.lower() in f.name.lower()
    ]
    if len(candidate_files) != 1:
        raise FileNotFoundError(
            f"Expected exactly one pairwise orthologues TSV matching focal="
            f"'{focal_keyword}' and target='{target_keyword}' under {focal_dir}, "
            f"found {len(candidate_files)}: {candidate_files}"
        )
    return candidate_files[0]


def parse_orthofinder_pairwise_orthologs(
    orthologues_dir: Path, logger: logging.Logger
) -> pd.DataFrame:
    """
    Parse the three focal-species pairwise Orthologues/*.tsv files
    (C. auris vs C. albicans / S. cerevisiae / S. pombe) to assign each
    C. auris gene its inferred ortholog(s) in each of the other three
    species. Coverage here is intentionally PARTIAL relative to the
    Orthogroup assignment above -- OrthoFinder ran on the full proteome, so
    a gene with no entry here has no confidently inferred tree-based
    orthologue relationship to that species. This is a genuine biological
    result, not a "not tested" gap (OrthoFinder has no subset-testing
    concept the way InterProScan/SignalP/DeepTMHMM/PredGPI do).

    Each C. auris gene appears in at most one row per pairwise file
    (confirmed by inspection), so no cross-row aggregation is required --
    but the code still groups defensively rather than assuming this holds
    forever.
    """
    require_dir(orthologues_dir, "OrthoFinder Orthologues directory")

    target_species = {
        SPECIES_KEYWORDS["candida_albicans"]: "Candida_albicans_orthologs",
        SPECIES_KEYWORDS["saccharomyces_cerevisiae"]: "Saccharomyces_cerevisiae_orthologs",
        SPECIES_KEYWORDS["schizosaccharomyces_pombe"]: "Schizosaccharomyces_pombe_orthologs",
    }

    per_species_frames = []
    for target_keyword, out_col in target_species.items():
        path = _find_orthologues_file(orthologues_dir, SPECIES_KEYWORDS["focal"], target_keyword)
        df = pd.read_csv(path, sep="\t", dtype=str, quoting=csv.QUOTE_NONE)

        focal_col = _find_orthofinder_species_column(df.columns, SPECIES_KEYWORDS["focal"])
        target_col = _find_orthofinder_species_column(df.columns, target_keyword)

        exploded = _explode_comma_column(df, focal_col, value_col=target_col)
        exploded = exploded.rename(columns={"value": "target_raw"})
        exploded["gene_id"] = normalize_id_series(
            exploded["raw_id"].tolist(), f"OrthoFinder orthologues ({out_col})", logger
        )

        report_duplicates(
            exploded, "gene_id", f"OrthoFinder orthologues ({out_col}, raw)", logger, expected=False
        )

        collapsed = (
            exploded.groupby("gene_id", sort=False)["target_raw"]
            .apply(unique_join)
            .rename(out_col)
            .reset_index()
        )
        assert_one_row_per_gene(collapsed, "gene_id", f"OrthoFinder orthologues ({out_col})")
        logger.info(
            f"OrthoFinder: {path.name} -> {len(collapsed)} C. auris gene(s) with an "
            f"inferred {out_col.replace('_', ' ')}."
        )
        per_species_frames.append(collapsed)

    out = per_species_frames[0]
    for frame in per_species_frames[1:]:
        out = out.merge(frame, on="gene_id", how="outer")
    return out


def parse_orthofinder(config: Dict[str, Path], logger: logging.Logger) -> pd.DataFrame:
    """Combine the Orthogroup assignment and the three pairwise ortholog tables."""
    orthogroups = parse_orthofinder_orthogroups(
        config["orthogroups_tsv"], config["orthogroups_unassigned_tsv"], logger
    )
    orthologs = parse_orthofinder_pairwise_orthologs(config["orthologues_dir"], logger)
    combined = orthogroups.merge(orthologs, on="gene_id", how="outer")
    assert_one_row_per_gene(combined, "gene_id", "OrthoFinder (combined)")
    logger.info(f"OrthoFinder: combined table has {len(combined)} row(s).")
    return combined


# =============================================================================
# SOURCE 4: SignalP 6.0 (up/down significant-gene subsets)
# =============================================================================


def _read_signalp_file(path: Path, logger: logging.Logger) -> pd.DataFrame:
    """Read one raw SignalP 6.0 prediction file (2 leading '#' header lines, then 5 fixed columns)."""
    require_file(path, f"SignalP predictions ({path.name})")
    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        names=SIGNALP_COLUMNS,
        dtype=str,
        engine="python",
        quoting=csv.QUOTE_NONE,
    )
    logger.info(f"SignalP: {len(df)} raw row(s) read from {path}")
    return df


def parse_signalp(up_path: Path, down_path: Path, logger: logging.Logger) -> pd.DataFrame:
    """
    Parse SignalP 6.0's up/down prediction files. One row per protein
    already (confirmed by inspection); still checked defensively.
    """
    up = _read_signalp_file(up_path, logger)
    down = _read_signalp_file(down_path, logger)
    raw = pd.concat([up, down], ignore_index=True)

    overlap = set(up["protein_accession"]) & set(down["protein_accession"])
    if overlap:
        logger.warning(
            f"SignalP: {len(overlap)} raw protein ID(s) appear in BOTH up and down "
            f"files (e.g. {sorted(overlap)[:5]})."
        )

    raw["gene_id"] = normalize_id_series(raw["protein_accession"].tolist(), "SignalP", logger)
    report_duplicates(raw, "gene_id", "SignalP (raw)", logger, expected=False)

    def collapse_group(g: pd.DataFrame) -> pd.Series:
        return pd.Series(
            {
                "SignalP_prediction": unique_join(g["prediction"]),
                "SignalP_SP_probability": unique_join(g["prob_sp"]),
                "SignalP_cleavage_site": unique_join(g["cs_position"]),
            }
        )

    collapsed = raw.groupby("gene_id", sort=False).apply(collapse_group).reset_index()
    assert_one_row_per_gene(collapsed, "gene_id", "SignalP (collapsed)")
    collapsed["SignalP_tested"] = True

    n_sp = (collapsed["SignalP_prediction"] == "SP").sum()
    logger.info(
        f"SignalP: {len(collapsed)} genes tested, {n_sp} predicted SP "
        f"(signal peptide), {len(collapsed) - n_sp} predicted OTHER."
    )
    return collapsed


# =============================================================================
# SOURCE 5: DeepTMHMM (up/down significant-gene subsets)
# =============================================================================


_TMHMM_LENGTH_RE = re.compile(r"^#\s*(\S+)\s+Length:\s*(\d+)\s*$")
_TMHMM_COUNT_RE = re.compile(r"^#\s*(\S+)\s+Number of predicted TMRs:\s*(\d+)\s*$")


def _parse_tmrs_gff3(path: Path, logger: logging.Logger) -> pd.DataFrame:
    """
    Parse DeepTMHMM's TMRs.gff3 (a DeepTMHMM-specific pseudo-GFF3, not
    standard GFF3): per protein, a '# ID Length: N' comment, a
    '# ID Number of predicted TMRs: N' comment (used directly as the TM
    count -- no need to count segments manually), then region rows typed
    inside/outside/TMhelix/signal, terminated by '//'.
    """
    require_file(path, f"DeepTMHMM TMRs.gff3 ({path.name})")

    records: Dict[str, Dict] = {}
    current_id: Optional[str] = None

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("##"):
                continue
            if line == "//":
                current_id = None
                continue

            m_len = _TMHMM_LENGTH_RE.match(line)
            if m_len:
                current_id, length = m_len.groups()
                records[current_id] = {"length": int(length), "tm_count": None, "regions": []}
                continue

            m_count = _TMHMM_COUNT_RE.match(line)
            if m_count:
                prot_id, count = m_count.groups()
                records.setdefault(prot_id, {"length": None, "tm_count": None, "regions": []})
                records[prot_id]["tm_count"] = int(count)
                continue

            if line.startswith("#"):
                continue  # any other comment line, ignored

            fields = line.split("\t")
            if len(fields) < 4:
                continue
            prot_id, region_type, start, stop = fields[0], fields[1], fields[2], fields[3]
            records.setdefault(prot_id, {"length": None, "tm_count": None, "regions": []})
            # Membrane-spanning segments are labelled 'TMhelix' for the usual
            # alpha-helical topology class, but 'Beta sheet' for DeepTMHMM's
            # BETA (beta-barrel / outer-membrane-type) class -- confirmed by
            # inspecting a real BETA-class protein whose 'Number of
            # predicted TMRs' count only matched when both were counted.
            # Both count towards 'Number of predicted TMRs', so both are
            # captured here (tagged with their region type so helix vs.
            # strand remains distinguishable in the output).
            if region_type in ("TMhelix", "Beta sheet"):
                records[prot_id]["regions"].append(f"{region_type}:{start}-{stop}")

    missing_count = [k for k, v in records.items() if v["tm_count"] is None]
    if missing_count:
        raise ValueError(
            f"{path}: {len(missing_count)} protein(s) never had a "
            f"'Number of predicted TMRs' line parsed. Example: {missing_count[0]}"
        )

    logger.info(f"DeepTMHMM: {len(records)} protein block(s) parsed from {path.name}")

    df = pd.DataFrame(
        [
            {
                "protein_accession": pid,
                "DeepTMHMM_TM_count": v["tm_count"],
                "DeepTMHMM_TM_regions": DELIM.join(v["regions"]) if v["regions"] else None,
            }
            for pid, v in records.items()
        ]
    )

    # Cross-check: the number of parsed TMhelix regions should equal the
    # protein's own reported TM count, for every protein.
    mismatch = df[
        df["DeepTMHMM_TM_regions"].apply(lambda s: 0 if pd.isna(s) else len(s.split(DELIM)))
        != df["DeepTMHMM_TM_count"]
    ]
    if len(mismatch):
        raise ValueError(
            f"{path}: {len(mismatch)} protein(s) have a 'Number of predicted TMRs' "
            f"count that does not match the number of TMhelix regions actually "
            f"parsed. Example: {mismatch.iloc[0].to_dict()}"
        )

    return df


def _parse_3line_classes(path: Path, logger: logging.Logger) -> pd.DataFrame:
    """
    Parse predicted_topologies.3line's FASTA-like headers
    ('>ID | CLASS') to extract DeepTMHMM's own overall topology
    classification (GLOB / TM / SP / SP+TM / BETA).
    """
    require_file(path, f"DeepTMHMM predicted_topologies.3line ({path.name})")

    records = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            if "|" not in header:
                raise ValueError(f"{path}: malformed 3line header (no '|'): {header!r}")
            prot_id, _, cls = header.partition("|")
            records.append({"protein_accession": prot_id.strip(), "DeepTMHMM_topology_class": cls.strip()})

    logger.info(f"DeepTMHMM: {len(records)} topology class label(s) parsed from {path.name}")
    return pd.DataFrame(records)


def parse_deeptmhmm(
    up_gff3: Path, down_gff3: Path, up_3line: Path, down_3line: Path, logger: logging.Logger
) -> pd.DataFrame:
    """Combine the TMRs.gff3 (counts/regions) and .3line (class label) parses."""
    gff3 = pd.concat(
        [_parse_tmrs_gff3(up_gff3, logger), _parse_tmrs_gff3(down_gff3, logger)],
        ignore_index=True,
    )
    threeline = pd.concat(
        [_parse_3line_classes(up_3line, logger), _parse_3line_classes(down_3line, logger)],
        ignore_index=True,
    )

    overlap = set(gff3["protein_accession"][gff3["protein_accession"].duplicated()])
    if overlap:
        logger.warning(
            f"DeepTMHMM: {len(overlap)} raw protein ID(s) appear in BOTH up and down "
            f"TMRs.gff3 files (e.g. {sorted(overlap)[:5]})."
        )

    gff3_genes = set(normalize_gene_id(x) for x in gff3["protein_accession"])
    threeline_genes = set(normalize_gene_id(x) for x in threeline["protein_accession"])
    if gff3_genes != threeline_genes:
        only_gff3 = gff3_genes - threeline_genes
        only_3line = threeline_genes - gff3_genes
        logger.warning(
            f"DeepTMHMM: TMRs.gff3 and predicted_topologies.3line do not cover the "
            f"exact same gene set ({len(only_gff3)} only in gff3, {len(only_3line)} "
            f"only in 3line). Proceeding with an outer merge; missing side will be NaN."
        )

    merged = gff3.merge(threeline, on="protein_accession", how="outer")
    merged["gene_id"] = normalize_id_series(
        merged["protein_accession"].tolist(), "DeepTMHMM", logger
    )

    report_duplicates(merged, "gene_id", "DeepTMHMM (raw)", logger, expected=False)

    def collapse_group(g: pd.DataFrame) -> pd.Series:
        return pd.Series(
            {
                "DeepTMHMM_TM_count": unique_join(g["DeepTMHMM_TM_count"].astype(str)),
                "DeepTMHMM_TM_regions": unique_join(g["DeepTMHMM_TM_regions"]),
                "DeepTMHMM_topology_class": unique_join(g["DeepTMHMM_topology_class"]),
            }
        )

    collapsed = merged.groupby("gene_id", sort=False).apply(collapse_group).reset_index()
    assert_one_row_per_gene(collapsed, "gene_id", "DeepTMHMM (collapsed)")

    # TM_count is numeric in the vast majority of cases (one raw ID per
    # gene); coerce back to a clean integer where unambiguous.
    def to_int_or_keep(v):
        if pd.isna(v):
            return None
        parts = str(v).split(DELIM)
        return int(parts[0]) if len(parts) == 1 else v

    collapsed["DeepTMHMM_TM_count"] = collapsed["DeepTMHMM_TM_count"].apply(to_int_or_keep)
    collapsed["DeepTMHMM_has_TM"] = collapsed["DeepTMHMM_TM_count"].apply(
        lambda v: None if pd.isna(v) else (isinstance(v, int) and v > 0)
    )
    collapsed["DeepTMHMM_tested"] = True

    n_tm = (collapsed["DeepTMHMM_has_TM"] == True).sum()  # noqa: E712
    logger.info(
        f"DeepTMHMM: {len(collapsed)} genes tested, {n_tm} have >=1 predicted TM region."
    )
    return collapsed


# =============================================================================
# SOURCE 6: PredGPI (up/down significant-gene subsets)
# =============================================================================


def _parse_predgpi_json(path: Path, logger: logging.Logger) -> pd.DataFrame:
    """
    Parse one PredGPI JSON file (a list of {accession, features, sequence}
    entries) into one row per raw protein ID. A non-empty 'features' list
    is a positive GPI-anchor prediction (first feature's begin/score are
    kept); an empty list is a genuine tested-and-negative result, not a
    missing value.
    """
    require_file(path, f"PredGPI JSON ({path.name})")
    with open(path) as fh:
        data = json.load(fh)

    if not isinstance(data, list):
        raise ValueError(f"{path}: expected a JSON list at the top level, got {type(data)}.")

    logger.info(f"PredGPI: {len(data)} entr(y/ies) read from {path}")

    records = []
    multi_feature_warned = 0
    for entry in data:
        accession = entry.get("accession")
        if accession is None:
            raise ValueError(f"{path}: an entry is missing its 'accession' field: {entry}")

        features = entry.get("features", [])
        if len(features) > 1:
            multi_feature_warned += 1
        if features:
            feat = features[0]
            records.append(
                {
                    "protein_accession": accession,
                    "PredGPI_prediction": "GPI-anchored",
                    "PredGPI_anchor_position": feat.get("begin"),
                    "PredGPI_score": feat.get("score"),
                }
            )
        else:
            records.append(
                {
                    "protein_accession": accession,
                    "PredGPI_prediction": "Not GPI-anchored",
                    "PredGPI_anchor_position": None,
                    "PredGPI_score": None,
                }
            )

    if multi_feature_warned:
        logger.warning(
            f"PredGPI: {multi_feature_warned} entr(y/ies) in {path.name} had MORE THAN "
            f"ONE feature; only the first was used. This did not occur in the data "
            f"originally inspected, so review these manually."
        )

    return pd.DataFrame(records)


def parse_predgpi(up_path: Path, down_path: Path, logger: logging.Logger) -> pd.DataFrame:
    """
    Parse PredGPI's up/down JSON outputs.

    'GPI-anchored' vs 'Not GPI-anchored' is a real tested-and-negative
    result whenever the accession is present in the JSON with an empty
    features list -- this is NOT the same as "not tested" (PredGPI_tested
    is set independently, from JSON membership, not from the prediction
    itself).
    """
    up = _parse_predgpi_json(up_path, logger)
    down = _parse_predgpi_json(down_path, logger)
    raw = pd.concat([up, down], ignore_index=True)

    overlap = set(up["protein_accession"]) & set(down["protein_accession"])
    if overlap:
        logger.warning(
            f"PredGPI: {len(overlap)} raw protein ID(s) appear in BOTH up and down "
            f"JSON files (e.g. {sorted(overlap)[:5]})."
        )

    raw["gene_id"] = normalize_id_series(raw["protein_accession"].tolist(), "PredGPI", logger)
    report_duplicates(raw, "gene_id", "PredGPI (raw)", logger, expected=False)

    def collapse_group(g: pd.DataFrame) -> pd.Series:
        return pd.Series(
            {
                "PredGPI_prediction": unique_join(g["PredGPI_prediction"]),
                "PredGPI_anchor_position": unique_join(g["PredGPI_anchor_position"].astype(str)),
                "PredGPI_score": unique_join(g["PredGPI_score"].astype(str)),
            }
        )

    collapsed = raw.groupby("gene_id", sort=False).apply(collapse_group).reset_index()
    assert_one_row_per_gene(collapsed, "gene_id", "PredGPI (collapsed)")
    collapsed["PredGPI_tested"] = True

    n_pos = (collapsed["PredGPI_prediction"] == "GPI-anchored").sum()
    logger.info(f"PredGPI: {len(collapsed)} genes tested, {n_pos} predicted GPI-anchored.")
    return collapsed


# =============================================================================
# MERGE ORCHESTRATION
# =============================================================================


def qc_before_merge(master: pd.DataFrame, step_name: str, logger: logging.Logger) -> int:
    """Log and return the master table's row count immediately before a merge step."""
    n = len(master)
    logger.info(f"[{step_name}] master row count before merge: {n}")
    return n


def qc_after_merge(
    master: pd.DataFrame, n_before: int, step_name: str, logger: logging.Logger
) -> None:
    """
    Confirm a merge step did not change the master table's row count and
    did not introduce a duplicated gene_id; raises immediately if either
    check fails, since a left-merge onto a one-row-per-gene backbone must
    never do either.
    """
    n_after = len(master)
    logger.info(f"[{step_name}] master row count after merge: {n_after}")
    if n_after != n_before:
        raise RuntimeError(
            f"[{step_name}] merge changed the master row count from {n_before} to "
            f"{n_after} -- a merge must NEVER change row count when merging onto a "
            f"backbone that already has one row per gene. Check the source table for "
            f"duplicate gene_id values."
        )
    n_dup = master["gene_id"].duplicated().sum()
    if n_dup:
        raise RuntimeError(
            f"[{step_name}] master table has {n_dup} duplicated gene_id value(s) "
            f"after merge."
        )


def clean_integer_display(value) -> Optional[str]:
    """
    Format a value that is semantically a single integer (a TM count, a
    residue position) for clean display. Left-joining a column that mixes
    real integers with missing values forces pandas to upcast the whole
    column to float64 (missing values have no integer NaN), which then
    renders every value with a spurious trailing '.0' (e.g. '0.0' instead
    of '0', '202.0' instead of '202') in the saved TSV/CSV. This corrects
    that cosmetic artifact WITHOUT touching genuinely multi-valued strings
    (e.g. a ';'-joined pair of differing values collapsed from the rare
    genes that appear in both the up and down subset files).
    """
    if pd.isna(value):
        return None
    text = str(value)
    if DELIM in text:
        return text  # genuinely multi-valued after collapsing; leave as-is
    try:
        return str(int(float(text)))
    except (TypeError, ValueError):
        return text


def report_source_coverage(
    master: pd.DataFrame, columns: Sequence[str], source_name: str, logger: logging.Logger
) -> None:
    """QC requirement: report how many / how few PCGs received annotation from each source."""
    any_annotated = master[list(columns)].notna().any(axis=1)
    n_annotated = int(any_annotated.sum())
    n_missing = len(master) - n_annotated
    logger.info(
        f"{source_name}: {n_annotated}/{len(master)} PCGs have at least one non-null "
        f"value among {list(columns)}; {n_missing} do not."
    )


def build_master_table(config: Dict[str, Path], logger: logging.Logger) -> Dict[str, pd.DataFrame]:
    """
    Run the full incremental parse -> QC -> merge -> QC pipeline. Returns a
    dict of {name: DataFrame} containing every collapsed intermediate table
    plus the final master table, so the caller can save all of them.
    """
    results: Dict[str, pd.DataFrame] = {}

    # ---- Backbone -----------------------------------------------------
    section(logger, "BACKBONE: C. auris proteome FASTA")
    master = build_gene_backbone(config["fasta_backbone"], logger)
    n_expected = len(master)

    # ---- 1: eggNOG ------------------------------------------------------
    section(logger, "SOURCE 1/6: eggNOG-mapper")
    eggnog = parse_eggnog(config["eggnog_annotations"], logger)
    results["eggnog_collapsed"] = eggnog
    n_before = qc_before_merge(master, "merge eggNOG", logger)
    master = master.merge(eggnog, on="gene_id", how="left")
    qc_after_merge(master, n_before, "merge eggNOG", logger)
    report_source_coverage(master, ["eggNOG_Description"], "eggNOG", logger)

    # ---- 2: InterProScan --------------------------------------------------
    section(logger, "SOURCE 2/6: InterProScan")
    interpro_collapsed, interpro_raw = parse_interpro(
        config["interpro_up"], config["interpro_down"], logger
    )
    interpro_tested_genes = set(interpro_raw["gene_id"])
    interpro_collapsed["InterPro_tested"] = True
    results["interpro_collapsed"] = interpro_collapsed
    n_before = qc_before_merge(master, "merge InterProScan", logger)
    master = master.merge(interpro_collapsed, on="gene_id", how="left")
    master["InterPro_tested"] = master["gene_id"].isin(interpro_tested_genes)
    qc_after_merge(master, n_before, "merge InterProScan", logger)
    report_source_coverage(master, ["InterPro_accessions"], "InterProScan", logger)

    # ---- 3: OrthoFinder -----------------------------------------------
    section(logger, "SOURCE 3/6: OrthoFinder")
    orthofinder = parse_orthofinder(config, logger)
    results["orthofinder_collapsed"] = orthofinder
    n_before = qc_before_merge(master, "merge OrthoFinder", logger)
    master = master.merge(orthofinder, on="gene_id", how="left")
    qc_after_merge(master, n_before, "merge OrthoFinder", logger)
    report_source_coverage(master, ["Orthogroup"], "OrthoFinder (Orthogroup)", logger)
    report_source_coverage(
        master, ["Candida_albicans_orthologs"], "OrthoFinder (C. albicans orthologs)", logger
    )
    report_source_coverage(
        master,
        ["Saccharomyces_cerevisiae_orthologs"],
        "OrthoFinder (S. cerevisiae orthologs)",
        logger,
    )
    report_source_coverage(
        master,
        ["Schizosaccharomyces_pombe_orthologs"],
        "OrthoFinder (S. pombe orthologs)",
        logger,
    )

    # ---- 4: SignalP -----------------------------------------------------
    section(logger, "SOURCE 4/6: SignalP 6.0")
    signalp = parse_signalp(config["signalp_up"], config["signalp_down"], logger)
    results["signalp_collapsed"] = signalp
    signalp_tested_genes = set(signalp["gene_id"])
    n_before = qc_before_merge(master, "merge SignalP", logger)
    master = master.merge(signalp, on="gene_id", how="left")
    master["SignalP_tested"] = master["gene_id"].isin(signalp_tested_genes)
    qc_after_merge(master, n_before, "merge SignalP", logger)
    report_source_coverage(master, ["SignalP_prediction"], "SignalP", logger)

    # ---- 5: DeepTMHMM -----------------------------------------------------
    section(logger, "SOURCE 5/6: DeepTMHMM")
    deeptmhmm = parse_deeptmhmm(
        config["deeptmhmm_up_gff3"],
        config["deeptmhmm_down_gff3"],
        config["deeptmhmm_up_3line"],
        config["deeptmhmm_down_3line"],
        logger,
    )
    results["deeptmhmm_collapsed"] = deeptmhmm
    deeptmhmm_tested_genes = set(deeptmhmm["gene_id"])
    n_before = qc_before_merge(master, "merge DeepTMHMM", logger)
    master = master.merge(deeptmhmm, on="gene_id", how="left")
    master["DeepTMHMM_tested"] = master["gene_id"].isin(deeptmhmm_tested_genes)
    qc_after_merge(master, n_before, "merge DeepTMHMM", logger)
    report_source_coverage(master, ["DeepTMHMM_topology_class"], "DeepTMHMM", logger)

    # ---- 6: PredGPI -----------------------------------------------------
    section(logger, "SOURCE 6/6: PredGPI")
    predgpi = parse_predgpi(config["predgpi_up"], config["predgpi_down"], logger)
    results["predgpi_collapsed"] = predgpi
    predgpi_tested_genes = set(predgpi["gene_id"])
    n_before = qc_before_merge(master, "merge PredGPI", logger)
    master = master.merge(predgpi, on="gene_id", how="left")
    master["PredGPI_tested"] = master["gene_id"].isin(predgpi_tested_genes)
    qc_after_merge(master, n_before, "merge PredGPI", logger)
    report_source_coverage(master, ["PredGPI_prediction"], "PredGPI", logger)

    # ---- Cosmetic cleanup: fix float-upcast display artifacts ---------
    # (see clean_integer_display() docstring for why this is needed)
    master["DeepTMHMM_TM_count"] = master["DeepTMHMM_TM_count"].apply(clean_integer_display)
    master["PredGPI_anchor_position"] = master["PredGPI_anchor_position"].apply(clean_integer_display)

    # ---- Final QC -----------------------------------------------------
    section(logger, "FINAL QC")
    if len(master) != n_expected:
        raise RuntimeError(
            f"Final master table has {len(master)} rows; expected {n_expected} "
            f"(the backbone gene count). This must never change across merges."
        )
    if master["gene_id"].duplicated().any():
        raise RuntimeError("Final master table has duplicated gene_id values.")
    if master["gene_id"].isna().any():
        raise RuntimeError("Final master table has null gene_id values.")
    logger.info(f"Final master table: {len(master)} rows, {master['gene_id'].nunique()} unique genes.")
    logger.info(f"Final master table: {len(master.columns)} columns: {list(master.columns)}")

    # Fixed, documented column order for the deliverable.
    column_order = [
        "gene_id", "transcript_id", "protein_id",
        "eggNOG_Description", "eggNOG_Preferred_name", "eggNOG_seed_ortholog",
        "eggNOG_COG_category", "eggNOG_GOs", "eggNOG_KEGG_ko", "eggNOG_KEGG_Pathway",
        "eggNOG_PFAMs",
        "InterPro_tested", "InterPro_accessions", "InterPro_descriptions",
        "InterPro_databases", "InterPro_signature_accessions",
        "InterPro_signature_descriptions", "InterPro_GOs", "InterPro_pathways",
        "Orthogroup", "Candida_albicans_orthologs", "Saccharomyces_cerevisiae_orthologs",
        "Schizosaccharomyces_pombe_orthologs",
        "SignalP_tested", "SignalP_prediction", "SignalP_SP_probability",
        "SignalP_cleavage_site",
        "DeepTMHMM_tested", "DeepTMHMM_has_TM", "DeepTMHMM_TM_count",
        "DeepTMHMM_TM_regions", "DeepTMHMM_topology_class",
        "PredGPI_tested", "PredGPI_prediction", "PredGPI_anchor_position", "PredGPI_score",
    ]
    missing_cols = set(column_order) - set(master.columns)
    if missing_cols:
        raise RuntimeError(f"Final master table is missing expected column(s): {missing_cols}")
    master = master[column_order]

    results["master"] = master
    return results


# =============================================================================
# MAIN
# =============================================================================


def parse_args() -> argparse.Namespace:
    """Parse the two command-line overrides (--project-root, --output-dir); see CONFIGURATION above."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path("/scratch/flljus002"),
        help="Root directory containing all annotation tool outputs.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory to write outputs to (default: <project-root>/PCG_functional_annotation).",
    )
    return parser.parse_args()


def main() -> None:
    """
    Entry point: resolve paths, set up logging, run build_master_table(),
    then save the master table (TSV + CSV) and every collapsed
    intermediate table. Re-raises any exception after logging its
    traceback, so a failed run exits non-zero rather than silently
    producing a partial/incorrect output.
    """
    args = parse_args()
    project_root: Path = args.project_root
    output_dir: Path = args.output_dir or (project_root / "PCG_functional_annotation")
    output_dir.mkdir(parents=True, exist_ok=True)

    logger = setup_logging(output_dir / "PCG_functional_annotation_QC.txt")
    section(logger, "Step 10: build master C. auris PCG functional annotation table")
    logger.info(f"project_root = {project_root}")
    logger.info(f"output_dir   = {output_dir}")

    config = build_config(project_root, output_dir)

    try:
        results = build_master_table(config, logger)
    except Exception:
        logger.exception("Step 10 FAILED -- see traceback above.")
        raise

    section(logger, "SAVING OUTPUTS")
    for name in [
        "eggnog_collapsed", "interpro_collapsed", "orthofinder_collapsed",
        "signalp_collapsed", "deeptmhmm_collapsed", "predgpi_collapsed",
    ]:
        out_path = output_dir / f"{name}.tsv"
        results[name].to_csv(out_path, sep="\t", index=False)
        logger.info(f"Saved {out_path} ({len(results[name])} rows)")

    master_tsv = output_dir / "PCG_functional_annotation.tsv"
    master_csv = output_dir / "PCG_functional_annotation.csv"
    results["master"].to_csv(master_tsv, sep="\t", index=False)
    results["master"].to_csv(master_csv, sep=",", index=False)
    logger.info(f"Saved {master_tsv} ({len(results['master'])} rows)")
    logger.info(f"Saved {master_csv} ({len(results['master'])} rows)")

    section(logger, "Step 10 COMPLETE")
    logger.info(f"QC report written to {output_dir / 'PCG_functional_annotation_QC.txt'}")


if __name__ == "__main__":
    main()
