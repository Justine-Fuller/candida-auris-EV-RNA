# *Candida auris* extracellular vesicle RNA and host RNA-binding proteins

Analysis code and supplementary tables for the BSc (Med) Honours thesis of Justine Fuller, Department of Molecular and Cell Biology, University of Cape Town (2026).

The scripts are organised by thesis section. They were run on the UCT HPC cluster (SLURM `.sbatch` jobs) and on a laptop (R and Python). File paths inside the scripts point to the original locations, so they need to be changed before the scripts are rerun elsewhere.

## Data used

| Data | Source |
|---|---|
| *C. auris* EV and cell RNA-seq (strains B8441 and MMC1, with and without caspofungin) | NCBI SRA SRP295539 (Munhoz da Rocha et al.) |
| *C. auris* genome and annotation | Candida Genome Database, B8441 version s02-m01-r01 |
| Human PBMCs exposed to *C. auris* (RNA-seq) | NCBI GEO GSE154911 (Bruno et al., 2020) |
| eCLIP peaks for human RNA-binding proteins (hg19) | ENCODE; accessions in `supplementary_tables/Supplementary_Table_S8_eCLIP_accessions.csv` |

Raw and processed sequencing files are not included because of their size.

## Folders

| Folder | Thesis section | What it does |
|---|---|---|
| `01_read_census` | 3.1 | Counts reads from each RNA class in EV and cell libraries |
| `02_cellular_DE` | 3.1, S2, S3 | Cellular mRNA differential expression (DESeq2), strand-aware recount, GO/KEGG enrichment and GSEA, lncRNA–gene pairs, protein annotation (SignalP, DeepTMHMM, PredGPI, InterProScan) |
| `03_RNA_abundance_and_full_length` | 3.2 | RNA class abundance and full-length versus fragment classification |
| `04_EV_vs_cell_bins` | 3.3 | Coverage in 100-nt bins and EV versus cell enrichment tests |
| `05_tRNA` | 3.4 | tRNA counting, 3′ CCA check and EV versus cell tRNA families |
| `06_sequence_tags` | 3.5 | k-mer, GGAC/DRACH and motif tests, and the 85 EV-selective regions |
| `07_codon_usage` | 3.6 | Codon usage scoring and null sets |
| `08_host_splicing` | 3.7 | LeafCutter differential splicing in PBMCs |
| `09_RBP_binding` | 3.7 | eCLIP peak overlap with spliced host genes and EV RNAs, FIMO motifs and RNA structure |
| `10_figures` | all | Builds the thesis figures |
| `supplementary_tables` | S3, S6 | Tables referred to in the thesis |

## Supplementary tables

- `EV_selective_regions_85.tsv` and `EV_selective_regions_85_sequences.fa`: the 85 EV-selective regions (Table S6.2).
- `lncRNA_gene_pairs_concordant_both_strains_164.csv`: the 164 lncRNA–protein-coding gene pairs that change in the same direction in both strains (section S3).
- `Supplementary_Table_S8_eCLIP_accessions.csv`: ENCODE eCLIP files used.

## Main software

R (DESeq2, edgeR, clusterProfiler), Python (pandas, scipy, statsmodels, gseapy, pyBigWig), STAR, samtools, bedtools, cutadapt, LeafCutter, MEME suite (STREME, FIMO), ViennaRNA, SignalP 6, DeepTMHMM 1.0.24, PredGPI, InterProScan and eggNOG-mapper.
