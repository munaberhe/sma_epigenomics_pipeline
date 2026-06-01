# SMA Epigenomics Pipeline

Snakemake WGBS pipeline for genome-wide epigenetic profiling of nusinersen (ASO1) and valproic acid (VPA) combination treatment in Spinal Muscular Atrophy.

MSc Bioinformatics thesis · Queen Mary University of London · 2026
Supervisor: Prof Radu Zabet (Zabet Lab, QMUL)
Collaborators: Prof Alberto Kornblihtt (IFIBYNE-UBA-CONICET), Dr Emilia Haberfeld, Dr Marcos Miretti

---

## Table of Contents

- [Features](#features)
- [Tech Stack](#tech-stack)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Configuration](#configuration)
- [Running the Pipeline](#running-the-pipeline)
- [Key Parameters](#key-parameters)
- [Project Structure](#project-structure)
- [References](#references)

---

## Features

- **Alignment** -- paired-end WGBS to hg38 via Bismark/Bowtie2, with PCR deduplication and CpG report generation
- **SMN1 masking** -- hard-masks the SMN1 paralog locus (chr5:70,924,941-70,953,015) with Ns so reads map unambiguously to SMN2
- **DMR calling** -- per-chromosome SLURM array jobs using DMRcaller; parameters benchmarked on chr1/chr6/chr13 permutation null; results combined into genome-wide GRanges objects
- **2x2 factorial design, five contrasts** -- ASO effect, VPA effect, combination vs baseline, plus two cross-background contrasts to test whether the drugs act independently
- **Genomic annotation** -- ChIPseeker annotation with GO/KEGG enrichment via clusterProfiler, run separately for hyper- and hypomethylated DMRs
- **MSigDB enrichment** -- gene set enrichment against neural, synaptic, chromatin and splicing collections for all five contrasts; used as an independent check on clusterProfiler GO results
- **TF motif enrichment** -- monaLisa/JASPAR2020 motif enrichment in ASO-specific DMR sequences vs GC-matched background
- **Splice junction proximity** -- Wilcoxon test for enrichment of ASO-specific DMRs near exon-intron boundaries
- **H3K9me2 validation** -- overlaps DMR loci with HEK293T H3K9me2 ChIP-seq signal (Marasco et al. 2022, GSE167762) to check the kinetic coupling model at SMN2
- **UpSet overlap** -- identifies DMRs shared across contrasts and isolates the ASO-specific DMR set
- **Low-resolution methylation profiles** -- sliding-window profiles for chr1, chrX and SMN2 zoom, all four conditions
- **QC** -- 12-sample PCA, M-bias plots, duplication rates, bisulfite conversion, correlation heatmap, coverage retention curves

---

## Tech Stack

| Component | Tool | Version |
|---|---|---|
| Trimming | Trim Galore | v0.6.11 |
| Alignment | Bismark | v0.25.1 |
| Aligner | Bowtie2 | v2.5.4 |
| Reference genome | hg38 GRCh38 Ensembl 109 | SMN1-masked |
| Deduplication | deduplicate_bismark | v0.25.1 |
| DMR calling | DMRcaller (Bioconductor) | v1.42.0 |
| Annotation | ChIPseeker | v1.46.1 |
| Pathway enrichment | clusterProfiler | v4.18.4 |
| MSigDB enrichment | msigdbr | -- |
| TF motif | monaLisa + JASPAR2020 | v1.16.0 |
| H3K9me2 validation | rtracklayer bigWig | GSE167762 |
| Visualisation | ggplot2, patchwork, UpSetR | v4.0.2 |
| Workflow manager | Snakemake | v9.17.2 |
| HPC scheduler | SLURM | v24.11.7 |
| Language | R | v4.5.1 |
| OS | Rocky Linux 9 | Apocrita HPC QMUL |

---

## Prerequisites

- Apocrita HPC account (QMUL) or equivalent SLURM cluster
- conda
- R 4.5.1+ via module or conda
- ~3TB scratch storage for intermediates
- ~500GB for final outputs

---

## Installation

```bash
git clone https://github.com/munaberhe/sma_epigenomics_pipeline.git
cd sma_epigenomics_pipeline
conda env create -f environment.yml
conda activate sma_epigenomics_pipeline
```

Install R packages:

```r
BiocManager::install(c(
  "DMRcaller", "GenomicRanges", "ChIPseeker", "clusterProfiler",
  "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "monaLisa", "SummarizedExperiment", "txdbmaker",
  "msigdbr", "rtracklayer", "UpSetR", "patchwork", "ggplot2"
))
```

Download hg38 and build the SMN1-masked Bismark index:

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -P data/reference/
gunzip data/reference/hg38.fa.gz
sbatch scripts/01_mask_and_index.sh
```

---

## Configuration

Edit `configs/config.yaml` to set sample names, paths and SLURM resource profiles.

| Variable | Description |
|---|---|
| RSCRIPT | Path to Rscript binary |
| R_LIBS | Path to R library directory |

Both are set automatically from the conda environment.

---

## Running the Pipeline

Full pipeline:

```bash
snakemake --profile configs/slurm_profile --jobs 12
```

Dry run:

```bash
snakemake -n
```

Step-by-step:

```bash
# 1. Alignment
sbatch scripts/02_bismark_align.sh
sbatch scripts/03_dedup_and_extract.sh
sbatch scripts/04_split_by_chr.sh

# 2. DMR calling
bash scripts/submit_dmr_by_chr.sh
Rscript scripts/06b_dmrcaller_combine_chr.R

# 3. SMN2 locus
Rscript scripts/05_smn2_locus_final.R

# 4. Annotation
Rscript scripts/07_dmr_annotate.R
Rscript scripts/07b_dmr_plots.R
Rscript scripts/09_top10_dmrs.R

# 5. QC
Rscript scripts/08_pca_chr1.R
Rscript scripts/08b_additional_qc.R
Rscript scripts/10_coverage_qc.R

# 6. Downstream
Rscript scripts/11_h3k9me2_overlap.R
Rscript scripts/12_upset_dmr_intersections.R
Rscript scripts/14_msigdb_enrichment_all.R
Rscript scripts/16_tf_motif_enrichment.R
Rscript scripts/17_splice_junction_proximity.R
```

---

## Key Parameters

Parameters used in this study, benchmarked against permutation null distributions on chr1, chr6 and chr13. For other datasets, binSize and minProportionDifference should be re-benchmarked using the permutation approach in scripts/archive/.

| Parameter | Value | Rationale |
|---|---|---|
| method | bins | Fixed-width windows, robust at pooled coverage |
| binSize | 300 bp | Largest signal-to-noise gap vs permutation null |
| minProportionDifference | 0.20 | Removes biologically marginal changes |
| pValueThreshold | 0.01 | Standard genome-wide threshold |
| minCytosinesCount | 4 | Prevents single-CpG noise calls |
| minReadsPerCytosine | 4 | Confirmed from benchmark scripts |
| minGap | 300 bp | Prevents iterative merge hang on dense VPA contrasts |
| test | score | Rao score test at ~27x pooled coverage |
| context | CG | CHG/CHH near-zero in human somatic cells |

Note: regionType gain = hypomethylated (proportion1 < proportion2); loss = hypermethylated.

Per-chromosome calling was necessary after genome-wide attempts timed out at 98GB RAM after 18-27 hours (SLURM jobs 10504170, 10591016).

---

## Project Structure

```
sma_epigenomics_pipeline/
├── Snakefile
├── README.md
├── environment.yml
├── configs/
│   ├── config.yaml
│   └── slurm_profile/
├── scripts/
│   ├── 01_mask_and_index.sh
│   ├── 02_bismark_align.sh
│   ├── 03_dedup_and_extract.sh
│   ├── 04_split_by_chr.sh
│   ├── 05_smn2_locus_final.R
│   ├── 06_dmrcaller_by_chr.R
│   ├── 06b_dmrcaller_combine_chr.R
│   ├── 07_dmr_annotate.R
│   ├── 07b_dmr_plots.R
│   ├── 08_pca_chr1.R
│   ├── 08b_additional_qc.R
│   ├── 09_top10_dmrs.R
│   ├── 10_coverage_qc.R
│   ├── 11_h3k9me2_overlap.R
│   ├── 12_upset_dmr_intersections.R
│   ├── 14_msigdb_enrichment_all.R
│   ├── 16_tf_motif_enrichment.R
│   ├── 17_splice_junction_proximity.R
│   ├── 19_smn2_h3k27ac_enhancer.R
│   ├── 20_dss_replicate_testing.R
│   ├── submit_dmr_by_chr.sh
│   └── archive/
├── data/
│   ├── reference/
│   └── reference_smn1_masked/
├── results/          (gitignored)
└── logs/             (gitignored)
```

---

## References

- Marasco et al. (2022) Cell 185:2057-2070 -- nusinersen, H3K9me2, kinetic coupling
- Catoni et al. (2018) Nucleic Acids Research 46:e114 -- DMRcaller
- Krueger and Andrews (2011) Bioinformatics 27:1571-1572 -- Bismark
- Yu et al. (2015) OMICS 19:284-287 -- ChIPseeker
- Wu et al. (2021) iMeta 1:e5 -- clusterProfiler 4.0
- Park and Wu (2016) Bioinformatics 32:1414-1416 -- DSS
- Finkel et al. (2017) NEJM 377:1723-1732 -- nusinersen clinical trial
- Gottlicher et al. (2001) EMBO J 20:6969-6978 -- VPA as HDAC inhibitor
