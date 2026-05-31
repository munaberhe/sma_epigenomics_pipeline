# SMA Epigenomics Pipeline

A Snakemake-based WGBS analysis pipeline for genome-wide epigenetic profiling of nusinersen (ASO) and valproic acid (VPA) combination treatment in Spinal Muscular Atrophy.

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

- **Alignment** — paired-end WGBS alignment to hg38 using Bismark and Bowtie2, with PCR deduplication and CpG report generation
- **SMN1 masking** — hard-masks the SMN1 paralog locus (chr5:70,924,941-70,953,015) with Ns before alignment so reads map unambiguously to SMN2
- **DMR calling** — per-chromosome parallelisation via SLURM array jobs using DMRcaller with locked parameters confirmed with supervisor; results combined into genome-wide GRanges objects
- **Five-contrast 2x2 factorial design** — ASO effect, VPA effect, combination vs baseline, and two cross-background contrasts that test whether the two drugs act independently
- **Genomic annotation** — ChIPseeker annotation of high-confidence DMRs with GO biological process and KEGG pathway enrichment via clusterProfiler, analysed separately for hyper and hypo DMRs
- **MSigDB enrichment** — independent gene set enrichment across neural, synaptic, chromatin and splicing gene sets for all five contrasts, providing a second orthogonal database for pathway validation
- **TF motif enrichment** — monaLisa enrichment of JASPAR2020 vertebrate TF motifs in ASO-specific DMR sequences vs matched background
- **Splice junction proximity** — tests whether ASO-specific DMRs are enriched near exon-intron boundaries relative to matched random background
- **H3K9me2 validation** — overlaps DMR loci with published HEK293T H3K9me2 ChIP-seq bigWig signal (Marasco et al. 2022, GSE167762) to validate kinetic coupling model at SMN2
- **UpSet overlap analysis** — identifies DMRs shared across contrasts and isolates ASO-specific high-confidence DMR set
- **Low-resolution genome browser tracks** — sliding-window methylation profiles for chr1, chrX and chr5/SMN2 zoom showing all four conditions
- **QC suite** — 12-sample PCA from per-replicate chr1 CpG methylation, M-bias plots, duplication rates, bisulfite conversion efficiency, sample correlation heatmap, coverage retention curves

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
| MSigDB enrichment | msigdbr | — |
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
- conda (docs.conda.io)
- R 4.5.1+ available via module or conda
- ~3TB scratch storage for intermediate files
- ~500GB storage for final outputs

---

## Installation

Clone the repository and create the conda environment:

    git clone https://github.com/munaberhe/sma_epigenomics_pipeline.git
    cd sma_epigenomics_pipeline
    conda env create -f environment.yml
    conda activate sma_epigenomics_pipeline

Install R packages:

    Rscript -e "
    BiocManager::install(c(
      'DMRcaller', 'GenomicRanges', 'ChIPseeker', 'clusterProfiler',
      'org.Hs.eg.db', 'TxDb.Hsapiens.UCSC.hg38.knownGene',
      'monaLisa', 'SummarizedExperiment', 'txdbmaker',
      'msigdbr', 'rtracklayer', 'UpSetR', 'patchwork', 'ggplot2'
    ))"

Download hg38 and build the SMN1-masked Bismark index:

    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -P data/reference/
    gunzip data/reference/hg38.fa.gz
    sbatch scripts/01_mask_and_index.sh

---

## Configuration

Edit configs/config.yaml to set sample names, paths and SLURM resource profiles.

Key environment variables used by the Snakefile:

| Variable | Description |
|---|---|
| RSCRIPT | Path to Rscript binary |
| R_LIBS | Path to R library directory |

These are set automatically from the conda environment. Override in configs/config.yaml if needed.

---

## Running the Pipeline

Full pipeline via Snakemake:

    snakemake --profile configs/slurm_profile --jobs 12

Dry run to check the DAG without executing:

    snakemake -n

Step-by-step execution:

    # 1. Alignment (unmasked and masked)
    sbatch scripts/02_bismark_align.sh
    sbatch scripts/03_dedup_and_extract.sh
    sbatch scripts/04_split_by_chr.sh

    # 2. DMR calling
    bash scripts/submit_dmr_by_chr.sh
    Rscript scripts/06b_dmrcaller_combine_chr.R

    # 3. SMN2 locus analysis
    Rscript scripts/05_smn2_locus_final.R
    Rscript scripts/05b_smn_locus_unmasked.R

    # 4. Annotation and enrichment
    Rscript scripts/07_dmr_annotate.R
    Rscript scripts/07b_dmr_plots.R
    Rscript scripts/07c_dmr_locus_plots.R
    Rscript scripts/09_top10_dmrs.R

    # 5. QC
    Rscript scripts/08_pca_chr1.R
    Rscript scripts/08b_additional_qc.R
    Rscript scripts/10_coverage_qc.R

    # 6. Downstream analyses
    Rscript scripts/11_h3k9me2_overlap.R
    Rscript scripts/12_upset_dmr_intersections.R
    Rscript scripts/14_msigdb_enrichment_all.R
    Rscript scripts/15_lowres_methylation_profile.R
    Rscript scripts/16_tf_motif_enrichment.R
    Rscript scripts/16b_tf_motif_plots.R
    Rscript scripts/17_splice_junction_proximity.R

---

## Key Parameters

DMR calling parameters confirmed with Prof Radu Zabet, 5 May 2026.

| Parameter | Value | Rationale |
|---|---|---|
| method | bins | Fixed-width windows, robust at pooled coverage |
| binSize | 300 bp | Benchmarked against permutation null on chr1/chr6/chr13 |
| minProportionDifference | 0.20 | Filters biologically trivial changes |
| pValueThreshold | 0.01 | Standard genome-wide threshold |
| minCytosinesCount | 4 | Prevents single-CpG noise calls |
| minReadsPerCytosine | 4 | Confirmed from benchmark scripts |
| minGap | 300 bp | One bin width -- prevents iterative merge hang on dense VPA contrasts |
| test | score | Rao score test, appropriate at ~27x pooled coverage |
| context | CG | CpG only -- CHG/CHH near-zero in human somatic cells |

regionType convention: gain = hypomethylated (proportion1 < proportion2); loss = hypermethylated. This is counter-intuitive -- documented explicitly to prevent misinterpretation.

Per-chromosome DMR calling was adopted after genome-wide calling proved computationally intractable. Attempts with up to 98GB RAM ran for 18-27 hours without completing the iterative merge step (SLURM jobs 10504170, 10591016).

---

## Project Structure

    sma_epigenomics_pipeline/
    |-- Snakefile                            Full pipeline DAG
    |-- README.md
    |-- environment.yml                      Conda environment
    |-- configs/
    |   |-- config.yaml                      Pipeline configuration
    |   `-- slurm_profile/                   Snakemake SLURM profile
    |-- scripts/
    |   |-- 01_mask_and_index.sh             Mask SMN1 locus, build Bismark index
    |   |-- 02_bismark_align.sh              WGBS alignment (SLURM array)
    |   |-- 03_dedup_and_extract.sh          Deduplication + methylation extraction
    |   |-- 04_split_by_chr.sh               Split CX reports by chromosome
    |   |-- 05_smn2_locus_final.R            SMN2 locus methylation plots (masked)
    |   |-- 05b_smn_locus_unmasked.R         SMN2 locus plots (unmasked)
    |   |-- 06_dmrcaller_by_chr.R            Per-chromosome DMR calling
    |   |-- 06b_dmrcaller_combine_chr.R      Combine per-chromosome results
    |   |-- 07_dmr_annotate.R                ChIPseeker annotation + GO/KEGG
    |   |-- 07b_dmr_plots.R                  Per-chr bar charts + meth diff histograms
    |   |-- 07c_dmr_locus_plots.R            Annotated locus overlay plots
    |   |-- 08_pca_chr1.R                    12-sample PCA from chr1 CpG methylation
    |   |-- 08b_additional_qc.R              M-bias, duplication rates, conversion QC
    |   |-- 09_top10_dmrs.R                  Top 10 hypo DMRs per contrast
    |   |-- 10_coverage_qc.R                 Coverage retention curves
    |   |-- 11_h3k9me2_overlap.R             H3K9me2 signal enrichment at DMR loci
    |   |-- 12_upset_dmr_intersections.R     UpSet overlap plots
    |   |-- 14_msigdb_enrichment_all.R       MSigDB gene set enrichment all contrasts
    |   |-- 15_lowres_methylation_profile.R  Low-res browser tracks
    |   |-- 16_tf_motif_enrichment.R         TF motif enrichment (monaLisa)
    |   |-- 16b_tf_motif_plots.R             TF motif result plots
    |   |-- 17_splice_junction_proximity.R   Splice junction proximity test
    |   |-- submit_dmr_by_chr.sh             Submit DMR SLURM array
    |   |-- submit_smn1_masked_pipeline.sh   Submit masked alignment pipeline
    |   `-- archive/                         Superseded scripts retained for reference
    |-- data/
    |   |-- reference/                       hg38 FASTA + Bismark index
    |   `-- reference_smn1_masked/           SMN1-masked reference + index
    |-- results/                             Pipeline outputs (gitignored)
    `-- logs/                                SLURM logs (gitignored)

---

## References

- Marasco et al. (2022) Cell 185:2057-2070 -- nusinersen kinetic coupling + H3K9me2 at SMN2
- Catoni et al. (2018) Nucleic Acids Research 46:e114 -- DMRcaller
- Krueger and Andrews (2011) Bioinformatics 27:1571-1572 -- Bismark
- Yu et al. (2015) OMICS 19:284-287 -- ChIPseeker
- Wu et al. (2021) iMeta 1:e5 -- clusterProfiler 4.0
- Finkel et al. (2017) NEJM 377:1723-1732 -- ENDEAR trial nusinersen
- Gottlicher et al. (2001) EMBO J 20:6969-6978 -- VPA as HDAC inhibitor
