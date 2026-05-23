# SMA Epigenomics Pipeline

A WGBS analysis pipeline for assessing genome-wide pleiotropic epigenetic effects of combined nusinersen (ASO) and valproic acid (VPA) treatment in Spinal Muscular Atrophy.
Built as an MSc Bioinformatics thesis project, Queen Mary University of London, 2026.

> For a full record of parameter decisions and rationale see FINAL_SCRIPTS.md

## Table of Contents

- [Background](#background)
- [Experiment Design](#experiment-design)
- [Pipeline Structure](#pipeline-structure)
- [Tech Stack](#tech-stack)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Running the Pipeline](#running-the-pipeline)
- [Key Parameters](#key-parameters)
- [Repository Layout](#repository-layout)
- [Key References](#key-references)

---

## Background

SMA is caused by loss-of-function mutations in SMN1. The SMN2 paralog compensates partially but skips exon 7 in most transcripts, producing very little functional protein. Nusinersen (Spinraza) corrects this splicing defect via antisense oligonucleotide targeting. Marasco et al. (2022, Cell) showed that nusinersen also deposits the repressive mark H3K9me2 at the SMN2 locus. Valproic acid (VPA), a broad HDAC inhibitor, counteracts this chromatin compaction and the combined ASO+VPA treatment outperforms ASO alone in patient fibroblasts and mouse models.

This pipeline investigates the genome-wide epigenetic effects of this combination therapy — specifically, whether VPA has off-target methylation effects beyond the SMN locus.

---

## Experiment Design

2x2 factorial WGBS experiment in motor neuron-like cells.

| Condition | Replicates | Description |
|-----------|-----------|-------------|
| ASO_CTRL | 3 | Nusinersen only (100 nM, saturating dose) |
| ASO_VPA | 3 | Nusinersen + VPA (combination treatment) |
| Scramble_CTRL | 3 | Scramble ASO — baseline control |
| Scramble_VPA | 3 | VPA only |

Three contrasts: ASO_VPA vs Scramble_CTRL (primary), Scramble_VPA vs Scramble_CTRL (VPA alone), ASO_CTRL vs Scramble_CTRL (negative control).

---

## Pipeline Structure

Raw FASTQ files are quality-trimmed with Trim Galore, then aligned to an SMN1-masked hg38 reference using Bismark. SMN1 is hard-masked with Ns so reads from both SMN1 and SMN2 map unambiguously to SMN2, resolving the paralog dropout problem. Aligned BAMs are deduplicated, methylation-extracted, and a genome-wide cytosine report is generated via coverage2cytosine. CX reports are split by chromosome (CpG only) for memory-efficient DMR calling. DMRcaller runs per-chromosome in parallel SLURM jobs and results are combined. Annotations use ChIPseeker and pathway enrichment uses clusterProfiler.

---

## Tech Stack

| Component | Tool |
|-----------|------|
| Alignment | Bismark v0.25.1 + Bowtie2 |
| Reference | hg38 (GRCh38), SMN1-masked |
| Deduplication | deduplicate_bismark |
| Methylation extraction | bismark_methylation_extractor + coverage2cytosine |
| DMR calling | DMRcaller (Bioconductor) |
| Annotation | ChIPseeker, TxDb.Hsapiens.UCSC.hg38.knownGene |
| Enrichment | clusterProfiler (GO + KEGG) |
| Workflow | SLURM array jobs (Apocrita HPC, QMUL) |
| Language | R 4.5.1, bash |

---

## Prerequisites

- Apocrita HPC account (QMUL) or equivalent SLURM cluster
- conda (docs.conda.io)
- R 4.5.1+ available via module load or conda
- ~3TB scratch storage for intermediate files
- ~2TB external storage for BAM backups

---

## Installation

Clone the repository and set up the conda environment:

    git clone https://github.com/munaberhe/sma_epigenomics_pipeline.git
    cd sma_epigenomics_pipeline
    conda env create -f environment.yml
    conda activate sma_epigenomics_pipeline

Install R packages:

    install.packages("BiocManager")
    BiocManager::install(c("DMRcaller","GenomicRanges","ChIPseeker",
      "clusterProfiler","org.Hs.eg.db",
      "TxDb.Hsapiens.UCSC.hg38.knownGene","R.utils"))

Download hg38 and build the SMN1-masked Bismark index:

    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -P data/reference/
    gunzip data/reference/hg38.fa.gz
    sbatch scripts/01_mask_and_index.sh

---

## Running the Pipeline

Step 1 — Alignment (all 12 samples as SLURM array):

    sbatch scripts/02_bismark_align.sh

Step 2 — Deduplication and methylation extraction:

    sbatch scripts/03_dedup_and_extract.sh

Step 3 — Split CX reports by chromosome:

    sbatch --dependency=afterok:<METH_JOB> scripts/04_split_by_chr.sh

Step 4 — Per-chromosome DMR calling (all 3 contrasts in parallel):

    bash scripts/submit_dmr_by_chr.sh
    Rscript scripts/dmrcaller_combine_chr.R

Step 5 — SMN2 locus masked analysis:

    Rscript scripts/dmrcaller_smn_locus_masked.R

Step 6 — Annotation and enrichment:

    Rscript scripts/dmr_annotate.R

---

## Key Parameters

All DMR parameters confirmed with supervisor Dr Radu Zabet, 5 May 2026.

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| method | bins | Fixed-width windows, robust at ~27x pooled coverage |
| binSize | 300 bp | Benchmarked against chr1 permutation null |
| minProportionDifference | 0.20 | Filters biologically trivial changes |
| pValueThreshold | 0.01 | Standard genome-wide threshold |
| minCytosinesCount | 4 | Prevents single-CpG noise calls |
| minReadsPerCytosine | 4 | Confirmed from benchmark scripts |
| minGap | 300 | One bin width — prevents merge hang on VPA contrasts |
| test | score | Rao test, appropriate at ~27x pooled coverage |
| context | CG | CpG only — CHG/CHH near-zero in human somatic cells |

See FINAL_SCRIPTS.md for full decision rationale.

---

## Repository Layout

    sma_epigenomics_pipeline/
    scripts/
        01_mask_and_index.sh              Mask SMN1, build Bismark index
        02_bismark_align.sh               WGBS alignment (SLURM array)
        03_dedup_and_extract.sh           Dedup + methylation extraction
        04_split_by_chr.sh                Split CX reports by chromosome
        submit_dmr_by_chr.sh              Submit per-chromosome DMR jobs
        dmrcaller_by_chr.R                Per-chromosome DMR calling
        dmrcaller_combine_chr.R           Combine per-chromosome results
        dmrcaller_genome_wide.R           Genome-wide DMR calling (single job)
        dmrcaller_smn_locus_masked.R      SMN2 locus profile (masked data)
        smn_locus_dmrcaller_comparisons.R SMN1/2 unmasked baseline
        coverage_4lines_per_condition.R   CpG coverage QC
        dmr_annotate.R                    DMR annotation + GO/KEGG enrichment
        check_scratch.sh                  Scratch space pre-flight check
    data/
        reference/                        hg38 FASTA
        reference_smn1_masked/            SMN1-masked reference + Bismark index
        processed/                        Trimmed FASTQ files
    results/
        alignments_smn1_masked/           Masked BAMs, dedup, CX reports, by_chr
        alignments/                       Original unmasked by_chr files
        dmr/                              DMR RDS, BED, summary files
        qc/                               Coverage plots, SMN locus plots
    logs/                                 SLURM and tool logs
    FINAL_SCRIPTS.md                      Parameter rationale for all scripts
    README.md

---

## Key References

- Marasco et al. (2022) Cell 185:2057-2070 — nusinersen + H3K9me2 at SMN2
- Catoni et al. (2018) Nucleic Acids Research 46:e114 — DMRcaller
- Kornblihtt et al. — gene looping and splicing regulation at SMN2
- Krueger & Andrews (2011) Bioinformatics — Bismark
- Yu et al. (2015) OMICS — ChIPseeker

---

Muna Berhe · bt25018 · MSc Bioinformatics, QMUL 2025-2026
Supervisor: Dr Radu Zabet
Collaborators: Prof Alberto Kornblihtt (IFIBYNE-UBA-CONICET), Dr Emilia Haberfeld, Dr Marcos Miretti
