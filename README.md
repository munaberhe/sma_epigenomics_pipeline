# SMA Epigenomics Pipeline

Whole-genome bisulfite sequencing (WGBS) analysis pipeline investigating genome-wide pleiotropic epigenetic effects of combined nusinersen (ASO1) and valproic acid (VPA) treatment in Spinal Muscular Atrophy.

MSc Bioinformatics thesis · Queen Mary University of London · 2026  
Supervisor: Prof Radu Zabet (Zabet Lab, QMUL)  
Collaborators: Prof Alberto Kornblihtt (IFIBYNE-UBA-CONICET), Dr Emilia Haberfeld, Dr Marcos Miretti

---

## Background

SMA is caused by loss-of-function mutations in SMN1. The SMN2 paralog compensates partially but skips exon 7 in ~90% of transcripts due to a C→T transition at position +6, producing insufficient functional SMN protein. Nusinersen (Spinraza) corrects this splicing defect via antisense oligonucleotide targeting of ISS-N1 in intron 7. Marasco et al. (2022, Cell) showed that nusinersen also deposits the repressive histone mark H3K9me2 at the SMN2 locus via a kinetic coupling mechanism. Valproic acid (VPA), a broad HDAC inhibitor, counteracts this chromatin compaction and combined ASO+VPA treatment outperforms ASO alone in preclinical models.

This pipeline characterises the genome-wide epigenetic consequences of this combination therapy in HEK293T cells, testing whether VPA and ASO act through independent epigenetic mechanisms and whether nusinersen produces off-target methylation changes at neural identity loci.

---

## Experimental Design

2x2 factorial WGBS experiment in HEK293T cells, n=3 biological replicates per condition.

| Condition | Replicates | Description |
|---|---|---|
| ASO_CTRL | 3 | Nusinersen 100nM (saturating dose) |
| ASO_VPA | 3 | Nusinersen + VPA (combination treatment) |
| Scramble_CTRL | 3 | Scramble ASO — baseline control |
| Scramble_VPA | 3 | VPA only (HDAC inhibitor) |

Five contrasts were analysed to fully characterise the 2x2 factorial design.

| Contrast | DMRs | Key finding |
|---|---|---|
| ASO_CTRL vs Scramble_CTRL | 3,423 | Bidirectional, chrX hotspot, neural pathways |
| Scramble_VPA vs Scramble_CTRL | 598,485 | Near-exclusively hypo, HDAC inhibitor signature |
| ASO_VPA vs Scramble_CTRL | 554,291 | VPA-dominated, same as VPA alone |
| ASO_VPA vs ASO_CTRL | 664,202 | VPA effect identical on ASO background |
| ASO_VPA vs Scramble_VPA | 23,669 | ASO neural signal persists on VPA background |

The parallel VPA profiles in contrasts 2 and 4 confirm that VPA and ASO act through completely independent epigenetic mechanisms.

---

## Pipeline Structure

    Raw FASTQ
        |
        v  Trim Galore (Q20, min 20bp, FastQC)
    Trimmed FASTQ
        |
        |---> Bismark + Bowtie2 -> unmasked hg38 (genome-wide DMR calling)
        |---> Bismark + Bowtie2 -> SMN1-masked hg38 (SMN2 locus analysis)
               chr5:70,924,941-70,953,015 replaced with Ns
        |
        v  deduplicate_bismark + bismark_methylation_extractor
    CpG reports (per sample, genome-wide)
        |
        v  Split by chromosome (CpG context only)
    Per-chromosome CpG reports (12 samples x 24 chromosomes)
        |
        v  DMRcaller per-chromosome (SLURM array)
        v  Combine chromosomes
    Genome-wide DMRs (5 contrasts, RDS + BED format)
        |
        |---> ChIPseeker annotation + GO/KEGG (clusterProfiler)
        |---> MSigDB enrichment (neural/synaptic/chromatin/splicing)
        |---> TF motif enrichment (monaLisa + JASPAR2020) -- negative result
        |---> Splice junction proximity analysis -- negative result
        |---> H3K9me2 overlap (Marasco et al. 2022, GSE167762)
        |---> UpSet overlap analysis (151 ASO-specific DMRs identified)
        |---> Low-resolution genome browser tracks (chr1, chrX, chr5/SMN2)

**SMN1 masking:** SMN1 and SMN2 share ~99% sequence identity causing paralog read misassignment. The SMN1 locus (chr5:70,924,941-70,953,015, 28,075 bp) is hard-masked with Ns so reads align unambiguously to SMN2.

**Per-chromosome DMR calling:** Genome-wide calling was computationally intractable. Attempts with up to 98GB RAM ran for 18-27 hours without completing the iterative merge step (SLURM jobs 10504170, 10591016). Per-chromosome parallelisation resolves this by giving DMRcaller a smaller merge problem per job.

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

## DMR Calling Parameters

All parameters confirmed with Prof Radu Zabet, 5 May 2026.

| Parameter | Value | Rationale |
|---|---|---|
| method | bins | Fixed-width windows, robust at ~27x pooled coverage |
| binSize | 300 bp | Benchmarked against chr1/chr6/chr13 permutation null |
| minProportionDifference | 0.20 | Filters biologically trivial changes |
| pValueThreshold | 0.01 | Standard genome-wide threshold |
| minCytosinesCount | 4 | Prevents single-CpG noise calls |
| minReadsPerCytosine | 4 | Confirmed from benchmark scripts |
| minGap | 300 bp | One bin width -- prevents iterative merge hang on VPA contrasts |
| test | score | Rao score test, appropriate at ~27x pooled coverage |
| context | CG | CpG only -- CHG/CHH near-zero in human somatic cells |

regionType convention: gain = hypomethylated (proportion1 < proportion2); loss = hypermethylated. This is counter-intuitive and opposite to some published conventions. All scripts and figures use this convention consistently.

---

## Installation

    git clone https://github.com/munaberhe/sma_epigenomics_pipeline.git
    cd sma_epigenomics_pipeline
    conda env create -f environment.yml
    conda activate sma_epigenomics_pipeline

Install R packages:

    BiocManager::install(c(
      "DMRcaller", "GenomicRanges", "ChIPseeker", "clusterProfiler",
      "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg38.knownGene",
      "monaLisa", "SummarizedExperiment", "txdbmaker",
      "msigdbr", "rtracklayer", "UpSetR", "patchwork", "ggplot2"
    ))

Download hg38 and build the SMN1-masked Bismark index:

    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -P data/reference/
    gunzip data/reference/hg38.fa.gz
    sbatch scripts/01_mask_and_index.sh

---

## Running the Pipeline

Full pipeline via Snakemake:

    snakemake --profile configs/slurm_profile --jobs 12

Dry run:

    snakemake -n

Step-by-step:

    # Alignment
    sbatch scripts/02_bismark_align.sh
    sbatch scripts/03_dedup_and_extract.sh
    sbatch scripts/04_split_by_chr.sh

    # DMR calling
    bash scripts/submit_dmr_by_chr.sh
    Rscript scripts/06b_dmrcaller_combine_chr.R

    # SMN2 locus
    Rscript scripts/05_smn2_locus_final.R
    Rscript scripts/05b_smn_locus_unmasked.R

    # Annotation and plots
    Rscript scripts/07_dmr_annotate.R
    Rscript scripts/07b_dmr_plots.R
    Rscript scripts/07c_dmr_locus_plots.R
    Rscript scripts/09_top10_dmrs.R

    # QC
    Rscript scripts/08_pca_chr1.R
    Rscript scripts/08b_additional_qc.R
    Rscript scripts/10_coverage_qc.R

    # Downstream analyses
    Rscript scripts/11_h3k9me2_overlap.R
    Rscript scripts/12_upset_dmr_intersections.R
    Rscript scripts/14_msigdb_enrichment_all.R
    Rscript scripts/15_lowres_methylation_profile.R
    Rscript scripts/16_tf_motif_enrichment.R
    Rscript scripts/16b_tf_motif_plots.R
    Rscript scripts/17_splice_junction_proximity.R

---

## Repository Layout

    sma_epigenomics_pipeline/
    ├── Snakefile                            Full pipeline DAG
    ├── README.md
    ├── environment.yml                      Conda environment specification
    ├── configs/
    │   └── config.yaml                      Pipeline configuration
    ├── scripts/
    │   ├── 01_mask_and_index.sh             Mask SMN1 locus, build Bismark index
    │   ├── 02_bismark_align.sh              WGBS alignment (SLURM array, 12 samples)
    │   ├── 03_dedup_and_extract.sh          PCR deduplication + methylation extraction
    │   ├── 04_split_by_chr.sh               Split CX reports by chromosome (CpG only)
    │   ├── 05_smn2_locus_final.R            SMN2 locus methylation plots (masked)
    │   ├── 05b_smn_locus_unmasked.R         SMN2 locus plots (unmasked, for comparison)
    │   ├── 06_dmrcaller_by_chr.R            Per-chromosome DMR calling
    │   ├── 06b_dmrcaller_combine_chr.R      Combine per-chromosome results
    │   ├── 07_dmr_annotate.R                ChIPseeker annotation + GO/KEGG enrichment
    │   ├── 07b_dmr_plots.R                  Per-chr bar charts + methylation diff histograms
    │   ├── 07c_dmr_locus_plots.R            Annotated locus overlay plots
    │   ├── 08_pca_chr1.R                    12-sample PCA from chr1 CpG methylation
    │   ├── 08b_additional_qc.R              M-bias, duplication rates, conversion efficiency
    │   ├── 09_top10_dmrs.R                  Top 10 hypo DMRs per contrast
    │   ├── 10_coverage_qc.R                 Coverage retention curves
    │   ├── 11_h3k9me2_overlap.R             H3K9me2 ChIP-seq signal at DMR loci
    │   ├── 12_upset_dmr_intersections.R     UpSet plots (3-contrast + 5-contrast)
    │   ├── 14_msigdb_enrichment_all.R       MSigDB enrichment all 5 contrasts
    │   ├── 15_lowres_methylation_profile.R  Low-res browser tracks (chr1/chrX/chr5)
    │   ├── 16_tf_motif_enrichment.R         TF motif enrichment (monaLisa, negative result)
    │   ├── 16b_tf_motif_plots.R             TF motif result visualisation
    │   ├── 17_splice_junction_proximity.R   Splice junction proximity test (negative result)
    │   ├── submit_dmr_by_chr.sh             Submit per-chromosome DMR SLURM array
    │   ├── submit_smn1_masked_pipeline.sh   Submit masked alignment pipeline
    │   └── archive/                         Superseded scripts retained for reference
    ├── data/
    │   ├── reference/                       hg38 FASTA + Bismark genome index
    │   └── reference_smn1_masked/           SMN1-masked reference + Bismark index
    ├── results/                             All pipeline outputs (gitignored)
    └── logs/                                SLURM and tool logs (gitignored)

---

## Key Results

- 3,423 ASO-specific DMRs with bidirectional methylation changes (66% hyper, 34% hypo)
- 598,485 VPA DMRs -- near-exclusively hypomethylated, HDAC inhibitor signature
- 151 high-confidence ASO-specific DMRs converging on axonogenesis, ISL1, ROBO2
- chrX hotspot: 620 DMRs (18% of all ASO DMRs), disproportionate to chromosome size
- SMN2: 0 DMRs called (delta-meth 1.4%, below threshold) but H3K9me2 2.2-fold increase confirmed
- VPA and ASO act independently: ASO_VPA vs ASO_CTRL GO terms identical to VPA alone
- TF motif enrichment: negative (min padj=0.18, 746 motifs) -- no TF binding site disruption
- Splice junction proximity: negative (p=0.939) -- not genome-wide kinetic coupling
- MSigDB validation: neural pathway enrichment confirmed by second independent database

---

## Key References

- Marasco et al. (2022) Cell 185:2057-2070 -- nusinersen kinetic coupling + H3K9me2 at SMN2
- Catoni et al. (2018) Nucleic Acids Research 46:e114 -- DMRcaller
- Krueger and Andrews (2011) Bioinformatics 27:1571-1572 -- Bismark
- Yu et al. (2015) OMICS 19:284-287 -- ChIPseeker
- Wu et al. (2021) iMeta 1:e5 -- clusterProfiler 4.0
- Finkel et al. (2017) NEJM 377:1723-1732 -- ENDEAR trial nusinersen
- Gottlicher et al. (2001) EMBO J 20:6969-6978 -- VPA as HDAC inhibitor
- Brichta et al. (2003) Hum Mol Genet 12:2481-2489 -- VPA increases SMN2 expression
