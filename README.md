# SMA Epigenomics Pipeline

MSc Bioinformatics thesis project — Muna Berhe, Queen Mary University of London, 2026

Supervisor: Professor Radu Zabet

---

## What this is

This pipeline was built for my thesis on assessing genome-wide pleiotropic effects of a combined ASO1+VPA treatment in Spinal Muscular Atrophy. The biological question is whether valproic acid — a broad HDAC inhibitor used alongside an antisense oligonucleotide targeting SMN2 exon 7 — has off-target transcriptional and epigenetic effects across the genome.

The pipeline runs on the Apocrita HPC cluster at QMUL and covers the full analysis from raw FASTQ files through to differential expression, differential methylation, and integration of both.

## Background

SMA is caused by loss-of-function mutations in SMN1. The SMN2 paralog compensates partially but skips exon 7 in most transcripts, producing very little functional protein. Nusinersen corrects this splicing but Marasco et al. (2022, Cell) showed it also deposits the repressive mark H3K9me2 at the SMN2 locus. VPA counteracts this chromatin compaction and the combined ASO1+VPA treatment outperforms ASO1 alone in patient fibroblasts and mouse models. This project looks at what else VPA is doing genome-wide.

## Pipeline structure

Raw FASTQ files go through QC (FastQC/MultiQC) and trimming (Trim Galore), then split into two branches — STAR alignment for RNA-seq and Bismark alignment for bisulfite sequencing. RNA-seq reads are counted with featureCounts and fed into DESeq2 for differential expression. Bisulfite alignments go through deduplication, methylation extraction, and DMRcaller for differential methylation. The integration step uses GenomicRanges to find overlapping DEGs and DMRs, followed by GO/KEGG pathway enrichment with clusterProfiler.

Everything is managed by Snakemake and runs on SLURM.

## Repository layout
```
configs/        pipeline parameters (config.yaml)
data/           raw FASTQ files, processed files, reference genome and indices
scripts/        R analysis scripts and SLURM submission scripts
results/        QC reports, BAM files, count matrices, differential results, figures
logs/           SLURM and tool logs
docs/           notes and documentation
Snakefile       main workflow definition
environment.yml conda environment
```

## Setup

Clone the repo and set up the conda environment:
```bash
git clone https://github.com/munaberhe/sma_epigenomics_pipeline.git
cd sma_epigenomics_pipeline
conda env create -f environment.yml
conda activate sma_epigenomics_pipeline
```

Install R packages:
```r
install.packages("BiocManager")
BiocManager::install(c(
  "DESeq2", "DMRcaller", "bsseq",
  "clusterProfiler", "org.Hs.eg.db",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "ChIPseeker", "EnhancedVolcano",
  "pheatmap", "GenomicRanges", "rtracklayer"
))
install.packages(c("tidyverse", "ggplot2", "RColorBrewer"))
```

Download the reference genome and build indices:
```bash
sbatch scripts/download_hg38.sh
sbatch scripts/build_star_index.sh
sbatch scripts/build_bismark_index.sh
```

Update `configs/config.yaml` with your sample paths before running.

## Running the pipeline

Dry run first to check everything resolves:
```bash
snakemake --dry-run --cores 1
```

Submit to SLURM:
```bash
snakemake --executor slurm --default-resources mem_mb=8000 runtime=120 --jobs 10
```

## Reference genome

hg38 (GRCh38) from UCSC, annotated with hg38.ensGene.gtf.

## Key references

- Marasco et al. (2022) Cell 185:2057-2070
- Catoni et al. (2018) Nucleic Acids Research 46:e114
- Love et al. (2014) Genome Biology 15:550
