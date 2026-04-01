# SMA Epigenomics Pipeline
MSc Bioinformatics thesis project — Muna Berhe, Queen Mary University of London, 2026
Supervisor: Professor Radu Zabet

---

## What this is

This pipeline was built for my thesis on assessing genome-wide pleiotropic epigenetic effects of a combined ASO1+VPA treatment in Spinal Muscular Atrophy. The biological question is whether valproic acid — a broad HDAC inhibitor used alongside an antisense oligonucleotide targeting SMN2 exon 7 — introduces off-target DNA methylation changes across the genome.

The pipeline runs on the Apocrita HPC cluster at QMUL and covers the full analysis from raw FASTQ files through to differential methylation analysis and pathway enrichment.

## Background

SMA is caused by loss-of-function mutations in SMN1. The SMN2 paralog compensates partially but skips exon 7 in most transcripts, producing very little functional protein. Nusinersen corrects this splicing but Marasco et al. (2022, Cell) showed it also deposits the repressive mark H3K9me2 at the SMN2 locus. VPA counteracts this chromatin compaction and the combined ASO1+VPA treatment outperforms ASO1 alone in patient fibroblasts and mouse models. This project characterises the genome-wide DNA methylation landscape of the combined treatment, identifying which genomic regions and biological pathways are epigenetically affected.

## Experimental design

Four conditions, three biological replicates each (12 samples total):

| Condition | Treatment | Purpose |
|---|---|---|
| ASO_CTRL | ASO1 only | ASO1 effect on methylation |
| ASO_VPA | ASO1 + VPA | Combined treatment effect |
| Scramble_CTRL | Scrambled ASO only | Negative control |
| Scramble_VPA | Scrambled ASO + VPA | VPA effect alone |

Sequencing: paired-end whole genome bisulfite sequencing (WGBS), 151bp, ~350M reads per sample, reference genome hg38.

## Pipeline structure

Raw FASTQ files go through QC (FastQC/MultiQC) and trimming (Trim Galore), then bisulfite alignment with Bismark. Aligned reads are deduplicated and methylation is extracted across all three cytosine contexts (CpG, CHG, CHH). Differential methylation regions are called with DMRcaller, annotated with ChIPseeker, and pathway enrichment is performed with clusterProfiler.

Everything is managed by Snakemake and runs on SLURM. All resource definitions (memory, runtime) are embedded in the Snakefile — no external flags needed.

## Current pipeline status

| Stage | Status |
|---|---|
| FastQC (raw reads) | Complete — all 24 files |
| MultiQC | Complete |
| Trim Galore | Complete — all 12 samples |
| FastQC (trimmed reads) | Complete — all 24 files |
| Bismark alignment | In progress — running 1 sample at a time |
| Bismark deduplication | Pending |
| Methylation extraction | Pending |
| DMRcaller | Pending — parameters TBD with supervisor |
| ChIPseeker annotation | Pending |
| clusterProfiler enrichment | Pending |

## Key QC findings

- Read quality: Phred 38-40 across all samples — excellent
- GC content: ~21-23% — expected for bisulfite sequencing due to C→T conversion
- Sequence length: uniform 151bp across all samples
- Duplication: 13-31% — within expected range for WGBS
- Adapter content: detected in raw reads, removed by Trim Galore

## Repository layout
```
configs/        pipeline parameters (config.yaml)
data/           raw FASTQ files, processed files, reference genome and indices
scripts/        R analysis scripts and SLURM submission scripts
results/        QC reports, BAM files, methylation data, differential results, figures
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
  "DMRcaller", "bsseq",
  "clusterProfiler", "org.Hs.eg.db",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "ChIPseeker", "GenomicRanges", "rtracklayer"
))
install.packages(c("tidyverse", "ggplot2", "RColorBrewer"))
```

Download the reference genome and build the Bismark index:
```bash
sbatch scripts/download_hg38.sh
sbatch scripts/build_bismark_index.sh
```

Update `configs/config.yaml` with your sample paths before running.

## Running the pipeline

Set this before running Snakemake on the login node to prevent thread exhaustion:
```bash
export OPENBLAS_NUM_THREADS=1
```

Dry run first to check everything resolves:
```bash
snakemake --dry-run --cores 1
```

Submit to SLURM:
```bash
snakemake --unlock

nohup snakemake \
  --executor slurm \
  --jobs 1 \
  --latency-wait 60 \
  --keep-going \
  --rerun-incomplete \
  > logs/snakemake_run.log 2>&1 &
```

Note: `--jobs 1` limits bismark to one sample at a time. Each sample generates ~240GB of intermediate temp files during alignment — running more than one simultaneously risks exceeding the Apocrita scratch quota. All resource definitions (memory, runtime) are embedded in the Snakefile.

## Bismark alignment parameters

Based on lab scripts provided by supervisor:
```
--bowtie2 -N 1 -L 20 --score_min L,0,-0.6 --parallel 4
```

- `-N 1` — allows 1 mismatch in seed
- `-L 20` — shorter seed length for more sensitive alignment
- `--score_min L,0,-0.6` — relaxed scoring to improve mapping efficiency

## Reference genome

hg38 (GRCh38), bisulfite index built with Bismark genome preparation.

## Key references

- Marasco et al. (2022) Cell 185:2057-2070
- Catoni et al. (2018) Nucleic Acids Research 46:e114
- Krueger & Andrews (2011) Bioinformatics 27:1571-1572
