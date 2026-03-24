# SMA Epigenomics Pipeline

**Assessing Genome-Wide Pleiotropic Effects of a Combined RNA Antisense Oligonucleotide Treatment in Spinal Muscular Atrophy**

MSc Bioinformatics Thesis — Muna Berhe | Queen Mary University of London | 2026

---

## Background

Spinal Muscular Atrophy (SMA) is a severe autosomal recessive motor neuron disease caused by loss-of-function mutations in the *SMN1* gene. The *SMN2* paralog provides a partial backup but produces only ~10% functional SMN protein due to constitutive skipping of exon 7.

Nusinersen (Spinraza) is an FDA-approved antisense oligonucleotide (ASO) therapy that corrects SMN2 exon 7 splicing, restoring SMN protein levels. However, Marasco et al. (2022, *Cell*) showed that the ASO also deposits the repressive histone mark H3K9me2 at the SMN2 locus, creating a roadblock to RNA Polymerase II elongation that partially counteracts its own splicing correction.

Valproic acid (VPA), an HDAC inhibitor, counteracts this chromatin compaction and cooperates with the ASO to enhance exon 7 inclusion. The ASO1+VPA combination outperforms ASO1 alone in SMA patient fibroblasts and mouse models.

**The open question this project addresses:** VPA is a genome-wide chromatin modifier. What are its pleiotropic transcriptional and epigenetic off-target effects across the genome, and what are their potential biological consequences?

---

## Research Aims

| Aim | Analysis | Tool |
|-----|----------|------|
| 1 | Differential gene expression: ASO1+VPA vs ASO1 | DESeq2 |
| 2 | Alternative splicing changes, including SMN2 exon 7 | rMATS |
| 3 | Differential DNA methylation (DMRs) integrated with expression | DMRcaller |

---

## Pipeline Overview
```
data/raw/*.fastq.gz
        |
        v
[ 1. QC & Trimming ]  FastQC -> MultiQC -> Trim Galore
        |
   -----+------
   |           |
   v           v
[ STAR ]   [ Bismark ]    2. Alignment (hg38)
 RNA-seq    BS-seq
   |           |
   v           v
[featureCounts] [Bismark extractor]    3. Quantification
   |           |
   v           v
[ DESeq2 ]  [ DMRcaller ]    4. Differential Analysis
  DEGs        DMRs
   |           |
   +-----+-----+
         |
         v
  [ Integration ]    5. GenomicRanges findOverlaps()
  DMR-DEG overlap
         |
         v
  [ Enrichment ]    6. clusterProfiler GO/KEGG
```

---

## Repository Structure
```
sma_epigenomics_pipeline/
├── configs/
│   └── config.yaml
├── data/
│   ├── raw/
│   ├── processed/
│   └── reference/
├── results/
│   ├── qc/
│   ├── alignments/
│   │   ├── rna/
│   │   └── bs/
│   ├── counts/
│   ├── differential/
│   └── figures/
├── scripts/
│   ├── deseq2.R
│   ├── dmrcaller.R
│   ├── integrate.R
│   ├── build_star_index.sh
│   ├── build_bismark_index.sh
│   └── download_hg38.sh
├── logs/
├── docs/
├── Snakefile
└── README.md
```

---

## Dependencies

### Bioinformatics Tools

| Tool | Version | Purpose |
|------|---------|---------|
| FastQC | 0.12.1 | Read quality control |
| Trim Galore | — | Adapter trimming |
| STAR | 2.7.11b | RNA-seq alignment |
| Bismark | 0.25.1 | Bisulfite sequencing alignment |
| samtools | 1.22.1 | BAM file processing |
| featureCounts | — | RNA-seq read counting |
| MultiQC | — | QC report aggregation |
| Snakemake | 9.17.2 | Workflow management |

### R / Bioconductor Packages

| Package | Purpose |
|---------|---------|
| DESeq2 | Differential expression analysis |
| DMRcaller | Differential methylation analysis |
| bsseq | Bisulfite sequencing data structures |
| clusterProfiler | GO and KEGG pathway enrichment |
| GenomicRanges | Genomic interval operations |
| rtracklayer | Genomic file import/export |
| EnhancedVolcano | Volcano plot visualisation |
| pheatmap | Heatmap visualisation |
| ChIPseeker | Genomic annotation of DMRs |
| tidyverse | Data manipulation and plotting |

---

## Installation & Setup

### 1. Clone the repository
```bash
git clone https://github.com/YOUR_USERNAME/sma_epigenomics_pipeline.git
cd sma_epigenomics_pipeline
```

### 2. Set up conda environment
```bash
conda create -n sma_epigenomics_pipeline python=3.11 -y
conda activate sma_epigenomics_pipeline
conda install -c bioconda -c conda-forge \
    fastqc trim-galore bismark star samtools subread multiqc bowtie2 snakemake -y
pip install snakemake-executor-plugin-slurm
```

### 3. Install R packages
```r
install.packages("BiocManager")
BiocManager::install(c(
  "DESeq2", "DMRcaller", "bsseq", "DSS",
  "clusterProfiler", "org.Hs.eg.db",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "annotatr", "ChIPseeker",
  "EnhancedVolcano", "pheatmap", "ComplexHeatmap"
))
install.packages(c("tidyverse", "ggplot2", "RColorBrewer"))
```

### 4. Download reference genome and build indices
```bash
sbatch scripts/download_hg38.sh
sbatch scripts/build_star_index.sh
sbatch scripts/build_bismark_index.sh
```

### 5. Configure the pipeline

Edit `configs/config.yaml` with your sample names and file paths.

---

## Running the Pipeline

### Dry run
```bash
snakemake --dry-run
```

### Run locally
```bash
snakemake --cores 8
```

### Run on HPC (SLURM)
```bash
snakemake \
  --executor slurm \
  --default-resources mem_mb=8000 runtime=120 \
  --jobs 10
```

---

## Key Outputs

| File | Description |
|------|-------------|
| results/qc/multiqc_report.html | Combined QC report |
| results/differential/deseq2_results.csv | DEG table with log2FC and padj |
| results/differential/dmrs.csv | DMR table with methylation differences |
| results/differential/dmr_deg_overlap.csv | DMR-DEG overlap table |
| results/figures/volcano_plot.png | Volcano plot |
| results/figures/pca_plot.png | PCA of RNA-seq samples |
| results/figures/heatmap_top50.png | Top 50 DEG heatmap |
| results/figures/dmr_deg_scatter.png | Methylation vs expression scatter |
| results/figures/dmr_annotation_pie.png | Genomic distribution of DMRs |

---

## Reference Genome

- Genome: hg38 (GRCh38) from UCSC
- Annotation: hg38.ensGene.gtf from UCSC

---

## Citation

If you use this pipeline please cite:

- Marasco et al. (2022) Cell 185:2057-2070
- Catoni et al. (2018) Nucleic Acids Research 46:e114
- Love et al. (2014) Genome Biology 15:550

---

## Author

Muna Berhe
MSc Bioinformatics, Queen Mary University of London
Supervisor: Professor Radu Zabet
Thesis submission: August 21, 2026

---

## Acknowledgements

This research utilised Queen Mary's Apocrita HPC facility, supported by QMUL Research-IT.
http://doi.org/10.5281/zenodo.438045
