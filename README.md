# SMA Epigenomics Pipeline

Whole-genome bisulfite sequencing (WGBS) analysis of nusinersen (ASO) and valproic acid (VPA) combination treatment in Spinal Muscular Atrophy cell lines.

**MSc Bioinformatics thesis · Queen Mary University of London · 2026**  
**Supervisor:** Prof Radu Zabet (Zabet Lab, QMUL)  
**Collaborators:** Prof Alberto Kornblihtt (IFIBYNE-UBA-CONICET), Dr Emilia Haberfeld, Dr Marcos Miretti

---

## Study Design

12 HEK293T samples in a 2x2 factorial design (3 replicates per condition):

| Condition | ASO | VPA |
|-----------|-----|-----|
| ASO_CTRL | Yes | No |
| Scramble_CTRL | No | No |
| ASO_VPA | Yes | Yes |
| Scramble_VPA | No | Yes |

Four pairwise contrasts isolating each treatment effect:

| Contrast | Biological question |
|----------|-------------------|
| ASO_CTRL vs Scramble_CTRL | ASO effect alone |
| Scramble_VPA vs Scramble_CTRL | VPA effect alone |
| ASO_VPA vs Scramble_VPA | ASO effect on VPA background |
| ASO_VPA vs ASO_CTRL | VPA effect on ASO background |

---

## Repository Structure

```
sma_epigenomics_pipeline/
├── pipeline/
│   ├── canonical/          # Trimming, alignment, dedup, extraction (unmasked)
│   ├── smn1_masked/        # SMN2 locus analysis using SMN1-masked reference
│   ├── smn2_predup/        # SMN2 coverage gap investigation (no deduplication)
│   └── utils/              # Shared colour palette
├── smn2/                   # SMN2-specific sensitive DMR scan
├── thesis_figures/         # Thesis figure generation
├── configs/                # SLURM profile and pipeline config
└── environment.yml         # Conda environment
```

---

## Pipeline Scripts

### Alignment — Canonical Unmasked (pipeline/canonical/)

| Script | Description |
|--------|-------------|
| `00_genome_prep.sh` | Download hg38 GRCh38, build Bismark bisulfite index |
| `01_trim_array.sh` | Trim Galore paired-end array job (12 samples, 6 at a time) |
| `02_align_array.sh` | Bismark paired-end alignment array job (12 samples, 3 at a time) |
| `03_dedup_extract_array.sh` | Deduplicate BAMs, extract methylation, split CX report by chromosome |

### Alignment — SMN1-Masked (pipeline/smn1_masked/)

| Script | Description |
|--------|-------------|
| `align_00_smn_reference.sh` | Download SMN1/SMN2 sequences and build masked reference |
| `align_01_mask_index.sh` | Build Bismark index on SMN1-masked hg38 |
| `align_02_bismark.sh` | Align all 12 samples to SMN1-masked reference |
| `align_03_dedup_extract.sh` | Deduplicate and extract methylation from masked BAMs |
| `align_04_split_chr.sh` | Split CX reports by chromosome, retain CpG context only |

### SMN2 Coverage Gap Investigation (pipeline/smn2_predup/)

| Script | Description |
|--------|-------------|
| `03_extract_array.sh` | Extract methylation without deduplication to assess coverage gap |

### Downstream Analysis (pipeline/)

| Script | Description |
|--------|-------------|
| `utils/00_sma_palette.R` | Canonical colour palette (single source of truth for all figures) |
| `01_qc.R` | PCA, M-bias, duplication rates, bisulfite conversion, coverage curves |
| `02_dmr_calling.R` | Per-chromosome DMR calling via DMRcaller (called by dmrcaller_by_chr.R) |
| `02_dmr_annotate.R` | ChIPseeker annotation + GO/KEGG enrichment for all four contrasts |
| `03_combine_tested_windows.R` | Combine per-chromosome tested windows for annotation background |
| `04_check_outputs.sh` | Verify all expected output files exist |
| `dmrcaller_by_chr.R` | DMRcaller SLURM array worker (one job per contrast per chromosome) |
| `06_grant_fig1a_sd_density.R` | SD density plot of methylation variability across conditions |
| `06_tss_heatmap.R` | TSS methylation heatmap |
| `07_grant_fig1cd_annotation_obsexp.R` | Obs/exp genomic annotation enrichment + 1000 permutation test |
| `08_grant_fig2a_motif_venn_logos.R` | TF motif Venn diagram and sequence logos |
| `10_dmr_significance_plots.R` | DMR count summaries and significance plots |
| `11_volcano_plots.R` | Methylation difference volcano plots per contrast |
| `12_circos_ideogram.R` | Circos ideogram and diverging bar chromosome plots |
| `13_mds_dmrsize.R` | MDS from pooled methylation correlation + DMR size distributions |
| `14_upset_pairwise.R` | UpSet plot of DMR overlap across four pairwise contrasts |
| `15_gokegg_pairwise.R` | GO/KEGG pathway enrichment panels for all four contrasts |
| `16_smn2_extended_igv.R` | SMN2 extended locus: methylation, DMRs, H3K27ac and regulatory elements |
| `17_pairwise_context_scan.R` | ASO/VPA context-dependent DMR scan and synergy scoring |
| `18_relevance_scoring.R` | Candidate gene relevance scoring and bubble chart |
| `20_master_locus_plots.R` | All candidate gene locus plots (18 genes, 4 contrasts each) |
| `21_smn2_predup_dmr.R` | DMR calling on pre-dedup masked data for SMN2 coverage analysis |
| `22_smn2_predup_locus_plots.R` | SMN2 locus plots using pre-dedup masked alignments |
| `smn2/07_smn2_sensitive.R` | Sensitive SMN2 DMR scan (minDiff=5%, minCyto=3) |

---

## DMR Calling Parameters

Parameters benchmarked on chr1 using two null models before locking.

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `method` | bins | Fixed-width windows, robust at pooled coverage |
| `binSize` | 300 bp | Largest signal-to-noise gap vs permutation null |
| `minProportionDifference` | 0.20 | Removes biologically marginal changes |
| `pValueThreshold` | 0.01 | Standard genome-wide threshold |
| `minCytosinesCount` | 4 | Prevents single-CpG noise calls |
| `minReadsPerCytosine` | 4 | Minimum coverage per cytosine |
| `minGap` | 300 bp | Prevents iterative merge on dense VPA contrasts |
| `test` | score | Rao score test at ~27x pooled coverage |
| `context` | CG | CHG/CHH near-zero in human somatic cells |

Note: regionType = gain means hypomethylated; loss means hypermethylated.

Per-chromosome calling was necessary after genome-wide attempts timed out at 98 GB RAM after 18-27 hours.

### Benchmarking Null Models

**Null 1 - read-count permutation:** shuffles readsM and readsN by a random position index applied identically to both conditions. Preserves coverage structure but destroys spatial methylation signal. 20 permutations per window size.

**Null 2 - label swap:** swaps ASO_VPA and ASO_CTRL conditions entirely. Most conservative - preserves spatial structure, coverage and methylation levels. Tests whether the method detects directionality rather than noise.

The optimal binSize maximises the signal-to-noise ratio z = (D_obs - mu_null) / sigma_null where D_obs is the observed DMR count and mu_null and sigma_null are the mean and SD across permutations. Both null models agreed on 300 bp as optimal.

---

## Tech Stack

| Component | Tool | Version |
|-----------|------|---------|
| Trimming | Trim Galore | v0.6.11 |
| Alignment | Bismark | v0.25.1 |
| Aligner | Bowtie2 | v2.5.4 |
| Reference genome | hg38 GRCh38 Ensembl 109 | SMN1-masked |
| Deduplication | deduplicate_bismark | v0.25.1 |
| DMR calling | DMRcaller (Bioconductor) | v1.42.0 |
| Annotation | ChIPseeker | v1.46.1 |
| Pathway enrichment | clusterProfiler | v4.18.4 |
| Visualisation | ggplot2, patchwork, UpSetR | -- |
| HPC scheduler | SLURM | v24.11.7 |
| Language | R | v4.5.1 |
| OS | Rocky Linux 9 | Apocrita HPC QMUL |

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
  "rtracklayer", "UpSetR", "patchwork", "ggplot2", "karyoploteR"
))
```

---

## Running the Pipeline

```bash
# 1. Prepare genome
sbatch pipeline/canonical/00_genome_prep.sh

# 2. Trim
sbatch pipeline/canonical/01_trim_array.sh

# 3. Align
sbatch pipeline/canonical/02_align_array.sh

# 4. Deduplicate and extract methylation
sbatch pipeline/canonical/03_dedup_extract_array.sh

# 5. Call DMRs per chromosome then combine
sbatch --array=1-24 pipeline/dmrcaller_by_chr.R
Rscript pipeline/03_combine_tested_windows.R

# 6. Downstream analysis
Rscript pipeline/02_dmr_annotate.R
Rscript pipeline/15_gokegg_pairwise.R
Rscript pipeline/20_master_locus_plots.R
```

---

## SMN2 Coverage Gap

WGBS reads at the SMN2 locus (chr5:70,049,638-70,078,522) show consistently low coverage (mean 0.13x, 6.2% of CpGs covered at 1x) across all conditions and even after SMN1 masking. Three independent analyses confirmed this is not an artefact of PCR deduplication. The gap reflects multi-mapping between SMN1/SMN2 paralogs combined with genuine low sequencing depth at this locus. No DMRs were detected at SMN2 under any contrast.

---

## References

- Marasco et al. (2022) Cell 185:2057-2070 -- nusinersen, H3K9me2, kinetic coupling
- Catoni et al. (2018) Nucleic Acids Research 46:e114 -- DMRcaller
- Krueger and Andrews (2011) Bioinformatics 27:1571-1572 -- Bismark
- Yu et al. (2015) OMICS 19:284-287 -- ChIPseeker
- Wu et al. (2021) iMeta 1:e5 -- clusterProfiler 4.0
- Finkel et al. (2017) NEJM 377:1723-1732 -- nusinersen clinical trial
- Gottlicher et al. (2001) EMBO J 20:6969-6978 -- VPA as HDAC inhibitor
# Last updated: August 2026
