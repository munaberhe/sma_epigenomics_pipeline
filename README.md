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
│   ├── regional/           # Regional/chromosome-level DMR density suite
│   └── utils/               # Shared colour palette
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
| `21_smn2_predup_dmr.R` | DMR calling on pre-dedup masked data for coverage-gap comparison |
| `22_smn2_predup_locus_plots.R` | SMN2 locus plots using pre-dedup masked alignments |

### DMR Calling (pipeline/)

| Script | Description |
|--------|-------------|
| `02_dmr_calling.R` | Per-chromosome DMR calling via DMRcaller (called by dmrcaller_by_chr.R) |
| `dmrcaller_by_chr.R` | DMRcaller SLURM array worker (one job per contrast per chromosome) |
| `03_combine_tested_windows.R` | Combine per-chromosome tested windows for annotation background |
| `04_check_outputs.sh` | Verify all expected output files exist |

### QC and Benchmarking (pipeline/)

| Script | Description |
|--------|-------------|
| `01_qc.R` | PCA, M-bias, duplication rates, bisulfite conversion, coverage curves |
| `09_benchmark_plots.R` | DMRcaller parameter selection: label-swap and read-count permutation nulls |
| `fig_correlation_heatmap.R` | Per-sample genome-wide methylation correlation heatmap |
| `fig_lowres_chr1.R` | chr1 low-resolution methylation profile with centromere annotation |
| `fig_qc_panel.R` | Assembles the multi-panel QC figure |
| `24_mingap_sensitivity.R` | DMR count sensitivity to the minGap threshold, genome-wide |
| `23_aso_vpa2_sensitivity.R` | Recomputes all four contrasts with the ASO_VPA_2 replicate excluded |

### Annotation and Enrichment (pipeline/)

| Script | Description |
|--------|-------------|
| `02_dmr_annotate.R` | ChIPseeker annotation + GO/KEGG enrichment for all four contrasts |
| `06_tss_heatmap.R` | TSS methylation heatmap |
| `07_grant_fig1cd_annotation_obsexp.R` | Obs/exp genomic annotation enrichment + permutation test |
| `10_dmr_significance_plots.R` | DMR count summaries and significance plots |
| `11_volcano_plots.R` | Methylation difference volcano plots per contrast |
| `12_circos_ideogram.R` | Circos ideogram and diverging bar chromosome plots |
| `13_mds_dmrsize.R` | MDS from pooled methylation correlation + DMR size distributions |
| `14_upset_pairwise.R` | UpSet plot of DMR overlap across four pairwise contrasts |
| `15_gokegg_pairwise.R` | GO/KEGG pathway enrichment panels for all four contrasts |
| `21_enhancer_enrichment_pairwise.R` | H9 enhancer overlap by positional permutation, all four contrasts |
| `22_enhancer_plot.R` | Fold-enrichment summary plot from the enhancer permutation results |
| `25_metagene_profile.R` | Scaled TSS-to-TES metagene methylation profile per condition |

### Regional and Chromosome-Level Analysis (pipeline/regional/)

| Script | Description |
|--------|-------------|
| `00_regional_common.R` | Shared configuration, coordinate reference and helper functions |
| `40_regional_hotspot_scan.R` | Coverage-weighted regional DMR density scan, genome-wide and per-chromosome |
| `41_smn2_adjacent_interaction.R` | Raw-count 2x2 test of the SMN2 flanking region and relaxed-scan positions |
| `42_chrX_dmr_composition.R` | chrX representation, directional composition and positioning among DMRs |

### Candidate Gene Analysis (pipeline/)

| Script | Description |
|--------|-------------|
| `17_pairwise_context_scan.R` | ASO/VPA context-dependent DMR scan and synergy scoring |
| `18_relevance_scoring.R` | Candidate gene relevance scoring and bubble chart |
| `20_master_locus_plots.R` | All candidate gene locus plots across the four pairwise contrasts |
| `genuine_synergy_scan.R` | Additive-null permutation screen (legacy methodology, documented for provenance) |
| `smn2/07_smn2_sensitive.R` | Sensitive SMN2 DMR scan (minDiff=5%, minCyto=3) |
| `16_smn2_extended_igv.R` | SMN2 extended locus: methylation, DMRs, H3K27ac and regulatory elements |

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

### minGap Sensitivity

DMR counts vary with the minGap threshold as expected, but the rank order and direction of the four pairwise contrasts are stable across the tested range (100 to 1000 bp). Reported alongside the locked parameter choice as a robustness check rather than a re-optimisation.

---

## Enhancer Enrichment

DMRs in all four pairwise contrasts overlap H9 predicted non-promoter enhancers more often than a coverage-matched genomic null (regioneR, randomizeRegions, 1000 permutations, seed 42, 50,000-region subsample for contrasts exceeding that count). Fold enrichment ranges from 1.63x to 2.33x depending on contrast, with the empirical p-value at the permutation floor in every case. H9 is a human embryonic stem cell line; the enhancer set is a cross-cell-type positional prior for the HEK293T background used here, not a validated enhancer annotation for these cells.

---

## Regional and Chromosome-Level Findings

A genome-wide scan for DMR density hotspots, weighted by covered-CpG background and tested against 1,000 genomic randomisations per contrast, was run for the first time under the pairwise framework (not previously performed under the legacy combined-contrast design).

**chrX** shows chromosome-level DMR enrichment specific to the ASO-alone contrast (observed/expected fold enrichment computed against covered-CpG share), with the two VPA-referenced contrasts instead depleted relative to the same background. Directional composition of chrX DMRs in the ASO-alone contrast does not differ from the autosomal composition, consistent with a variance-driven rather than a directional methylation effect. Baseline methylation on chrX is substantially below the autosomal mean, plausibly limiting the detectability of a fixed absolute-difference threshold on that chromosome; this accounts for part but not all of the observed depletion pattern.

**chr13** shows regional DMR density elevation confined to a sub-chromosomal interval in the VPA-referenced contrasts, but no genome-wide chromosome-level excess. Baseline methylation in that interval is elevated relative to the rest of the genome, consistent with the same detection-threshold mechanism operating in the opposite direction from chrX.

---

## SMN2 Locus

WGBS coverage at the SMN2 locus (chr5:70,049,638-70,078,522) is low relative to genome-wide mean depth, reflecting multi-mapping between the SMN1/SMN2 paralogs. A flanking region approximately 9.7 kb downstream of the annotated SMN2 3' end (chr5:70,088,223-70,088,522) shows a substantial pooled methylation loss under the combined ASO+VPA condition relative to each single-treatment arm (raw counts and interaction estimate reported in the thesis text); permutation testing at this locus has limited resolution given the small number of covered cytosines and replicates, and the finding is reported descriptively alongside its confidence limits rather than as a validated significant interaction. No DMR is called within the SMN2 gene body itself under the locked genome-wide parameters in any contrast.

---

## Non-Additive Survivor Loci

Seven loci originally flagged under an additive-null permutation screen (predicted combination effect = single-treatment effects summed; deviation from that prediction tested against a permutation null) were re-examined for consistency with the pairwise contrast framework. Three of the seven (RELL2, DDIT4L, TCEAL4) show a called DMR in the pairwise contrast matching their originally claimed context-dependent category and no called DMR in the corresponding single-treatment comparator. The remaining four either show no called DMR at this locus under the locked parameters, or are called in both the defining contrast and its comparator and therefore do not meet the context-dependent criterion as originally stated.

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
| Permutation enrichment | regioneR | v1.42.0 |
| Visualisation | ggplot2, patchwork, UpSetR, ComplexHeatmap | -- |
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
  "rtracklayer", "UpSetR", "patchwork", "ggplot2", "karyoploteR",
  "regioneR", "ComplexHeatmap", "aod"
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

# 7. Regional and chromosome-level scan
Rscript pipeline/regional/40_regional_hotspot_scan.R
Rscript pipeline/regional/41_smn2_adjacent_interaction.R
Rscript pipeline/regional/42_chrX_dmr_composition.R
```

---

## References

- Marasco et al. (2022) Cell 185:2057-2070 -- nusinersen, H3K9me2, kinetic coupling
- Stigliano et al. (2025) bioRxiv 10.1101/2025.10.10.681673 -- ASO-promoted gene looping at SMN2
- Catoni et al. (2018) Nucleic Acids Research 46:e114 -- DMRcaller
- Krueger and Andrews (2011) Bioinformatics 27:1571-1572 -- Bismark
- Yu et al. (2015) OMICS 19:284-287 -- ChIPseeker
- Wu et al. (2021) iMeta 1:e5 -- clusterProfiler 4.0
- Gel and Malinverni (2019) Bioinformatics 35:1064-1065 -- regioneR
- Finkel et al. (2017) NEJM 377:1723-1732 -- nusinersen clinical trial
- Gottlicher et al. (2001) EMBO J 20:6969-6978 -- VPA as HDAC inhibitor

# Last updated: August 2026
