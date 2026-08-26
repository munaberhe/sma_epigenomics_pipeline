# Changelog

## 2026-06-01 - DSS validation and enhancer analysis

- Add DSS replicate-level DMR validation (20_dss_replicate_testing.R)
  688 high-confidence DMRs at n=3 per group confirm neural pathway enrichment
- Add SMN2 H3K27ac enhancer analysis (19_smn2_h3k27ac_enhancer.R)
  Cross-referenced GSE246399 (Calandrelli et al.); no peaks at SMN2 3 prime end
- Add intron 6/7 specific H3K27ac check and chr13 hotspot annotation
  (21_intron67_enhancer_chr13_annotation.R)
- Add sensitive local DMR calling at SMN2 (22_smn2_local_dmr_sensitive.R)
  minDiff=0.05, binSize=100bp focused on chr5 SMN2 region
- Fix UpSet plot black rectangle (mb.ratio and text_scale adjustments)
- Fix Fig4 H3K9me2 boxplot to use saved h3k9me2_plot_df.rds

## 2026-05-31 - Thesis figures finalised

- Add MSigDB enrichment across all 5 contrasts (14_msigdb_enrichment_all.R)
  Independently confirms neural pathway signal using C2/C5/C3 gene set collections
- Fix Fig5 negative results colour coding (hardcoded hex values)
- Regenerate thesis combined figures (18_thesis_combined_figures.R)
- Add UpSet plot for ASO-specific DMR set (12_upset_dmr_intersections.R)
  151 ASO-specific DMRs not present in any VPA contrast

## 2026-05-28 - SMN1 masking complete

- Complete SMN1 masked realignment for all 12 samples
  Reference: chr5:70,924,941-70,953,015 replaced with Ns in hg38
- Add SMN2 locus methylation plots from masked chr5 CX reports
  plotLocalMethylationProfile, WIN_SIZE=300bp, FLANK=2000bp
- Fix exon labelling following Alberto Kornblihtt correction
  E2 split into E2a/E2b; penultimate exon correctly labelled E7 in red
- Update Snakefile to include all 5 contrasts and 8 new downstream rules

## 2026-05-20 - DMR annotation and pathway enrichment

- Add ChIPseeker annotation for all 5 contrasts (07_dmr_annotate.R)
  Promoter defined as +/-2kb from TSS; separate hyper/hypo GO enrichment
- Add per-chromosome DMR bar charts and methylation diff histograms (07b_dmr_plots.R)
- Add locus overlay plots for top hits (07c_dmr_locus_plots.R)
- Add top 10 hypo DMRs per contrast (09_top10_dmrs.R)
- Identify RNA45SN2 (chr21) and MTA1-DT as top ASO off-target loci
- Add H3K9me2 overlap analysis (11_h3k9me2_overlap.R)
  Validates kinetic coupling model using GSE167762 (Marasco et al. 2022)
- Add TF motif enrichment (16_tf_motif_enrichment.R) - result: min p.adj=0.18
- Add splice junction proximity test (17_splice_junction_proximity.R) - result: p=0.939

## 2026-05-15 - DMR calling complete

- Lock DMR parameters after benchmarking on chr1/chr6/chr13 permutation null:
  binSize=300bp, pValueThreshold=0.01, minProportionDifference=0.20,
  minCytosinesCount=4, minReadsPerCytosine=4, minGap=300bp, test=score
- Per-chromosome approach adopted after genome-wide attempts timed out
  SLURM jobs 10504170 and 10591016 ran 27h and 18h at 98GB without completing
- DMR calling complete across all 5 contrasts and 24 chromosomes (120 SLURM jobs)
  ASO_CTRL: 3,423 DMRs; Scramble_VPA: 598,485 DMRs; ASO_VPA: 554,291 DMRs

## 2026-05-10 - QC complete

- 12-sample PCA from per-replicate chr1 CpG methylation (08_pca_chr1.R)
  4 groups cluster cleanly; within-group correlations 0.845-0.93
- M-bias, duplication rates, bisulfite conversion QC (08b_additional_qc.R)
  Conversion rate 96-96.8% consistent across all 12 samples
- Coverage retention curves (10_coverage_qc.R)
  53.8% CpGs at >=10x after pooling; mean depth 23.5x

## 2026-05-05 - Meeting with Kornblihtt lab

- Meeting with Prof Alberto Kornblihtt, Dr Emilia Haberfeld, Dr Marcos Miretti
- SMN1 masking approach agreed (vs merge alternative)
- DMR parameters reviewed and confirmed with Radu Zabet
- Exon labelling correction noted by Alberto

## 2026-04-15 - Alignment complete

- Bismark alignment complete for all 12 samples to unmasked hg38
  Alignment rate ~70-75%; CpG methylation ~69-70% (expected for HEK293T)
- Deduplication and methylation extraction complete
- CX reports split by chromosome for per-chromosome DMR calling
- Per-replicate coverage ~9x; replicates pooled to ~27x for DMR calling

## 2026-04-01 - Initial QC and preliminary analysis

- Raw FASTQ QC with FastQC; adapter trimming with Trim Galore v0.6.11
  Q20 quality threshold; >99% adapter removal
- Preliminary Bismark alignment to hg38 GRCh38 Ensembl 109
- Coverage analysis: mean ~9x per replicate
- Preliminary methylation profiles show VPA global hypomethylation vs no
  difference between ASO_CTRL and Scramble_CTRL at chromosome scale

## 2026-03-10 - Project start

- Repository initialised, conda environment set up on Apocrita HPC
  (Rocky Linux 9, SLURM 24.11.7)
- hg38 GRCh38 Ensembl 109 reference genome downloaded
- Bismark genome index built
- Initial Snakefile covering QC, trimming and alignment stages
- configs/config.yaml with pipeline parameters
