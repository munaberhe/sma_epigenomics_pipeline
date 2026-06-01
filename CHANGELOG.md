# Changelog

## 2026-06-01
- Add DSS replicate-level validation (20_dss_replicate_testing.R)
- Add SMN2 H3K27ac enhancer analysis (19_smn2_h3k27ac_enhancer.R, GSE246399)
- Add sensitive local DMR calling at SMN2 (22_smn2_local_dmr_sensitive.R)
- Fix UpSet plot black rectangle
- Fix Fig4 H3K9me2 boxplot — load from saved RDS, proper two-panel figure
- Add chr13 hotspot annotation and intron 6/7 H3K27ac check

## 2026-05-31
- Add MSigDB enrichment all contrasts (14_msigdb_enrichment_all.R)
- Fix Fig5 negative results colour coding
- Regenerate thesis combined figures

## 2026-05-28
- Complete SMN1 masked realignment, all 12 samples
- Add SMN2 locus plots from masked chr5 CX reports
- Fix exon labelling: E2 split into E2a/E2b, E7 in red (Alberto Kornblihtt correction)
- Add Snakefile rules for all 5 contrasts

## 2026-05-15
- Lock DMR parameters: binSize=300, pVal=0.01, minDiff=0.20, minCyto=4
- Per-chromosome approach adopted after genome-wide attempts timed out
  (SLURM jobs 10504170, 10591016 -- 27h and 18h at 98GB, no completion)
- DMR calling complete for all 5 contrasts, 24 chromosomes

## 2026-05-05
- Initial meeting with Alberto Kornblihtt, Emilia Haberfeld, Marcos Miretti
- SMN1 masking approach agreed
- DMR parameters reviewed and confirmed with Radu

## 2026-04-01
- Initial WGBS QC and preliminary alignment
- Coverage: ~9x per replicate, 53.8% CpGs at >=10x pooled
- Bismark alignment rate ~70-75%

## 2026-03-10
- Project start, conda environment on Apocrita HPC
- hg38 GRCh38 Ensembl 109 reference downloaded
