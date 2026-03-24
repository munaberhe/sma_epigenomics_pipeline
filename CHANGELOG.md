# Changelog

All notable changes to the SMA Epigenomics Pipeline are documented here.

---

## [0.3.0] - 2026-03-24

### Added
- Full R analysis scripts: deseq2.R, dmrcaller.R, integrate.R
- Snakemake rules for DESeq2, DMRcaller, featureCounts, Bismark extraction and DMR-DEG integration
- Test alignment scripts for STAR and Bismark validation
- README.md with full installation, usage and pipeline documentation

---

## [0.2.0] - 2026-03-10

### Added
- hg38 reference genome download script
- STAR genome index build script
- Bismark bisulfite genome index build script
- Test data download and subsetting scripts

---

## [0.1.0] - 2026-03-10

### Added
- Initial project structure and folder layout
- Conda environment: sma_epigenomics_pipeline (Python 3.11)
- Full bioinformatics tool stack installed on Apocrita HPC
- Snakefile covering QC, trimming and alignment stages
- configs/config.yaml with pipeline parameters
- .gitignore configured for large reference and data files
