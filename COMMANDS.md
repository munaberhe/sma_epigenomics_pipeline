# SMA Epigenomics Pipeline — Command Record
# Muna Berhe, bt25018, QMUL Zabet Lab
# Thesis: Genome-Wide Pleiotropic Epigenetic Effects of Combined ASO1 and VPA Treatment in SMA

## 0. Environment setup
module load miniforge/24.7.1
conda activate sma_epigenomics_pipeline
module unload miniforge/24.7.1
module load R/4.5.1

## 1. Reference preparation
sbatch scripts/align_00_smn_reference.sh

## 2. Alignment
sbatch scripts/align_01_mask_index.sh
sbatch scripts/align_02_bismark.sh
sbatch scripts/align_03_dedup_extract.sh
sbatch scripts/align_04_split_chr.sh
sbatch scripts/align_05_cleanup.sh
sbatch scripts/align_06_smn_merged.sh

## 3. QC
sbatch --wrap="Rscript scripts/01_qc.R"
sbatch --wrap="Rscript scripts/coverage_correlation_qc.R"

## 4. DMR calling
# Genome-wide (all chromosomes)
sbatch scripts/dmr_00_submit_masked.sh
sbatch scripts/dmr_01_submit_by_chr.sh
# Combine per-chromosome results
sbatch --wrap="Rscript scripts/03_dmr_combine.R"

## 5. DMR annotation and enrichment
sbatch --mem=64G --time=3:00:00 --partition=compute \
  --wrap="Rscript scripts/02_dmr_annotate.R"
sbatch --mem=32G --time=1:00:00 --partition=compute \
  --wrap="Rscript scripts/03_enrichment.R"

## 6. TF motif enrichment (PWMEnrich)
sbatch --mem=64G --time=8:00:00 --cpus-per-task=4 --partition=compute \
  --wrap="Rscript scripts/04_tf_motif.R"

## 7. Locus plots
sbatch --mem=64G --time=2:00:00 --cpus-per-task=2 --partition=compute \
  --wrap="Rscript scripts/05_locus_plots.R"
sbatch --mem=32G --time=2:00:00 --partition=compute \
  --wrap="Rscript scripts/05_smn2_locus_final.R"

## 8. SMN2 sensitive DMR analysis
sbatch --mem=32G --time=1:00:00 --partition=compute \
  --wrap="Rscript scripts/07_smn2_sensitive.R"

## 9. TSS heatmap
sbatch --mem=32G --time=1:00:00 --partition=compute \
  --wrap="Rscript scripts/06_tss_heatmap.R"

## 10. Manhattan plot
sbatch --mem=32G --time=30:00 --partition=compute \
  --wrap="Rscript scripts/09_manhattan.R"

## 11. DSS validation
sbatch --mem=64G --time=2:00:00 --partition=compute \
  --wrap="Rscript scripts/08_dss_validation.R"

## 12. Bock-style annotation (CpG islands + genomic features)
sbatch --mem=64G --time=2:00:00 --cpus-per-task=2 --partition=compute \
  --wrap="Rscript scripts/bock_style_annotation.R"

## 13. Hotspot screen
sbatch --mem=16G --time=15:00 --partition=compute \
  --wrap="Rscript scripts/lowres_hotspot_screen.R"

## 14. Parameter benchmark
# Label swap (array)
sbatch --array=1-24 --mem=32G --time=6:00:00 --cpus-per-task=16 --partition=compute \
  --wrap="Rscript scripts/benchmark_labelswap_array.R"
# Read count permutation (array)
sbatch --array=1-6 --mem=64G --time=12:00:00 --cpus-per-task=32 --partition=compute \
  --wrap="Rscript scripts/benchmark_readcount_array.R"
# MinGap sweep (array)
sbatch --array=1-12 --mem=64G --time=36:00:00 --cpus-per-task=32 --partition=compute \
  --wrap="Rscript scripts/benchmark_mingap_array.R"
# Combine results and plot
Rscript scripts/benchmark_combine_results.R
Rscript scripts/plot_radu_panels.R

## 15. Sync results to Mac
bash ~/Desktop/sync_sma_results.sh
