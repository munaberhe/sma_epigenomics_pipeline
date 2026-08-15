.libPaths(c("~/R/library", .libPaths()))
library(pdftools)

setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")
out_dir <- "results/thesis_figures"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

pdfs <- c(
  "results/dmr_benchmark_final_plots/labelswap_real_corrected.pdf",
  "results/dmr_benchmark_final_plots/readcount_real_corrected.pdf",
  "results/dmr_benchmark_final_plots/readcount_real_coverage.pdf"
)

pdf_combine(pdfs, output="results/thesis_figures/Fig4.3_DMR_benchmark_panel.pdf")
message("Saved: results/thesis_figures/Fig4.3_DMR_benchmark_panel.pdf")
