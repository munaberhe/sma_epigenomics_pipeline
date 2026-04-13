.libPaths("~/R/library")
# dmrcaller.R
# Differentially Methylated Region Analysis
# SMA Epigenomics Pipeline — Muna Berhe, QMUL
# Per-chromosome loading strategy as recommended by Prof. Radu Zabet
# Run after all 12 methylation extractions are complete

library(DMRcaller)
library(BiocParallel)

# Config
COV_DIR  <- "results/alignments/bs"
OUT_DIR  <- "results/dmr"
MIN_COV  <- 10     # minimum coverage per CpG — confirm with Radu
CONTEXT  <- "CpG"  # CpG context only — confirm with Radu
CHROMOSOMES <- paste0("chr", c(1:22, "X", "Y"))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "per_chromosome"), recursive = TRUE, showWarnings = FALSE)

# Sample metadata
SAMPLES <- data.frame(
  name  = c("Scramble_CTRL_1", "Scramble_CTRL_2", "Scramble_CTRL_3",
            "Scramble_VPA_1",  "Scramble_VPA_2",  "Scramble_VPA_3",
            "ASO_CTRL_1",      "ASO_CTRL_2",      "ASO_CTRL_3",
            "ASO_VPA_1",       "ASO_VPA_2",       "ASO_VPA_3"),
  group = c(rep("Scramble_CTRL", 3), rep("Scramble_VPA", 3),
            rep("ASO_CTRL", 3),      rep("ASO_VPA", 3)),
  stringsAsFactors = FALSE
)

# Contrasts to run — confirm priorities with Radu
CONTRASTS <- list(
  ASO_vs_Scramble_CTRL     = list(treatment = "ASO_CTRL",     control = "Scramble_CTRL"),
  VPA_vs_Scramble_CTRL     = list(treatment = "Scramble_VPA", control = "Scramble_CTRL"),
  ASO_VPA_vs_Scramble_CTRL = list(treatment = "ASO_VPA",      control = "Scramble_CTRL"),
  ASO_VPA_vs_ASO_CTRL      = list(treatment = "ASO_VPA",      control = "ASO_CTRL"),
  ASO_VPA_vs_Scramble_VPA  = list(treatment = "ASO_VPA",      control = "Scramble_VPA")
)

load_cx_chromosome <- function(sample_name, chrom, min_cov = MIN_COV) {
  path <- file.path(COV_DIR,
    paste0(sample_name, "_bismark.deduplicated.bismark.cov.gz"))
  if (!file.exists(path)) {
    message("Missing: ", path)
    return(NULL)
  }
  tryCatch({
    readBismark(path,
                chromosome = chrom,
                context    = CONTEXT,
                coverage   = min_cov)
  }, error = function(e) {
    message("Error loading ", sample_name, " chr ", chrom, ": ", e$message)
    NULL
  })
}

run_dmr_chromosome <- function(contrast_name, contrast, chrom) {
  message("  Chromosome: ", chrom)
  treatment_samples <- SAMPLES$name[SAMPLES$group == contrast$treatment]
  control_samples   <- SAMPLES$name[SAMPLES$group == contrast$control]
  treatment_data <- lapply(treatment_samples, load_cx_chromosome, chrom = chrom)
  control_data   <- lapply(control_samples,   load_cx_chromosome, chrom = chrom)
  treatment_data <- Filter(Negate(is.null), treatment_data)
  control_data   <- Filter(Negate(is.null), control_data)
  if (length(treatment_data) < 2 || length(control_data) < 2) {
    message("  Skipping ", chrom, " — insufficient samples")
    return(NULL)
  }
  tryCatch({
    computeDMRs(treatment_data,
                control_data,
                regions        = NULL,
                method         = "bins",
                binSize        = 100,
                test           = "betareg",
                minCoverage    = MIN_COV,
                minProportionDiff = 0.1,
                pValueThreshold   = 0.05,
                cores          = 4)
  }, error = function(e) {
    message("  DMR error on ", chrom, ": ", e$message)
    NULL
  })
}

for (contrast_name in names(CONTRASTS)) {
  contrast <- CONTRASTS[[contrast_name]]
  message("\nRunning contrast: ", contrast_name,
          " (", contrast$treatment, " vs ", contrast$control, ")")
  contrast_out <- file.path(OUT_DIR, contrast_name)
  dir.create(contrast_out, recursive = TRUE, showWarnings = FALSE)
  all_dmrs <- list()
  for (chrom in CHROMOSOMES) {
    dmrs <- run_dmr_chromosome(contrast_name, contrast, chrom)
    if (!is.null(dmrs) && length(dmrs) > 0) {
      all_dmrs[[chrom]] <- dmrs
      saveRDS(dmrs, file.path(contrast_out,
                              paste0(contrast_name, "_", chrom, ".rds")))
    }
  }
  if (length(all_dmrs) > 0) {
    message("  Concatenating ", length(all_dmrs), " chromosomes...")
    all_dmrs_combined <- do.call(c, all_dmrs)
    saveRDS(all_dmrs_combined,
            file.path(contrast_out, paste0(contrast_name, "_all_chr.rds")))
    dmr_df <- as.data.frame(all_dmrs_combined)
    write.csv(dmr_df,
              file.path(contrast_out, paste0(contrast_name, "_DMRs.csv")),
              row.names = FALSE)
    message("  Total DMRs: ", length(all_dmrs_combined))
    message("  Hypermethylated: ",
            sum(dmr_df$methylationDifference > 0, na.rm = TRUE))
    message("  Hypomethylated: ",
            sum(dmr_df$methylationDifference < 0, na.rm = TRUE))
    chr_counts <- sapply(all_dmrs, length)
    chr_df <- data.frame(chromosome = names(chr_counts), n_dmrs = chr_counts)
    chr_df$chromosome <- factor(chr_df$chromosome,
                                levels = paste0("chr", c(1:22, "X", "Y")))
    library(ggplot2)
    p_chr <- ggplot(chr_df, aes(x = chromosome, y = n_dmrs)) +
      geom_col(fill = "#065A82") +
      labs(title    = paste0("DMRs per Chromosome — ", contrast_name),
           subtitle = paste0(contrast$treatment, " vs ", contrast$control),
           x = "Chromosome", y = "Number of DMRs") +
      theme_bw(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(contrast_out,
                     paste0(contrast_name, "_DMRs_per_chromosome.pdf")),
           p_chr, width = 12, height = 6)
    p_hist <- ggplot(dmr_df, aes(x = methylationDifference)) +
      geom_histogram(bins = 50, fill = "#1C7293", colour = "white") +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
      labs(title    = paste0("Methylation Difference — ", contrast_name),
           subtitle = paste0(contrast$treatment, " vs ", contrast$control),
           x = "Methylation difference", y = "Count") +
      theme_bw(base_size = 11)
    ggsave(file.path(contrast_out,
                     paste0(contrast_name, "_methylation_difference.pdf")),
           p_hist, width = 8, height = 6)
  } else {
    message("  No DMRs found for ", contrast_name)
  }
}

message("\nDone. Outputs saved to: ", OUT_DIR)
