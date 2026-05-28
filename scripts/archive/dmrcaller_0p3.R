
.libPaths("~/R/library")
# dmrcaller.R
# Differentially Methylated Region Analysis
# SMA Epigenomics Pipeline — Muna Berhe, QMUL
# Parameters follow DMRcaller vignette (Zabet et al.) sections 3.2-3.9

library(DMRcaller)
library(ggplot2)
library(GenomicRanges)

# Config
COV_DIR     <- "results/alignments/bs/by_chr"
OUT_DIR     <- "results/dmr_0p3"
QC_DIR      <- "results/qc/dmrcaller"
MIN_COV     <- 4
CONTEXT     <- "CG"
CHROMOSOMES <- paste0("chr", c(1:22, "X"))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR,  recursive = TRUE, showWarnings = FALSE)

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

# Contrasts
CONTRASTS <- list(
  ASO_vs_Scramble_CTRL     = list(treatment = "ASO_CTRL",     control = "Scramble_CTRL"),
  VPA_vs_Scramble_CTRL     = list(treatment = "Scramble_VPA", control = "Scramble_CTRL"),
  ASO_VPA_vs_Scramble_CTRL = list(treatment = "ASO_VPA",      control = "Scramble_CTRL"),
  ASO_VPA_vs_ASO_CTRL      = list(treatment = "ASO_VPA",      control = "ASO_CTRL"),
  ASO_VPA_vs_Scramble_VPA  = list(treatment = "ASO_VPA",      control = "Scramble_VPA")
)

# Load and pool replicates for one group on one chromosome
load_group_chromosome <- function(sample_names, chrom, min_cov = MIN_COV) {
  paths <- file.path(COV_DIR,
    paste0(sample_names, "_", chrom, ".CpG_report.txt.gz"))
  missing <- !file.exists(paths)
  if (any(missing)) {
    message("    Missing files: ", paste(paths[missing], collapse=", "))
    paths <- paths[!missing]
  }
  if (length(paths) < 2) {
    message("    Insufficient samples for ", chrom)
    return(NULL)
  }
  tryCatch({
    dat <- readBismarkPool(paths)
    dat <- dat[dat$readsN >= min_cov]
    if (length(dat) == 0) return(NULL)
    dat
  }, error = function(e) {
    message("    Error loading ", chrom, ": ", e$message)
    NULL
  })
}

# Section 3.2 — Low resolution methylation profiles
# Run once for all groups on chr1 as QC overview
message("\n=== Section 3.2: Low resolution profiles (chr1) ===")
tryCatch({
  ref_sample_ctrl <- SAMPLES$name[SAMPLES$group == "Scramble_CTRL"]
  ref_sample_vpa  <- SAMPLES$name[SAMPLES$group == "ASO_VPA"]
  ctrl_chr1 <- load_group_chromosome(ref_sample_ctrl, "chr1")
  vpa_chr1  <- load_group_chromosome(ref_sample_vpa,  "chr1")
  if (!is.null(ctrl_chr1) && !is.null(vpa_chr1)) {
    regions_chr1 <- GRanges(seqnames = Rle("chr1"),
                            ranges   = IRanges(1, 248956422))
    profileCtrl <- computeMethylationProfile(ctrl_chr1, regions_chr1,
                                             windowSize = 10000, context = "CG")
    profileVPA  <- computeMethylationProfile(vpa_chr1,  regions_chr1,
                                             windowSize = 10000, context = "CG")
    profileList <- GRangesList("Scramble_CTRL" = profileCtrl,
                               "ASO_VPA"       = profileVPA)
    pdf(file.path(QC_DIR, "low_resolution_profile_chr1.pdf"), width=14, height=5)
    par(mar=c(4,4,3,1)+0.1)
    plotMethylationProfile(profileList,
                           autoscale = FALSE,
                           labels    = NULL,
                           title     = "CG methylation — chr1 (10kb windows)",
                           col       = c("#065A82","#F59E0B"),
                           pch       = c(1, 0),
                           lty       = c(4, 1))
    dev.off()
    message("  Saved: low_resolution_profile_chr1.pdf")
  }
}, error = function(e) message("  Low res profile error: ", e$message))

# Main DMR calling loop
for (contrast_name in names(CONTRASTS)) {
  contrast <- CONTRASTS[[contrast_name]]
  message("\n=== Running contrast: ", contrast_name,
          " (", contrast$treatment, " vs ", contrast$control, ") ===")

  contrast_out <- file.path(OUT_DIR, contrast_name)
  dir.create(contrast_out, recursive = TRUE, showWarnings = FALSE)

  all_dmrs <- list()

  for (chrom in CHROMOSOMES) {
    message("  Chromosome: ", chrom)

    treatment_samples <- SAMPLES$name[SAMPLES$group == contrast$treatment]
    control_samples   <- SAMPLES$name[SAMPLES$group == contrast$control]

    treatment_data <- load_group_chromosome(treatment_samples, chrom)
    control_data   <- load_group_chromosome(control_samples,   chrom)

    if (is.null(treatment_data) || is.null(control_data)) {
      message("    Skipping ", chrom, " — insufficient data")
      next
    }

    # Explicit GRanges region as per vignette section 3.5
    # Use chromosome length that covers all data
    chrom_max <- max(end(treatment_data), end(control_data))
    regions <- GRanges(seqnames = Rle(chrom),
                       ranges   = IRanges(1, chrom_max))

    dmrs <- tryCatch({
      # Section 3.5 — computeDMRs with vignette parameters
      # noise_filter uses minGap=0 (vignette default for this method)
      computeDMRs(treatment_data,
                  control_data,
                  regions                 = regions,
                  context                 = CONTEXT,
                  method                  = "noise_filter",
                  windowSize              = 100,
                  kernelFunction          = "triangular",
                  test                    = "score",
                  pValueThreshold         = 0.01,
                  minCytosinesCount       = 4,
                  minProportionDifference = 0.3,
                  minGap                  = 0,
                  minSize                 = 50,
                  minReadsPerCytosine     = 4,
                  parallel                = TRUE,
                  cores                   = 4)
    }, error = function(e) {
      message("    DMR error on ", chrom, ": ", e$message)
      NULL
    })

    if (is.null(dmrs) || length(dmrs) == 0) {
      message("    No DMRs on ", chrom)
      next
    }


    saveRDS(dmrs, file.path(contrast_out,
                            paste0(contrast_name, "_", chrom, ".rds")))
    all_dmrs[[chrom]] <- dmrs
  }

  if (length(all_dmrs) == 0) {
    message("  No DMRs found for ", contrast_name)
    next
  }

  # Combine all chromosomes
  message("  Concatenating ", length(all_dmrs), " chromosomes...")
  all_dmrs_combined <- do.call(c, unname(all_dmrs))
  saveRDS(all_dmrs_combined,
          file.path(contrast_out, paste0(contrast_name, "_all_chr.rds")))

  dmr_df <- as.data.frame(all_dmrs_combined)
  write.csv(dmr_df,
            file.path(contrast_out, paste0(contrast_name, "_DMRs.csv")),
            row.names = FALSE)

  message("  Total DMRs: ",         length(all_dmrs_combined))
  message("  Hypermethylated: ",    sum(dmr_df$regionType == "gain", na.rm=TRUE))
  message("  Hypomethylated: ",     sum(dmr_df$regionType == "loss", na.rm=TRUE))

  # Chromosome distribution plot
  chr_counts <- sapply(all_dmrs, length)
  chr_df <- data.frame(chromosome = names(chr_counts), n_dmrs = chr_counts)
  chr_df$chromosome <- factor(chr_df$chromosome,
                              levels = paste0("chr", c(1:22, "X")))
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

  # Methylation difference histogram
  p_hist <- ggplot(dmr_df, aes(x = proportion1 - proportion2)) +
    geom_histogram(bins = 50, fill = "#1C7293", colour = "white") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
    labs(title    = paste0("Methylation Difference — ", contrast_name),
         subtitle = paste0(contrast$treatment, " vs ", contrast$control),
         x = "Methylation difference (treatment - control)", y = "Count") +
    theme_bw(base_size = 11)
  ggsave(file.path(contrast_out,
                   paste0(contrast_name, "_methylation_difference.pdf")),
         p_hist, width = 8, height = 6)

  # Section 3.9 — plotLocalMethylationProfile for top 3 DMRs
  message("  Plotting local methylation profiles for top DMRs...")
  top_dmrs <- head(all_dmrs_combined[order(all_dmrs_combined$pValue)], 3)

  for (i in seq_along(top_dmrs)) {
    dmr_i   <- top_dmrs[i]
    chrom_i <- as.character(seqnames(dmr_i))
    start_i <- max(1, start(dmr_i) - 5000)
    end_i   <- end(dmr_i) + 5000
    region_i <- GRanges(seqnames = Rle(chrom_i),
                        ranges   = IRanges(start_i, end_i))

    tryCatch({
      treatment_samples <- SAMPLES$name[SAMPLES$group == contrast$treatment]
      control_samples   <- SAMPLES$name[SAMPLES$group == contrast$control]
      t_dat <- load_group_chromosome(treatment_samples, chrom_i)
      c_dat <- load_group_chromosome(control_samples,   chrom_i)

      if (!is.null(t_dat) && !is.null(c_dat)) {
        pdf(file.path(contrast_out,
                      paste0(contrast_name, "_DMR", i, "_locus_",
                             chrom_i, "_", start(dmr_i), ".pdf")),
            width = 14, height = 6)
        par(mar=c(4,4,3,1)+0.1)
        plotLocalMethylationProfile(
          t_dat, c_dat,
          region_i,
          GRangesList(DMRs = all_dmrs_combined),
          conditionsNames = c(contrast$treatment, contrast$control),
          windowSize = 300,
          main = paste0("CG methylation — ", contrast_name,
                        " — DMR ", i, " (", chrom_i, ":",
                        start(dmr_i), "-", end(dmr_i), ")"))
        dev.off()
        message("    Saved locus plot ", i)
      }
    }, error = function(e) message("    Locus plot error: ", e$message))
  }
}

message("\nDone. Outputs saved to: ", OUT_DIR)
