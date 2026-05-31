#!/usr/bin/env Rscript
# 07b_dmr_plots.R
# DMRcaller diagnostic and locus plots for the SMA epigenomics project

suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
  library(ggplot2)
})

.libPaths(c("~/R/library", .libPaths()))

DMR_DIR    <- "results/dmr"
BY_CHR_DIR <- "results/alignments/bs/by_chr"
OUT_DIR    <- "results/dmr/plots"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CHROMS <- paste0("chr", c(1:22, "X", "Y"))

# regionType in DMRcaller is counterintuitive:
#   "gain" means the treatment is hypomethylated (proportion1 < proportion2)
#   "loss" means the treatment is hypermethylated (proportion1 > proportion2)
# This was confirmed by checking chr13 DMRs where proportion1 (ASO_VPA) was
# consistently lower than proportion2 (Scramble_CTRL) yet regionType = "gain"

CONTRASTS <- list(
  list(
    name   = "ASO_CTRL_vs_Scramble_CTRL",
    cond_a = "ASO_CTRL",
    cond_b = "Scramble_CTRL",
    label  = "ASO alone vs Scramble CTRL\n(nusinersen effect, CTRL background)",
    colour = c(ASO_CTRL = "#E41A1C", Scramble_CTRL = "#377EB8")
  ),
  list(
    name   = "ASO_VPA_vs_Scramble_CTRL",
    cond_a = "ASO_VPA",
    cond_b = "Scramble_CTRL",
    label  = "ASO+VPA vs Scramble CTRL\n(combination vs baseline)",
    colour = c(ASO_VPA = "#984EA3", Scramble_CTRL = "#377EB8")
  ),
  list(
    name   = "ASO_VPA_vs_ASO_CTRL",
    cond_a = "ASO_VPA",
    cond_b = "ASO_CTRL",
    label  = "ASO+VPA vs ASO alone\n(VPA effect on ASO background)",
    colour = c(ASO_VPA = "#984EA3", ASO_CTRL = "#E41A1C")
  ),
  list(
    name   = "ASO_VPA_vs_Scramble_VPA",
    cond_a = "ASO_VPA",
    cond_b = "Scramble_VPA",
    label  = "ASO+VPA vs VPA alone\n(ASO effect on VPA background)",
    colour = c(ASO_VPA = "#984EA3", Scramble_VPA = "#FF7F00")
  ),
  list(
    name   = "Scramble_VPA_vs_Scramble_CTRL",
    cond_a = "Scramble_VPA",
    cond_b = "Scramble_CTRL",
    label  = "VPA alone vs Scramble CTRL\n(HDAC inhibitor effect)",
    colour = c(Scramble_VPA = "#FF7F00", Scramble_CTRL = "#377EB8")
  )
)

CONDITIONS <- list(
  ASO_CTRL      = paste0("ASO_CTRL_",      1:3),
  ASO_VPA       = paste0("ASO_VPA_",       1:3),
  Scramble_CTRL = paste0("Scramble_CTRL_", 1:3),
  Scramble_VPA  = paste0("Scramble_VPA_",  1:3)
)

# Loci to plot with DMR overlays.
# SMN1 coordinates are masked (chr5:70,924,941-70,953,015) so reads from
# SMN1 are not misassigned to SMN2. The SMN2 window here is wider than the
# gene body to capture promoter and downstream DMRs.
# RNA45SN2 coordinates from the top10 hypo annotation -- UPDATE these once
# you have confirmed the exact position from your annotated CSV.
KEY_LOCI <- list(
  list(
    name  = "SMN2",
    chr   = "chr5",
    start = 70040000,
    end   = 70090000,
    label = "SMN2 locus (masked alignment)"
  ),
  list(
    name  = "RNA45SN2_promoter",
    chr   = "chr21",
    start = 8158909,
    end   = 8259208,
    label = "RNA45SN2 promoter (top hit, p=1.56e-48, 109 CpGs, 41% drop)"
  ),
  list(
    name  = "MTA1-DT",
    chr   = "chr14",
    start = 105365123,
    end   = 105465422,
    label = "MTA1-DT locus (top ASO hypo hit, intron, p=1.72e-19)"
  ),
  list(
    name  = "MYO1D",
    chr   = "chr17",
    start = 32772454,
    end   = 32872753,
    label = "MYO1D locus (97 CpGs, p=7.83e-11)"
  )
)

# Step 1: Load pooled methylation data
# plotDMRs() needs the raw methylation GRanges objects, not just the DMR calls.
# Loading all 12 samples takes a while so I cache the pooled result to RDS.
# On the first run this takes ~20-30 min; after that it loads in seconds.

meth_cache <- file.path(DMR_DIR, "meth_pooled_cache.rds")

if (file.exists(meth_cache)) {
  message("Loading cached pooled methylation data...")
  meth_pooled <- readRDS(meth_cache)
  message("  Conditions loaded: ", paste(names(meth_pooled), collapse = ", "))
} else {
  message("No cache found -- pooling from CpG reports (this takes ~20-30 min)...")

  meth_pooled <- list()

  for (cond in names(CONDITIONS)) {
    message("  Pooling: ", cond)
    reps <- CONDITIONS[[cond]]

    # Load each replicate chromosome by chromosome then concatenate.
    # readBismark() reads the bismark cytosine report format which includes
    # zero-coverage positions -- these are needed for correct DMRcaller stats.
    glist <- GRangesList(lapply(reps, function(s) {
      chr_files <- file.path(BY_CHR_DIR,
        paste0(s, "_", CHROMS, ".CpG_report.txt.gz"))
      chr_files <- chr_files[file.exists(chr_files)]
      if (length(chr_files) == 0) stop("No chr files found for sample: ", s)
      do.call(c, lapply(chr_files, readBismark))
    }))

    # poolMethylationDatasets sums read counts across the three replicates
    # at each CpG position, giving ~27x effective depth per condition.
    # We pool rather than use computeDMRsReplicates because per-replicate
    # coverage (~9x) is too low to reliably estimate between-replicate variance.
    meth_pooled[[cond]] <- poolMethylationDatasets(glist)
    message("    CpGs after pooling: ",
            format(length(meth_pooled[[cond]]), big.mark = ","))
  }

  saveRDS(meth_pooled, meth_cache)
  message("  Cached to: ", meth_cache)
}

# Step 2: Load the DMR GRanges objects
# These were produced by 06c_dmrcaller_combine_chr.R and are the locked
# parameter results from May 23. Nothing needs to be recomputed here.

message("\nLoading DMR results...")
dmr_results <- list()
for (ct in CONTRASTS) {
  rds_path <- file.path(DMR_DIR, paste0("dmr_", ct$name, ".rds"))
  if (!file.exists(rds_path)) {
    message("  Missing: ", ct$name, " -- skipping")
    next
  }
  dmr_results[[ct$name]] <- readRDS(rds_path)
  message("  ", ct$name, ": ", length(dmr_results[[ct$name]]), " DMRs")
}

# Step 3: Per-chromosome DMR counts
# Mirrored bar chart with hypo below the axis and hyper above, one bar per
# chromosome. Useful for spotting if a particular chromosome is driving the
# signal or if the effect is uniformly distributed (expected for VPA).

message("\nGenerating per-chromosome DMR plots...")

for (ct in CONTRASTS) {
  if (is.null(dmr_results[[ct$name]])) next
  dmrs <- dmr_results[[ct$name]]

  chr_counts <- as.data.frame(table(
    chr  = as.character(seqnames(dmrs)),
    type = dmrs$regionType
  ))
  chr_counts$chr <- factor(chr_counts$chr,
    levels = paste0("chr", c(1:22, "X", "Y")))
  chr_counts$direction <- ifelse(chr_counts$type == "gain", "Hypo", "Hyper")

  # Make hypo negative so it mirrors below zero on the y axis
  chr_counts$Freq_signed <- ifelse(
    chr_counts$direction == "Hypo",
    -chr_counts$Freq,
     chr_counts$Freq
  )

  p <- ggplot(chr_counts, aes(x = chr, y = Freq_signed, fill = direction)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey30") +
    scale_fill_manual(
      values = c(Hypo = "#4393C3", Hyper = "#D6604D"),
      name   = "Direction"
    ) +
    scale_y_continuous(
      labels = function(x) format(abs(x), big.mark = ",")
    ) +
    labs(
      title    = paste("DMRs per chromosome --", ct$name),
      subtitle = ct$label,
      x        = "Chromosome",
      y        = "Number of DMRs (hypo below zero, hyper above)",
      caption  = paste0(
        "bins 300bp | minDiff 0.20 | p 0.01 | minCpGs 4 | ",
        "minGap 300bp | score test | pooled replicates\n",
        "Total: ", format(length(dmrs), big.mark = ","),
        "  Hypo: ", format(sum(dmrs$regionType == "gain"), big.mark = ","),
        "  Hyper: ", format(sum(dmrs$regionType == "loss"), big.mark = ",")
      )
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
      plot.title       = element_text(face = "bold"),
      plot.caption     = element_text(size = 7, colour = "grey40"),
      legend.position  = "top",
      panel.grid.minor = element_blank()
    )

  out_path <- file.path(OUT_DIR,
    paste0(ct$name, "_DMRs_per_chromosome.pdf"))
  ggsave(out_path, p, width = 12, height = 5)
  message("  Saved: ", out_path)
}

# Step 4: Methylation difference distributions
# Histogram of proportion1 - proportion2 across all DMRs.
# Negative values = hypomethylated in the treatment condition.
# The dotted lines mark the +/-20% threshold so you can see how far
# beyond the minimum threshold most DMRs actually fall.

message("\nGenerating methylation difference plots...")

for (ct in CONTRASTS) {
  if (is.null(dmr_results[[ct$name]])) next
  dmrs <- dmr_results[[ct$name]]

  df <- data.frame(
    methDiff  = dmrs$proportion1 - dmrs$proportion2,
    direction = ifelse(dmrs$regionType == "gain", "Hypo", "Hyper")
  )

  p <- ggplot(df, aes(x = methDiff, fill = direction)) +
    geom_histogram(bins = 80, colour = "white", linewidth = 0.15) +
    geom_vline(xintercept = 0, linewidth = 0.6, linetype = "dashed",
               colour = "grey20") +
    geom_vline(xintercept = c(-0.20, 0.20), linewidth = 0.4,
               linetype = "dotted", colour = "grey50") +
    scale_fill_manual(
      values = c(Hypo = "#4393C3", Hyper = "#D6604D"),
      name   = "Direction"
    ) +
    scale_x_continuous(
      limits = c(-1, 1),
      breaks = seq(-1, 1, 0.2),
      labels = scales::percent_format(accuracy = 1)
    ) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      title    = paste("Methylation difference distribution --", ct$name),
      subtitle = ct$label,
      x        = "Methylation difference (treatment minus reference)",
      y        = "Number of DMRs",
      caption  = paste0(
        "Dotted lines mark the +/-20% minimum threshold | n = ",
        format(nrow(df), big.mark = ","), " DMRs"
      )
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold"),
      plot.caption     = element_text(size = 7, colour = "grey40"),
      legend.position  = "top",
      panel.grid.minor = element_blank()
    )

  out_path <- file.path(OUT_DIR,
    paste0(ct$name, "_methylation_difference.pdf"))
  ggsave(out_path, p, width = 8, height = 5)
  message("  Saved: ", out_path)
}

# Step 5: DMR locus overlay plots
# This is what was missing from the earlier SMN locus plots -- the shaded
# rectangles showing where DMRcaller actually called significant regions.
# plotDMRs() draws the methylation proportion for both conditions and then
# overlays the DMR calls as coloured boxes below the trace (blue = hypo,
# red = hyper). This is the standard way to present DMRcaller results.

message("\nGenerating DMR locus overlay plots...")

for (ct in CONTRASTS) {
  if (is.null(dmr_results[[ct$name]])) next
  dmrs   <- dmr_results[[ct$name]]
  meth_a <- meth_pooled[[ct$cond_a]]
  meth_b <- meth_pooled[[ct$cond_b]]

  for (locus in KEY_LOCI) {

    # Pull out only the DMRs that fall within this window so plotDMRs
    # doesn't try to render the entire genome
    locus_gr   <- GRanges(locus$chr, IRanges(locus$start, locus$end))
    dmrs_locus <- subsetByOverlaps(dmrs, locus_gr)
    n_dmrs     <- length(dmrs_locus)

    out_path <- file.path(OUT_DIR, paste0(
      ct$name, "_", locus$name, "_locus_overlay.pdf"
    ))

    pdf(out_path, width = 10, height = 5, bg="white")
    par(bg="white", col.axis="black", col.lab="black", col.main="black", fg="black")

    tryCatch({
      region_gr <- GRanges(locus$chr, IRanges(locus$start, locus$end))
      plotLocalMethylationProfile(
        methylationData1 = meth_a,
        methylationData2 = meth_b,
        region           = region_gr,
        DMRs             = NULL,
        conditionsNames  = c(ct$cond_a, ct$cond_b),
        windowSize       = 300,
        context          = "CG",
        main             = paste0(
          locus$name, " -- ", ct$name,
          "\n", locus$label,
          "\n(", n_dmrs, " DMRs in window)"
        )
      )
    }, error = function(e) {
      message("    plotLocalMethylationProfile error for ", locus$name,
              " / ", ct$name, ": ", conditionMessage(e))
      plot.new()
      text(0.5, 0.5,
           paste0("plotLocalMethylationProfile failed\n\n",
                  conditionMessage(e),
                  "\n\nLocus: ", locus$name,
                  "\nContrast: ", ct$name),
           cex = 0.9, col = "red")
    })
    dev.off()
    message("  Saved: ", out_path, "  (", n_dmrs, " DMRs in window)")
  }
}

message("\nAll plots saved to: ", OUT_DIR)
plots <- list.files(OUT_DIR, pattern = "\\.pdf$", full.names = FALSE)
for (f in sort(plots)) message("  ", f)
