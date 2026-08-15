#!/usr/bin/env Rscript
# 07f_pwmenrich_barplots.R
# Barplots of top PWMEnrich TF motif hits for all four contrasts.
# Reads MotifEnrichmentReport RDS objects directly -- no re-running motifEnrichment.
# Muna Berhe -- bt25018 -- QMUL MSc Bioinformatics

suppressPackageStartupMessages({
  library(PWMEnrich)
  library(ggplot2)
})
.libPaths(c("~/R/library", .libPaths()))

OUT_DIR <- "results/tf_motif/plots"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

# ---------------------------------------------------------------------------
# Report objects and display labels
# ---------------------------------------------------------------------------
REPORTS <- list(
  list(
    rds   = "results/tf_motif/pwmenrich_report.rds",
    name  = "ASO_specific_151dmrs",
    title = "ASO-specific DMRs (n=151)\nPWMEnrich top motifs",
    color = "#02C39A"
  ),
  list(
    rds   = "results/tf_motif/pwmenrich_report_ASO_VPA_5k.rds",
    name  = "ASO_VPA_vs_Scramble_CTRL_5k",
    title = "ASO+VPA vs Scramble_CTRL (5k DMR set)\nPWMEnrich top motifs",
    color = "#F59E0B"
  ),
  list(
    rds   = "results/tf_motif/pwmenrich_report_VPA_5k.rds",
    name  = "VPA_vs_Scramble_CTRL_5k",
    title = "VPA vs Scramble_CTRL (5k DMR set)\nPWMEnrich top motifs",
    color = "#1C7293"
  ),
  list(
    rds   = "results/tf_motif/pwmenrich_report_synergy_only_5k.rds",
    name  = "synergy_only_5k",
    title = "Synergy-only / combination-emergent DMRs\nPWMEnrich top motifs",
    color = "#D4A017"
  )
)

TOP_N <- 15

# ---------------------------------------------------------------------------
# Helper: extract top N rows and clean up target labels
# ---------------------------------------------------------------------------
extract_top <- function(report, n=TOP_N) {
  df <- report@d
  # p.value column name varies slightly -- find it
  pval_col <- grep("p.val|pval|p_val", names(df), value=TRUE, ignore.case=TRUE)[1]
  if (is.na(pval_col)) stop("Cannot find p-value column in report")
  df$pval <- df[[pval_col]]
  # Take top n by rank (already sorted)
  df <- head(df, n)
  # Clean target label: strip database prefix noise
  df$label <- gsub("Hsapiens-jolma2013-", "", df$target)
  df$label <- gsub("Hsapiens-", "", df$label)
  df$label <- gsub("_[0-9]+$", "", df$label)
  # log10 p-value for x axis
  df$neglog10p <- -log10(df$pval)
  # Preserve rank order for plotting (rank 1 at top)
  df$label <- factor(df$label, levels=rev(df$label))
  df
}

# ---------------------------------------------------------------------------
# Plot one report
# ---------------------------------------------------------------------------
plot_report <- function(cfg) {
  if (!file.exists(cfg$rds)) {
    message("Skipping -- file not found: ", cfg$rds)
    return(invisible(NULL))
  }
  message("Plotting: ", cfg$name)
  report <- readRDS(cfg$rds)
  df     <- tryCatch(extract_top(report),
                     error=function(e) { message("  ERROR: ", e$message); return(NULL) })
  if (is.null(df)) return(invisible(NULL))

  p <- ggplot(df, aes(x=neglog10p, y=label)) +
    geom_col(fill=cfg$color, width=0.7) +
    geom_text(aes(label=sprintf("p=%.1e", pval)),
              hjust=-0.05, size=3.0, colour="grey30") +
    scale_x_continuous(
      expand=expansion(mult=c(0, 0.25)),
      name="-log10(p-value)"
    ) +
    labs(
      title=cfg$title,
      y=NULL
    ) +
    theme_bw(base_size=11) +
    theme(
      plot.title   = element_text(size=11, face="bold", lineheight=1.2),
      axis.text.y  = element_text(size=9, face="bold"),
      axis.text.x  = element_text(size=8),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.margin  = margin(10, 20, 10, 10)
    )

  outfile <- file.path(OUT_DIR, paste0(cfg$name, "_barplot.pdf"))
  ggsave(outfile, p, width=7, height=5.5, device="pdf")
  message("  Saved: ", basename(outfile))
}

# ---------------------------------------------------------------------------
# Combined four-panel plot (for the presentation slide)
# ---------------------------------------------------------------------------
plot_combined <- function() {
  message("\nBuilding combined four-panel plot...")
  panels <- list()
  for (cfg in REPORTS) {
    if (!file.exists(cfg$rds)) next
    report <- readRDS(cfg$rds)
    df <- tryCatch(extract_top(report, n=10),
                   error=function(e) NULL)
    if (is.null(df)) next
    # Shorter title for panel
    short_title <- sub("\n.*", "", cfg$title)
    p <- ggplot(df, aes(x=neglog10p, y=label)) +
      geom_col(fill=cfg$color, width=0.7) +
      scale_x_continuous(expand=expansion(mult=c(0, 0.2))) +
      labs(title=short_title, x="-log10(p)", y=NULL) +
      theme_bw(base_size=9) +
      theme(
        plot.title         = element_text(size=8, face="bold"),
        axis.text.y        = element_text(size=7.5),
        axis.text.x        = element_text(size=7),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        plot.margin        = margin(6, 12, 6, 6)
      )
    panels[[length(panels)+1]] <- p
  }
  if (length(panels) == 0) { message("No panels built"); return(invisible(NULL)) }

  # Use patchwork if available, otherwise cowplot, otherwise gridExtra
  combined <- NULL
  if (requireNamespace("patchwork", quietly=TRUE)) {
    library(patchwork)
    combined <- Reduce("+", panels) + plot_layout(ncol=2)
  } else if (requireNamespace("cowplot", quietly=TRUE)) {
    library(cowplot)
    combined <- plot_grid(plotlist=panels, ncol=2)
  } else if (requireNamespace("gridExtra", quietly=TRUE)) {
    library(gridExtra)
    outfile <- file.path(OUT_DIR, "all_contrasts_combined.pdf")
    pdf(outfile, width=12, height=10)
    do.call(grid.arrange, c(panels, ncol=2))
    dev.off()
    message("  Saved: all_contrasts_combined.pdf")
    return(invisible(NULL))
  } else {
    message("Install patchwork, cowplot, or gridExtra for the combined plot")
    return(invisible(NULL))
  }

  outfile <- file.path(OUT_DIR, "all_contrasts_combined.pdf")
  ggsave(outfile, combined, width=12, height=10, device="pdf")
  message("  Saved: all_contrasts_combined.pdf")
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
for (cfg in REPORTS) plot_report(cfg)
plot_combined()

message("\nAll done. Outputs in: ", OUT_DIR)
