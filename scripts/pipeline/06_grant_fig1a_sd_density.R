#!/usr/bin/env Rscript
.libPaths("~/R/library")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(bsseq)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# SD density plot of per-CpG methylation, Grant et al 2026 Fig 1A style.
# three overlays: CpGs in DMRs, CpGs in tested-but-not-DMR windows, all tested CpGs.
# input : results/bsseq/bsseq_filtered.rds (BSseq object, 10x in >=2 reps per group)
#         results/dmr/dmr_<contrast>.rds
# output: results/plots/grant_fig1a_sd_density_<contrast>.{pdf,png}

OUT_DIR <- "results/plots"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

CONTRASTS_V <- c(
  "ASO_VPA_vs_Scramble_CTRL",
  "Scramble_VPA_vs_Scramble_CTRL",
  "ASO_CTRL_vs_Scramble_CTRL",
  "ASO_VPA_vs_ASO_CTRL",
  "ASO_VPA_vs_Scramble_VPA"
)

# load filtered BSseq once
bs_path <- "results/bsseq/bsseq_filtered.rds"
if (!file.exists(bs_path)) stop("missing: ", bs_path)
bs <- readRDS(bs_path)

# per-CpG SD of methylation beta across all 12 samples
M <- getCoverage(bs, type="M")
Cov <- getCoverage(bs, type="Cov")
beta <- M / pmax(Cov, 1)  # avoid div by 0
sd_all <- apply(beta, 1, sd, na.rm=TRUE)
cpg_gr <- granges(bs)

for (contrast in CONTRASTS_V) {
  rds_path <- file.path("results/dmr", paste0("dmr_", contrast, ".rds"))
  if (!file.exists(rds_path)) { message("missing: ", rds_path); next }
  dmrs <- readRDS(rds_path)
  if (length(dmrs) == 0) { message("no DMRs for: ", contrast); next }

  # split CpGs into 3 groups
  in_dmr <- overlapsAny(cpg_gr, dmrs)
  # tested-but-not-DMR == all CpGs in tested windows that are not inside a DMR
  # if a "tested windows" file exists use it, otherwise approximate by all CpGs
  tested_path <- file.path("results/dmr", paste0("tested_windows_", contrast, ".rds"))
  if (file.exists(tested_path)) {
    tested <- readRDS(tested_path)
    in_tested <- overlapsAny(cpg_gr, tested)
  } else {
    in_tested <- rep(TRUE, length(cpg_gr))
  }

  df <- bind_rows(
    data.frame(sd = sd_all[in_dmr],                  group = "DMR CpGs"),
    data.frame(sd = sd_all[in_tested & !in_dmr],     group = "non-DMR tested CpGs"),
    data.frame(sd = sd_all,                          group = "all tested CpGs")
  ) %>% filter(is.finite(sd), sd > 0)

  # log10 for x like Grant
  df$log10sd <- log10(df$sd)

  p <- ggplot(df, aes(x=log10sd, fill=group, colour=group)) +
    geom_density(alpha=0.35, linewidth=0.6, adjust=1.1) +
    scale_fill_manual(values = c(
      "DMR CpGs"            = "#20808D",
      "non-DMR tested CpGs" = "#B2182B",
      "all tested CpGs"     = "#7A7974"
    )) +
    scale_colour_manual(values = c(
      "DMR CpGs"            = "#20808D",
      "non-DMR tested CpGs" = "#B2182B",
      "all tested CpGs"     = "#7A7974"
    )) +
    labs(
      x = expression(log[10]~"(standard deviation of methylation beta)"),
      y = "density",
      title = paste0("per-CpG methylation SD: ", contrast),
      fill = NULL, colour = NULL
    ) +
    theme_classic(base_size=11) +
    theme(
      legend.position = "top",
      plot.title = element_text(size=11, face="plain"),
      axis.title = element_text(size=10),
      axis.text  = element_text(size=9)
    )

  out_pdf <- file.path(OUT_DIR, paste0("grant_fig1a_sd_density_", contrast, ".pdf"))
  out_png <- file.path(OUT_DIR, paste0("grant_fig1a_sd_density_", contrast, ".png"))
  ggsave(out_pdf, p, width=6, height=4)
  ggsave(out_png, p, width=6, height=4, dpi=300)
  message("wrote: ", out_pdf)
}

message("done.")
