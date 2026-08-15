#!/usr/bin/env Rscript
# fig_lowres_chr1.R
# Low-resolution methylation profile across chr1 at 1Mb windows
# Four conditions, canonical colours, solid/dashed linetype distinction

.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
  library(ggplot2)
  library(dplyr)
})

setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")
OUT <- "results/thesis_figures"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

COND_COLS <- c(
  Scramble_CTRL = "#6B7280",
  ASO_CTRL      = "#1F3A5F",
  Scramble_VPA  = "#F0A500",
  ASO_VPA       = "#C0392B"
)

# solid = Scramble (reference), dashed = ASO (treatment)
COND_LTY <- c(
  Scramble_CTRL = "solid",
  ASO_CTRL      = "dashed",
  Scramble_VPA  = "solid",
  ASO_VPA       = "dashed"
)

WIN_SIZE  <- 1000000  # 1 Mb bins
CHR       <- "chr1"
MIN_READS <- 3

message("Loading methylation cache...")
meth <- readRDS("results/dmr/meth_pooled_cache.rds")

message("Computing profiles...")
profiles <- lapply(names(COND_COLS), function(cond) {
  m <- meth[[cond]]
  m <- m[as.character(seqnames(m)) == CHR & m$context == "CG" & m$readsN >= MIN_READS]
  
  # bin by 1Mb windows
  pos  <- start(m)
  bin  <- floor(pos / WIN_SIZE) * WIN_SIZE + WIN_SIZE/2  # bin midpoint
  
  df <- data.frame(
    pos   = pos,
    bin   = bin,
    readsM = m$readsM,
    readsN = m$readsN
  )
  
  # coverage-weighted mean per bin
  agg <- df %>%
    group_by(bin) %>%
    summarise(
      sumM  = sum(readsM),
      sumN  = sum(readsN),
      nCpGs = n(),
      .groups = "drop"
    ) %>%
    filter(sumN >= 10) %>%  # minimum reads per bin
    mutate(
      meth      = 100 * sumM / sumN,
      condition = cond
    )
  agg
})

prof_df <- do.call(rbind, profiles)
prof_df$condition <- factor(prof_df$condition, levels=names(COND_COLS))

p <- ggplot(prof_df, aes(x=bin/1e6, y=meth,
                          colour=condition, linetype=condition)) +
  geom_line(linewidth=0.7, alpha=0.9) +
  scale_colour_manual(values=COND_COLS, name=NULL) +
  scale_linetype_manual(values=COND_LTY, name=NULL) +
  scale_x_continuous(
    breaks=seq(0, 250, by=50),
    labels=function(x) paste0(x, " Mb"),
    expand=c(0.01, 0)
  ) +
  scale_y_continuous(
    limits=c(0, 100),
    breaks=seq(0, 100, by=20),
    labels=function(x) paste0(x, "%")
  ) +
  labs(
    x = "chr1 position",
    y = "Mean CpG methylation (%)"
  ) +
  theme_classic(base_size=12) +
  theme(
    legend.position    = "right",
    legend.text        = element_text(size=10),
    panel.grid.major.y = element_line(colour="grey92", linewidth=0.3),
    axis.text          = element_text(colour="#1A1A1A"),
    axis.title         = element_text(colour="#1A1A1A")
  )

out <- file.path(OUT, "Fig_lowres_chr1_1Mb.pdf")
ggsave(out, p, width=12, height=4, device=cairo_pdf)
message("Saved: ", out)
