#!/usr/bin/env Rscript
.libPaths("~/R/library")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(karyoploteR)
  library(circlize)
})
setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")

OUT <- "results/figures/genomic_distribution"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

HYPO  <- "#1F3A5F"  # navy
HYPER <- "#C0392B"  # red
CHR_LEVELS <- paste0("chr", c(1:22, "X", "Y"))

CONTRASTS <- list(
  list(name="ASO_CTRL_vs_Scramble_CTRL",    title="ASO alone"),
  list(name="Scramble_VPA_vs_Scramble_CTRL", title="VPA alone"),
  list(name="ASO_VPA_vs_Scramble_VPA",       title="ASO in VPA context"),
  list(name="ASO_VPA_vs_ASO_CTRL",           title="VPA in ASO context")
)


for (ct in CONTRASTS) {
  rds <- paste0("results/dmr/dmr_", ct$name, ".rds")
  if (!file.exists(rds)) next
  dmrs <- readRDS(rds)
  dmrs <- dmrs[seqnames(dmrs) %in% CHR_LEVELS]

  hypo  <- dmrs[dmrs$regionType == "loss"]
  hyper <- dmrs[dmrs$regionType == "gain"]

  out_pdf <- file.path(OUT, paste0("ideogram_linear_", ct$name, ".pdf"))
  pdf(out_pdf, width=14, height=8, bg="white")
  kp <- plotKaryotype(genome="hg38", plot.type=1,
                      main=paste("DMR distribution:", ct$title),
                      cex=0.7)
  kpPlotDensity(kp, data=hypo,  window.size=5e6,
                col=HYPO,  border=NA, r0=0.55, r1=1.0)
  kpPlotDensity(kp, data=hyper, window.size=5e6,
                col=HYPER, border=NA, r0=0.0,  r1=0.45)
  kpAddLabels(kp, "Hypo",  r0=0.55, r1=1.0,  cex=0.7, col=HYPO)
  kpAddLabels(kp, "Hyper", r0=0.0,  r1=0.45, cex=0.7, col=HYPER)
  dev.off()
  message("Saved linear ideogram: ", basename(out_pdf))
}


for (ct in CONTRASTS) {
  rds <- paste0("results/dmr/dmr_", ct$name, ".rds")
  if (!file.exists(rds)) next
  dmrs <- as.data.frame(readRDS(rds))
  colnames(dmrs)[1] <- "chr"
  dmrs <- dmrs[dmrs$chr %in% CHR_LEVELS, ]

  hypo  <- dmrs[dmrs$regionType == "loss",  ]
  hyper <- dmrs[dmrs$regionType == "gain", ]

  out_pdf <- file.path(OUT, paste0("circos_", ct$name, ".pdf"))
  pdf(out_pdf, width=8, height=8, bg="white")

  circos.clear()
  circos.par(
    "start.degree" = 90,
    "gap.degree"   = 1.5,
    "track.height" = 0.1,
    cell.padding   = c(0, 0, 0, 0)
  )

  # initialise with hg38 chromosome sizes
  circos.initializeWithIdeogram(
    species        = "hg38",
    chromosome.index = CHR_LEVELS,
    plotType       = c("ideogram", "labels"),
    ideogram.height = 0.03
  )

  # Hypo density track
  if (nrow(hypo) > 0) {
    circos.genomicDensity(
      hypo[, c("chr","start","end")],
      col          = adjustcolor(HYPO, 0.8),
      track.height = 0.12,
      baseline     = 0,
      bg.border    = NA
    )
  }

  # Hyper density track
  if (nrow(hyper) > 0) {
    circos.genomicDensity(
      hyper[, c("chr","start","end")],
      col          = adjustcolor(HYPER, 0.8),
      track.height = 0.12,
      baseline     = 0,
      bg.border    = NA
    )
  }

  # Title and legend
  title(main=paste("DMR distribution:", ct$title), cex.main=1.1, font.main=2)
  legend("bottomright",
         legend = c("Hypomethylated", "Hypermethylated"),
         fill   = c(HYPO, HYPER),
         bty    = "n", cex = 0.9)

  circos.clear()
  dev.off()
  message("Saved circos: ", basename(out_pdf))
}


out_combined <- file.path(OUT, "circos_4contrasts_combined.pdf")
pdf(out_combined, width=16, height=16, bg="white")
par(mfrow=c(2,2), mar=c(1,1,2,1))

for (ct in CONTRASTS) {
  rds <- paste0("results/dmr/dmr_", ct$name, ".rds")
  if (!file.exists(rds)) next
  dmrs <- as.data.frame(readRDS(rds))
  colnames(dmrs)[1] <- "chr"
  dmrs <- dmrs[dmrs$chr %in% CHR_LEVELS, ]

  hypo  <- dmrs[dmrs$regionType == "loss", ]
  hyper <- dmrs[dmrs$regionType == "gain", ]

  circos.clear()
  circos.par("start.degree"=90, "gap.degree"=1.5,
             "track.height"=0.1, cell.padding=c(0,0,0,0))
  circos.initializeWithIdeogram(
    species          = "hg38",
    chromosome.index = CHR_LEVELS,
    plotType         = c("ideogram", "labels"),
    ideogram.height  = 0.03
  )
  if (nrow(hypo) > 0)
    circos.genomicDensity(hypo[,c("chr","start","end")],
                          col=adjustcolor(HYPO,0.8),
                          track.height=0.12, baseline=0, bg.border=NA)
  if (nrow(hyper) > 0)
    circos.genomicDensity(hyper[,c("chr","start","end")],
                          col=adjustcolor(HYPER,0.8),
                          track.height=0.12, baseline=0, bg.border=NA)
  title(ct$title, cex.main=1.2, font.main=2)
  circos.clear()
}
dev.off()
message("Saved combined circos: ", basename(out_combined))


all_bar_plots <- lapply(CONTRASTS, function(ct) {
  rds <- paste0("results/dmr/dmr_", ct$name, ".rds")
  if (!file.exists(rds)) return(NULL)
  dmrs <- as.data.frame(readRDS(rds))
  colnames(dmrs)[1] <- "chr"

  counts <- dmrs %>%
    filter(chr %in% CHR_LEVELS, chr != "chrY") %>%
    mutate(direction = ifelse(regionType == "gain", "Hypo", "Hyper")) %>%
    count(chr, direction) %>%
    mutate(chr      = factor(chr, levels=CHR_LEVELS),
           n_signed = ifelse(direction == "Hypo", -n, n),
           chr = factor(chr, levels=CHR_LEVELS[CHR_LEVELS != "chrY"]))

  ggplot(counts, aes(x=chr, y=n_signed, fill=direction)) +
    geom_bar(stat="identity") +
    geom_hline(yintercept=0, linewidth=0.4, colour="grey30") +
    scale_fill_manual(values=c(Hypo=HYPO, Hyper=HYPER)) +
    scale_y_continuous(labels=function(x) ifelse(x<0, paste0("-",format(abs(x),big.mark=",")), format(x,big.mark=","))) +
    labs(x=NULL, y="DMR count", fill=NULL,
         title=sprintf("%s vs %s", ct$cond_b, ct$cond_a)) +
    theme_classic(base_size=12) +
    theme(axis.text.x     = element_text(angle=45, hjust=1, size=7),
          plot.title      = element_text(face="bold", size=10),
          legend.position = "top")
})

bar_plots_filtered <- Filter(Negate(is.null), all_bar_plots)
for (i in seq_along(bar_plots_filtered)) {
  bar_plots_filtered[[i]] <- bar_plots_filtered[[i]] +
    labs(tag=c("A","B","C","D")[i]) +
    theme(plot.tag=element_text(face="bold", size=14))
}
combined_bars <- wrap_plots(bar_plots_filtered, ncol=2) +
  plot_annotation(theme=theme(plot.title=element_text(face="bold", size=13)))

ggsave(file.path(OUT, "chr_dmr_diverging_4contrasts.pdf"),
       combined_bars, width=16, height=10, device=cairo_pdf)
ggsave(file.path(OUT, "chr_dmr_diverging_4contrasts.png"),
       combined_bars, width=16, height=10, dpi=150)
message("Saved diverging bars")
message("All done.")
