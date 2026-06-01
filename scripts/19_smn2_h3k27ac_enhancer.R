#!/usr/bin/env Rscript
# 19_smn2_h3k27ac_enhancer.R
# H3K27ac peak analysis at SMN2 locus — addresses Alberto Kornblihtt's
# question about transcriptional enhancers at the 3' end of SMN2.
# Data: GSE246399 (HEK293T, CTRL and VPA-treated, Calandrelli et al.)
# Cross-referenced with our WGBS DMR results at the SMN2 locus.

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(patchwork)
  library(rtracklayer)
})
.libPaths(c("~/R/library", .libPaths()))

setwd("/data/home/bt25018/sma_epigenomics_pipeline")
OUT <- "results/smn2_enhancer"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

# SMN2 gene body and extended window
SMN2_START  <- 70049638
SMN2_END    <- 70078522
WINDOW_START <- 69950000
WINDOW_END   <- 70150000
CHR <- "chr5"

# Exon boundaries (Alberto's convention)
EXONS <- data.frame(
  exon  = c("E1","E2a","E2b","E3","E4","E5","E6","E7","E8"),
  start = c(70049638,70053107,70056229,70063044,70069090,
            70069235,70070641,70076521,70077019),
  end   = c(70050437,70053264,70056357,70063153,70069186,
            70069330,70070751,70076574,70077592)
)

# Load narrowPeak files
load_peaks <- function(file, label) {
  df <- read.table(gzfile(file), header=FALSE, sep="\t",
    col.names=c("chr","start","end","name","score",
                "strand","fc","pval","qval","summit"))
  df <- df[df$chr==CHR & df$start>=WINDOW_START & df$end<=WINDOW_END,]
  df$label <- label
  df$mid <- (df$start + df$end) / 2
  df
}

peaks_ctrl1 <- load_peaks(
  "data/external/h3k27ac_gse246399/H3K27ac_CTRL_Rep1.narrowPeak.gz",
  "CTRL Rep1")
peaks_ctrl2 <- load_peaks(
  "data/external/h3k27ac_gse246399/H3K27ac_CTRL_Rep2.narrowPeak.gz",
  "CTRL Rep2")
peaks_vpa1  <- load_peaks(
  "data/external/h3k27ac_gse246399/H3K27ac_VPA_Rep1.narrowPeak.gz",
  "VPA Rep1")

all_peaks <- rbind(peaks_ctrl1, peaks_ctrl2, peaks_vpa1)
all_peaks$label <- factor(all_peaks$label,
  levels=c("CTRL Rep1","CTRL Rep2","VPA Rep1"))

message("Peaks found:")
print(table(all_peaks$label))

# Save peak table
write.csv(all_peaks[,c("chr","start","end","fc","qval","label")],
  file.path(OUT, "smn2_h3k27ac_peaks.csv"), row.names=FALSE)
message("Saved: smn2_h3k27ac_peaks.csv")

# ── FIGURE: H3K27ac peaks at SMN2 locus ──────────────────────
COLS <- c("CTRL Rep1"="#1D6FA4", "CTRL Rep2"="#2E9B6F", "VPA Rep1"="#F0A500")

# Panel A: peak map showing position and strength
p_peaks <- ggplot(all_peaks,
    aes(xmin=start/1e6, xmax=end/1e6,
        ymin=as.numeric(label)-0.35,
        ymax=as.numeric(label)+0.35,
        fill=label)) +
  geom_rect(alpha=0.8) +
  # Gene body
  annotate("rect", xmin=SMN2_START/1e6, xmax=SMN2_END/1e6,
           ymin=0.3, ymax=0.7, fill="grey80", colour="grey50", linewidth=0.3) +
  annotate("text", x=(SMN2_START+SMN2_END)/2/1e6, y=0.5,
           label="SMN2", size=3, fontface="bold") +
  # Exons
  geom_rect(data=EXONS,
    aes(xmin=start/1e6, xmax=end/1e6, ymin=0.2, ymax=0.8),
    inherit.aes=FALSE, fill="#1B4F8A", colour="white", linewidth=0.2) +
  # E7 in red
  annotate("rect",
    xmin=70076521/1e6, xmax=70076574/1e6,
    ymin=0.15, ymax=0.85, fill="#D94F3D", colour="white", linewidth=0.2) +
  annotate("text", x=70076547/1e6, y=0.05,
           label="E7\n(ASO target)", size=2, colour="#D94F3D") +
  scale_fill_manual(values=COLS, name=NULL) +
  scale_y_continuous(breaks=1:3,
    labels=c("CTRL Rep1","CTRL Rep2","VPA Rep1")) +
  scale_x_continuous(labels=function(x) paste0(x, " Mb")) +
  coord_cartesian(xlim=c(WINDOW_START/1e6, WINDOW_END/1e6)) +
  labs(title="(A) H3K27ac ChIP-seq peaks at SMN2 locus",
       subtitle=paste0("GSE246399 (HEK293T, Calandrelli et al.). ",
                       "No peaks at SMN2 3' end in untreated cells."),
       x="chr5 position", y=NULL) +
  theme_bw(base_size=10) +
  theme(legend.position="none",
        plot.title=element_text(face="bold", size=9),
        plot.subtitle=element_text(size=8),
        panel.grid.minor=element_blank())

# Panel B: fold enrichment at each peak coloured by condition
p_fc <- ggplot(all_peaks,
    aes(x=mid/1e6, y=fc, colour=label, shape=label)) +
  geom_point(size=3, alpha=0.9) +
  geom_vline(xintercept=c(SMN2_START/1e6, SMN2_END/1e6),
             linetype="dashed", colour="grey60", linewidth=0.4) +
  annotate("text", x=SMN2_START/1e6, y=max(all_peaks$fc)*0.95,
           label="SMN2\n5'", size=2.5, hjust=1.1, colour="grey40") +
  annotate("text", x=SMN2_END/1e6, y=max(all_peaks$fc)*0.95,
           label="SMN2\n3'", size=2.5, hjust=-0.1, colour="grey40") +
  annotate("rect", xmin=70076521/1e6, xmax=70076574/1e6,
           ymin=-Inf, ymax=Inf, fill="#D94F3D", alpha=0.1) +
  scale_colour_manual(values=COLS, name=NULL) +
  scale_shape_manual(values=c(16,17,15), name=NULL) +
  scale_x_continuous(labels=function(x) paste0(x, " Mb")) +
  coord_cartesian(xlim=c(WINDOW_START/1e6, WINDOW_END/1e6)) +
  labs(title="(B) H3K27ac fold enrichment by peak position",
       subtitle="Red shading = E7 (ASO target). Dashed lines = SMN2 gene boundaries.",
       x="chr5 position", y="Fold enrichment over input") +
  theme_bw(base_size=10) +
  theme(legend.position="top",
        plot.title=element_text(face="bold", size=9),
        plot.subtitle=element_text(size=8),
        panel.grid.minor=element_blank())

fig <- (p_peaks / p_fc) +
  plot_annotation(
    title="SMN2 locus H3K27ac enhancer analysis",
    subtitle=paste0(
      "H3K27ac peaks mark the SMN2 promoter region and an upstream element (~70.025 Mb) ",
      "in untreated HEK293T cells.\n",
      "No H3K27ac peaks detected at the SMN2 3' end in baseline conditions. ",
      "VPA treatment gains additional gene body peaks\n",
      "consistent with HDAC inhibitor-driven chromatin opening, ",
      "but not at the 3' end specifically."),
    theme=theme(plot.title=element_text(face="bold", size=11),
                plot.subtitle=element_text(size=8, colour="grey30"))
  )

ggsave(file.path(OUT, "SMN2_H3K27ac_enhancer_analysis.pdf"),
       fig, width=14, height=10)
ggsave(file.path(OUT, "SMN2_H3K27ac_enhancer_analysis.png"),
       fig, width=14, height=10, dpi=300)
message("Saved: SMN2_H3K27ac_enhancer_analysis")

# Summary table for thesis
summary_df <- all_peaks[,c("chr","start","end","fc","qval","label")]
summary_df$position <- ifelse(summary_df$end < SMN2_START, "upstream",
                       ifelse(summary_df$start > SMN2_END, "downstream",
                       ifelse(summary_df$start < SMN2_START+5000, "promoter/5-end",
                       "gene body")))
write.csv(summary_df, file.path(OUT, "smn2_h3k27ac_peak_summary.csv"),
  row.names=FALSE)
message("\nPeak summary:")
print(table(summary_df$label, summary_df$position))
message("\nDone. Results in: ", OUT)
