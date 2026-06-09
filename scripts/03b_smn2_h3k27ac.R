.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(ggplot2)
  library(patchwork)
})
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
  annotate("text", x=70.025, y=3.6,
           label="SERF1B", size=2.5, colour="grey40", fontface="italic") +
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
  coord_cartesian(xlim=c(70.0, 70.085)) +
  labs(title="(A) H3K27ac ChIP-seq peaks at SMN2 locus",
       subtitle=NULL,
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
  coord_cartesian(xlim=c(70.0, 70.085)) +
  labs(title="(B) H3K27ac fold enrichment by peak position",
       subtitle=NULL,
       x="chr5 position", y="Fold enrichment over input") +
  theme_bw(base_size=10) +
  theme(legend.position="top",
        plot.title=element_text(face="bold", size=9),
        plot.subtitle=element_text(size=8),
        panel.grid.minor=element_blank())

fig <- (p_peaks / p_fc) +
  plot_annotation(
    title="SMN2 locus H3K27ac enhancer analysis",
    subtitle=NULL,
    theme=theme(plot.title=element_text(face="bold", size=11),
                plot.subtitle=element_text(size=8, colour="grey30"))
  )

ggsave(file.path(OUT, "SMN2_H3K27ac_enhancer_analysis.pdf"),
       fig, width=14, height=10)
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


# ── ITEM 9: SMN2 introns 6 and 7 H3K27ac check ───────────────────────────────
message('=== ITEM 9: SMN2 intron 6-7 H3K27ac analysis ===')

# SMN2 exon coordinates (Alberto's convention, E7 = penultimate exon)
# Intron 6 = between E6 and E7
# Intron 7 = between E7 and E8
EXONS <- data.frame(
  exon  = c("E1","E2a","E2b","E3","E4","E5","E6","E7","E8"),
  start = c(70049638,70053107,70056229,70063044,70069090,
            70069235,70070641,70076521,70077019),
  end   = c(70050437,70053264,70056357,70063153,70069186,
            70069330,70070751,70076574,70077592)
)

# Intron boundaries
intron6_start <- EXONS$end[EXONS$exon=="E6"]   # after E6
intron6_end   <- EXONS$start[EXONS$exon=="E7"] # before E7
intron7_start <- EXONS$end[EXONS$exon=="E7"]   # after E7
intron7_end   <- EXONS$start[EXONS$exon=="E8"] # before E8

message("Intron 6: chr5:", intron6_start, "-", intron6_end,
        " (", intron6_end-intron6_start, "bp)")
message("Intron 7: chr5:", intron7_start, "-", intron7_end,
        " (", intron7_end-intron7_start, "bp)")

# Load H3K27ac peaks
load_peaks <- function(file, label) {
  df <- read.table(gzfile(file), header=FALSE, sep='\t',
    col.names=c('chr','start','end','name','score',
                'strand','fc','pval','qval','summit'))
  df$label <- label
  df
}

peaks_ctrl1 <- load_peaks(
  'data/external/h3k27ac_gse246399/H3K27ac_CTRL_Rep1.narrowPeak.gz', 'CTRL_Rep1')
peaks_ctrl2 <- load_peaks(
  'data/external/h3k27ac_gse246399/H3K27ac_CTRL_Rep2.narrowPeak.gz', 'CTRL_Rep2')
peaks_vpa1  <- load_peaks(
  'data/external/h3k27ac_gse246399/H3K27ac_VPA_Rep1.narrowPeak.gz',  'VPA_Rep1')
all_peaks <- rbind(peaks_ctrl1, peaks_ctrl2, peaks_vpa1)

# Check intron 6
message('\n--- Intron 6 (E6-E7 boundary, flanks ASO target) ---')
int6 <- all_peaks[all_peaks$chr=='chr5' &
                  all_peaks$start >= intron6_start &
                  all_peaks$end   <= intron6_end, ]
if (nrow(int6) > 0) {
  message('PEAKS FOUND in intron 6:')
  print(int6[,c('chr','start','end','fc','qval','label')])
} else {
  message('No H3K27ac peaks in intron 6 in any condition')
}

# Check intron 7
message('\n--- Intron 7 (E7-E8 boundary, 3 prime end) ---')
int7 <- all_peaks[all_peaks$chr=='chr5' &
                  all_peaks$start >= intron7_start &
                  all_peaks$end   <= intron7_end, ]
if (nrow(int7) > 0) {
  message('PEAKS FOUND in intron 7:')
  print(int7[,c('chr','start','end','fc','qval','label')])
} else {
  message('No H3K27ac peaks in intron 7 in any condition')
}

# Extended check — 5kb around E7 (ASO target region)
message('\n--- 5kb window around E7 (ASO target ±5kb) ---')
e7_window <- all_peaks[all_peaks$chr=='chr5' &
                        all_peaks$start >= (EXONS$start[EXONS$exon=="E7"] - 5000) &
                        all_peaks$end   <= (EXONS$end[EXONS$exon=="E7"]   + 5000), ]
if (nrow(e7_window) > 0) {
  message('PEAKS near E7:')
  print(e7_window[,c('chr','start','end','fc','qval','label')])
} else {
  message('No H3K27ac peaks within 5kb of E7')
}

# Save intron results
intron_results <- data.frame(
  region = c('Intron 6', 'Intron 7', 'E7 ±5kb'),
  coordinates = c(
    paste0('chr5:', intron6_start, '-', intron6_end),
    paste0('chr5:', intron7_start, '-', intron7_end),
    paste0('chr5:', EXONS$start[EXONS$exon=="E7"]-5000, '-',
                    EXONS$end[EXONS$exon=="E7"]+5000)
  ),
  peaks_CTRL_Rep1 = c(
    sum(int6$label=='CTRL_Rep1'),
    sum(int7$label=='CTRL_Rep1'),
    sum(e7_window$label=='CTRL_Rep1')
  ),
  peaks_CTRL_Rep2 = c(
    sum(int6$label=='CTRL_Rep2'),
    sum(int7$label=='CTRL_Rep2'),
    sum(e7_window$label=='CTRL_Rep2')
  ),
  peaks_VPA_Rep1 = c(
    sum(int6$label=='VPA_Rep1'),
    sum(int7$label=='VPA_Rep1'),
    sum(e7_window$label=='VPA_Rep1')
  )
)
write.csv(intron_results,
  'results/smn2_enhancer/SMN2_intron67_H3K27ac_summary.csv',
  row.names=FALSE)
message('\nSaved: SMN2_intron67_H3K27ac_summary.csv')

