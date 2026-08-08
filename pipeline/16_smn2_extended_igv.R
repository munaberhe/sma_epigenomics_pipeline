#!/usr/bin/env Rscript
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(DMRcaller)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")

OUT <- "results/thesis_figures/smn2_extended_igv"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)


SMN2_START <- 70049638
SMN2_END   <- 70078522
FLANK      <- 50000
PLOT_START <- SMN2_START - FLANK
PLOT_END   <- SMN2_END   + FLANK
CHR        <- "chr5"

REGION <- GRanges(CHR, IRanges(PLOT_START, PLOT_END))


COND_COLS <- c(
  ASO_CTRL      = "#1F3A5F",
  Scramble_CTRL = "#6B7280",
  ASO_VPA       = "#C0392B",
  Scramble_VPA  = "#F0A500"
)
CONTRAST_COLS <- c(
  "ASO alone"  = "#1F3A5F",
  "VPA alone"  = "#F0A500",
  "ASO in VPA" = "#C0392B",
  "VPA in ASO" = "#8E44AD"
)
ENH_H9_COL   <- "#D4820A"
ENH_CCRE_COL <- "#0072B2"
CGI_COL      <- "#2196A6"
H3K27_CTRL   <- "#1F3A5F"
H3K27_VPA    <- "#6B2D8B"

message("Loading annotation tracks...")


cgi_raw <- read.table("data/reference/cpg_islands_hg38.bed",
                      header=FALSE, sep="\t", stringsAsFactors=FALSE)
cgi <- GRanges(cgi_raw[,1], IRanges(cgi_raw[,2], cgi_raw[,3]))
cgi <- as.data.frame(subsetByOverlaps(cgi, REGION))
message("  CpG islands: ", nrow(cgi))


h9_raw <- read.table(gzfile("data/reference/H9_predicted_non_promoter_non_fragments.bed.gz"),
                     header=TRUE, sep="\t", stringsAsFactors=FALSE)
h9_enh <- GRanges(h9_raw$seqnames, IRanges(h9_raw$start, h9_raw$end))
h9_enh <- as.data.frame(subsetByOverlaps(h9_enh, REGION))
message("  H9 enhancers: ", nrow(h9_enh))


ccre_raw <- read.table("data/reference/encode_cCREs_hg38.bed",
                       header=FALSE, sep="\t", stringsAsFactors=FALSE)
ccre <- GRanges(ccre_raw[,1], IRanges(ccre_raw[,2], ccre_raw[,3]))
ccre <- as.data.frame(subsetByOverlaps(ccre, REGION))
message("  ENCODE cCREs: ", nrow(ccre))


h3k27_ctrl <- tryCatch({
  r1 <- import("data/external/h3k27ac_gse246399/H3K27ac_CTRL_Rep1.narrowPeak.gz", format="narrowPeak")
  r2 <- import("data/external/h3k27ac_gse246399/H3K27ac_CTRL_Rep2.narrowPeak.gz", format="narrowPeak")
  as.data.frame(subsetByOverlaps(c(r1,r2), REGION))
}, error=function(e) { message("H3K27ac CTRL not loaded: ", e$message); data.frame() })

h3k27_vpa <- tryCatch({
  r1 <- import("data/external/h3k27ac_gse246399/H3K27ac_VPA_Rep1.narrowPeak.gz", format="narrowPeak")
  as.data.frame(subsetByOverlaps(r1, REGION))
}, error=function(e) { message("H3K27ac VPA not loaded: ", e$message); data.frame() })

message("  H3K27ac CTRL peaks: ", nrow(h3k27_ctrl))
message("  H3K27ac VPA peaks:  ", nrow(h3k27_vpa))


CONTRASTS <- list(
  list(name="ASO_CTRL_vs_Scramble_CTRL",   label="ASO alone"),
  list(name="Scramble_VPA_vs_Scramble_CTRL", label="VPA alone"),
  list(name="ASO_VPA_vs_Scramble_VPA",     label="ASO in VPA"),
  list(name="ASO_VPA_vs_ASO_CTRL",         label="VPA in ASO")
)

dmr_df <- do.call(rbind, lapply(CONTRASTS, function(ct) {
  rds <- paste0("results/dmr/dmr_", ct$name, ".rds")
  if (!file.exists(rds)) return(NULL)
  d <- as.data.frame(subsetByOverlaps(readRDS(rds), REGION))
  if (nrow(d) == 0) return(NULL)
  d$contrast <- ct$label
  d
}))

if (!is.null(dmr_df) && nrow(dmr_df) > 0) {
  dmr_df$contrast <- factor(dmr_df$contrast, levels=names(CONTRAST_COLS))
  message("  Total DMRs in region: ", nrow(dmr_df))
} else {
  message("  No DMRs in region")
  dmr_df <- data.frame(seqnames=character(), start=numeric(), end=numeric(), contrast=character())
}


# Load sensitive DMRs (unmasked canonical data, all contrasts)
sens_dmr <- tryCatch({
  csv <- "results/smn2_local_dmr/SMN2_sensitive_DMRs_all_contrasts.csv"
  if (file.exists(csv)) {
    d <- read.csv(csv, stringsAsFactors=FALSE)
    d <- d[d$seqnames == "chr5" &
           d$start >= PLOT_START &
           d$end   <= PLOT_END, ]
    if (nrow(d) > 0) d else NULL
  } else NULL
}, error=function(e) { message("Sensitive DMRs not loaded: ", e$message); NULL })
message("  Sensitive DMRs in region: ", if (!is.null(sens_dmr)) nrow(sens_dmr) else 0)

message("Loading methylation cache...")
meth_cache <- readRDS("results/dmr/meth_pooled_cache.rds")

WIN <- 1000
meth_df <- do.call(rbind, lapply(names(COND_COLS), function(cond) {
  tryCatch({
    prof <- computeMethylationProfile(meth_cache[[cond]], REGION,
                                      windowSize=WIN, context="CG")
    df <- as.data.frame(prof)
    df <- df[!is.na(df$sumReadsM) & df$sumReadsN >= 3, ]
    if (nrow(df) == 0) return(NULL)
    df$meth <- df$sumReadsM / df$sumReadsN
    df$pos  <- (df$start + df$end) / 2
    df$condition <- cond
    df[, c("pos","meth","condition")]
  }, error=function(e) {
    message("  meth profile failed for ", cond, ": ", e$message)
    NULL
  })
}))
meth_df$condition <- factor(meth_df$condition, levels=names(COND_COLS))


message("Parsing exon structure...")
gtf_cmd <- paste0("zcat data/reference/Homo_sapiens.GRCh38.109.chr.gtf.gz | ",
                  "grep -w SMN2 | awk '$3==\"exon\"'")
gtf_raw <- read.table(pipe(gtf_cmd), sep="\t", stringsAsFactors=FALSE)
exons <- data.frame(
  start = gtf_raw[,4],
  end   = gtf_raw[,5]
)
# deduplicate exons
exons <- unique(exons)
exons <- exons[exons$start >= PLOT_START & exons$end <= PLOT_END, ]

# intron backbone
intron_df <- data.frame(
  x    = PLOT_START,
  xend = PLOT_END,
  y    = 0.5,
  yend = 0.5
)

message("  Exons: ", nrow(exons))


theme_track <- function() {
  theme_classic(base_size=14) +
  theme(
    axis.title.x    = element_blank(),
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank(),
    axis.line.x     = element_blank(),
    plot.margin     = margin(2,10,2,10),
    legend.position = "none"
  )
}

mb <- function(x) paste0(round(x/1e6, 2), " Mb")

# Track 1: Gene model
p_gene <- ggplot() +
  geom_segment(data=intron_df, aes(x=x, xend=xend, y=y, yend=y),
               colour="grey40", linewidth=0.8) +
  annotate("segment", x=SMN2_END-2000, xend=SMN2_END+1000,
           y=0.5, yend=0.5, colour="#2C3E50", linewidth=1.2,
           arrow=arrow(length=unit(0.2,"cm"), type="closed")) +
  annotate("text", x=(SMN2_START+SMN2_END)/2, y=0.75,
           label="SMN2", fontface="bold", size=5, colour="#2C3E50") +
  geom_vline(xintercept=c(SMN2_START, SMN2_END),
             linetype="dashed", colour="grey50", linewidth=0.4) +
  scale_x_continuous(limits=c(PLOT_START, PLOT_END)) +
  scale_y_continuous(limits=c(0,1.1)) +
  labs(y="Gene") +
  theme_track() +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.title.y=element_text(size=11, colour="grey30"))

# Track 2: Regulatory annotations
annot_df <- rbind(
  if (nrow(cgi)>0)    data.frame(start=cgi$start,    end=cgi$end,    type="CpG island", y=3) else NULL,
  if (nrow(h9_enh)>0) data.frame(start=h9_enh$start, end=h9_enh$end, type="H9 enhancer", y=2) else NULL,
  if (nrow(ccre)>0)   data.frame(start=ccre$start,   end=ccre$end,   type="ENCODE cCRE", y=1) else NULL
)
annot_cols <- c("CpG island"=CGI_COL, "H9 enhancer"=ENH_H9_COL, "ENCODE cCRE"=ENH_CCRE_COL)

p_annot <- ggplot() +
  { if (!is.null(annot_df) && nrow(annot_df)>0)
      geom_rect(data=annot_df,
                aes(xmin=start, xmax=end, ymin=y-0.35, ymax=y+0.35, fill=type),
                colour=NA)
    else geom_blank() } +
  geom_vline(xintercept=c(SMN2_START, SMN2_END),
             linetype="dashed", colour="grey50", linewidth=0.4) +
  scale_fill_manual(values=annot_cols) +
  scale_x_continuous(limits=c(PLOT_START, PLOT_END)) +
  scale_y_continuous(limits=c(0.5, 3.5),
                     breaks=c(1,2,3),
                     labels=c("ENCODE cCRE","H9 enhancer","CpG island")) +
  labs(y=NULL) +
  theme_track() +
  theme(axis.text.y=element_text(size=9, colour="grey30"),
        axis.ticks.y=element_blank())

# Track 3: H3K27ac CTRL
p_h3k_ctrl <- ggplot() +
  { if (nrow(h3k27_ctrl)>0)
      geom_rect(data=h3k27_ctrl,
                aes(xmin=start, xmax=end, ymin=0, ymax=1),
                fill=H3K27_CTRL, colour=NA, alpha=0.8)
    else geom_blank() } +
  geom_vline(xintercept=c(SMN2_START, SMN2_END),
             linetype="dashed", colour="grey50", linewidth=0.4) +
  scale_x_continuous(limits=c(PLOT_START, PLOT_END)) +
  scale_y_continuous(limits=c(0,1)) +
  labs(y="H3K27ac\nCTRL") +
  theme_track() +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.title.y=element_text(size=10, colour=H3K27_CTRL))

# Track 4: H3K27ac VPA
p_h3k_vpa <- ggplot() +
  { if (nrow(h3k27_vpa)>0)
      geom_rect(data=h3k27_vpa,
                aes(xmin=start, xmax=end, ymin=0, ymax=1),
                fill=H3K27_VPA, colour=NA, alpha=0.8)
    else geom_blank() } +
  geom_vline(xintercept=c(SMN2_START, SMN2_END),
             linetype="dashed", colour="grey50", linewidth=0.4) +
  scale_x_continuous(limits=c(PLOT_START, PLOT_END)) +
  scale_y_continuous(limits=c(0,1)) +
  labs(y="H3K27ac\nVPA") +
  theme_track() +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.title.y=element_text(size=10, colour=H3K27_VPA))

# Track 5: DMR tracks
p_dmr <- ggplot() +
  { if (nrow(dmr_df)>0)
      geom_rect(data=dmr_df,
                aes(xmin=start, xmax=end,
                    ymin=as.numeric(contrast)-0.35,
                    ymax=as.numeric(contrast)+0.35,
                    fill=contrast),
                colour=NA, alpha=0.85)
    else geom_blank() } +
  geom_vline(xintercept=c(SMN2_START, SMN2_END),
             linetype="dashed", colour="grey50", linewidth=0.4) +
  scale_fill_manual(values=CONTRAST_COLS) +
  scale_x_continuous(limits=c(PLOT_START, PLOT_END)) +
  scale_y_continuous(limits=c(0.5, 4.5),
                     breaks=1:4,
                     labels=names(CONTRAST_COLS)) +
  labs(y=NULL) +
  theme_track() +
  theme(axis.text.y=element_text(size=9, colour="grey30"),
        axis.ticks.y=element_blank())

# Track 6: Sensitive DMRs (2% threshold)
p_sens <- ggplot() +
  { if (!is.null(sens_dmr) && nrow(sens_dmr) > 0)
      geom_rect(data=sens_dmr,
                aes(xmin=start, xmax=end, ymin=0.1, ymax=0.9),
                fill="#7B2D8B", colour=NA, alpha=0.85)
    else geom_blank() } +
  geom_vline(xintercept=c(SMN2_START, SMN2_END),
             linetype="dashed", colour="grey50", linewidth=0.4) +
  scale_x_continuous(limits=c(PLOT_START, PLOT_END)) +
  scale_y_continuous(limits=c(0,1)) +
  labs(y="Sensitive
DMRs (2%)") +
  theme_track() +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.title.y=element_text(size=10, colour="#7B2D8B"))

# Track 7: Methylation profiles (filled area)
p_meth <- ggplot(meth_df, aes(x=pos, y=meth, colour=condition)) +
  geom_smooth(method="loess", span=0.15, se=FALSE, linewidth=1.2) +
  geom_vline(xintercept=c(SMN2_START, SMN2_END),
             linetype="dashed", colour="grey50", linewidth=0.4) +
  scale_colour_manual(values=COND_COLS, name=NULL) +
  scale_x_continuous(limits=c(PLOT_START, PLOT_END),
                     labels=function(x) paste0(round(x/1e6,2)," Mb")) +
  scale_y_continuous(limits=c(0,1), breaks=c(0,0.5,1),
                     labels=c("0","0.5","1")) +
  labs(y="CpG\nmethylation", x=NULL) +
  theme_classic(base_size=14) +
  theme(
    plot.margin     = margin(2,10,2,10),
    legend.position = "right",
    legend.text     = element_text(size=11),
    axis.title.y    = element_text(size=11, colour="grey30"),
    axis.text.x     = element_blank(),
    axis.ticks.x    = element_blank()
  )

# Track 7: Exon structure at bottom (IGV-style)
p_exon <- ggplot() +
  geom_segment(data=intron_df, aes(x=x, xend=xend, y=y, yend=y),
               colour="grey40", linewidth=0.6) +
  geom_rect(data=exons,
            aes(xmin=start, xmax=end, ymin=0.15, ymax=0.85,
                fill=is_e7), colour=NA) +
  scale_fill_manual(values=c("FALSE"="#2C3E50", "TRUE"="#E31A1C")) +
  geom_text(data=exons,
            aes(x=(start+end)/2, y=1.05, label=exon),
            size=3, colour=ifelse(exons$is_e7, "#E31A1C", "grey30"),
            fontface=ifelse(exons$is_e7, "bold", "plain")) +
  geom_vline(xintercept=c(SMN2_START, SMN2_END),
             linetype="dashed", colour="grey50", linewidth=0.4) +
  scale_x_continuous(limits=c(PLOT_START, PLOT_END),
                     labels=function(x) paste0(round(x/1e6,2)," Mb")) +
  scale_y_continuous(limits=c(0,1.2)) +
  labs(y="Exons", x="chr5 position") +
  theme_classic(base_size=14) +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_text(size=11, colour="grey30"),
    axis.text.x  = element_text(size=10),
    legend.position = "none",
    plot.margin  = margin(2,10,5,10)
  )


layout <- "
A
B
C
D
E
F
G
H
"

p_final <- p_gene + p_annot + p_h3k_ctrl + p_h3k_vpa + p_dmr + p_sens + p_meth + p_exon +
  plot_layout(design=layout,
              heights=c(1.5, 1.5, 0.8, 0.8, 2, 0.8, 3, 1)) +
  plot_annotation(
    title    = "SMN2 locus (±50kb): CpG methylation, DMRs, H3K27ac and regulatory elements",
    subtitle = paste0("chr5:", format(PLOT_START, big.mark=","),
                      "-", format(PLOT_END, big.mark=",")),
    theme = theme(
      plot.title    = element_text(size=16, face="bold", hjust=0.5),
      plot.subtitle = element_text(size=11, hjust=0.5, colour="grey40")
    )
  )

outfile <- file.path(OUT, "SMN2_extended_IGV_50kb.pdf")
ggsave(outfile, p_final, width=16, height=14, device=cairo_pdf)
message("Saved: ", outfile)
message("All done.")
