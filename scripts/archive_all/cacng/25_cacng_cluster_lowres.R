.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(DMRcaller)
  library(ggplot2)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")
OUT_DIR <- "results/lowres_profiles"

COND_COLOURS <- c(ASO_CTRL="#1F3A5F", Scramble_CTRL="#6B7280",
                  ASO_VPA="#C0392B", Scramble_VPA="#D4A017")
CONDITIONS <- c("ASO_CTRL","Scramble_CTRL","ASO_VPA","Scramble_VPA")
BY_CHR <- "results/alignments/bs/by_chr"

# CACNG cluster: chr17:66.2-67.1 Mb (covers PRKCA, CACNG5, CACNG4, CACNG1)
REGION_START <- 66800000
REGION_END   <- 67100000
WIN_SIZE     <- 300

# Gene annotation for plot
GENES <- data.frame(
  name  = c("CACNG5",  "CACNG4",  "CACNG1"),
  start = c(66835117,  66964707,  67044554),
  end   = c(66894751,  67033398,  67056797)
)

message("Loading chr17 CpG data...")
pooled <- lapply(CONDITIONS, function(cond) {
  files <- file.path(BY_CHR,
    sprintf("%s_%d_chr17.CpG_report.txt.gz", cond, 1:3))
  files <- files[file.exists(files)]
  if (length(files)==0) return(NULL)
  grs <- lapply(files, function(f) {
    d <- read.table(gzfile(f), header=FALSE, sep="\t",
      col.names=c("chr","pos","strand","M","U","ctx","tri"))
    d <- d[d$ctx=="CG",]
    GRanges(d$chr, IRanges(d$pos,d$pos), strand=d$strand,
            readsM=d$M, readsN=d$M+d$U, context=d$ctx,
            trinucleotide_context=d$tri)
  })
  poolMethylationDatasets(GRangesList(grs))
})
names(pooled) <- CONDITIONS
pooled <- Filter(Negate(is.null), pooled)

region <- GRanges("chr17", IRanges(REGION_START, REGION_END))

df_all <- do.call(rbind, lapply(names(pooled), function(cond) {
  prof <- computeMethylationProfile(pooled[[cond]], region, WIN_SIZE, "CG")
  df <- as.data.frame(prof)
  df$meth <- df$sumReadsM / df$sumReadsN
  df$pos  <- (df$start + df$end) / 2
  df$cond <- cond
  df[!is.na(df$meth) & df$sumReadsN >= 3, ]
}))
df_all$cond <- factor(df_all$cond, levels=CONDITIONS)

# Gene track panel (top)
p_genes <- ggplot(GENES) +
  geom_rect(aes(xmin=start/1e6, xmax=end/1e6, ymin=0.2, ymax=0.8),
    fill="#1F3A5F", colour="grey30", linewidth=0.3) +
  geom_text(aes(x=(start+end)/2/1e6, y=0.5, label=name),
    size=3.5, colour="white", fontface="bold") +
  scale_x_continuous(limits=c(REGION_START/1e6, REGION_END/1e6),
    labels=function(x) paste0(round(x,2)," Mb")) +
  scale_y_continuous(limits=c(0,1)) +
  theme_void(base_size=11) +
  theme(axis.text.x=element_text(size=8, colour="grey40"),
        plot.margin=margin(2,5,0,5))

# Methylation panel (bottom)
p <- ggplot(df_all, aes(x=pos/1e6, y=meth, colour=cond, linetype=cond)) +
  geom_line(linewidth=0.4, alpha=0.2) +
  geom_smooth(aes(group=cond), method="loess", span=0.15,
              se=FALSE, linewidth=1.2) +
  scale_colour_manual(values=COND_COLOURS, name=NULL) +
  scale_linetype_manual(
    values=c(ASO_CTRL="solid", Scramble_CTRL="dashed",
             ASO_VPA="solid", Scramble_VPA="dashed"), name=NULL) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2)) +
  scale_x_continuous(labels=function(x) paste0(round(x,2)," Mb")) +
  theme_classic(base_size=13) +
  theme(legend.position="right",
        plot.title=element_text(face="bold", size=12),
        plot.subtitle=element_text(colour="grey30", size=9),
        panel.grid.major.y=element_line(colour="grey92")) +
  labs(x="chr17 position (Mb)", y="CpG methylation (proportion)")

library(patchwork)
final <- p_genes / p +
  plot_layout(heights=c(0.15, 1)) +
  plot_annotation(
    title="chr17 CACNG cluster: combined therapy hotspot",
    subtitle=paste0("CACNG5, CACNG4, CACNG1 — AMPA receptor auxiliary subunits | ",
                    "chr17:66.8-67.1 Mb | 300bp bins + loess smoothing\n",
                    "Solid=ASO conditions, dashed=Scramble controls | ",
                    "ASO+VPA (red) shows coordinated hypomethylation across cluster"),
    theme=theme(plot.title=element_text(face="bold", size=13),
                plot.subtitle=element_text(colour="grey30", size=9)))
ggsave(file.path(OUT_DIR, "lowres_chr17_CACNG_cluster_300bp.pdf"),
       final, width=12, height=6, device=cairo_pdf)
message("saved: lowres_chr17_CACNG_cluster_300bp.pdf")
message("Done.")
