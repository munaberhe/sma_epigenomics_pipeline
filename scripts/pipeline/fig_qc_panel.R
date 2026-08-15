#!/usr/bin/env Rscript
# fig_qc_panel.R
# Fig 5.1 QC panel: A=coverage B=violin C=correlation D=chr1 profile

.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
  library(patchwork)
  library(dplyr)
  library(GenomicRanges)
  library(DMRcaller)
})

setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")
OUT <- "results/thesis_figures"

COND_COLS <- c(
  Scramble_CTRL="#6B7280",
  ASO_CTRL="#1F3A5F",
  Scramble_VPA="#F0A500",
  ASO_VPA="#C0392B"
)

# Panel A - coverage retention (convert PDF to PNG via magick, load as raster)
message("Panel A: coverage...")
# use existing PNG
library(png)
cov_img <- readPNG("results/qc/coverage_4lines/coverage_4lines_per_condition.png")
pA <- ggplot() +
  annotation_raster(cov_img, xmin=0, xmax=1, ymin=0, ymax=1) +
  annotate("text", x=0.02, y=0.97, label="A", fontface="bold", size=7,
           hjust=0, vjust=1) +
  coord_fixed(ratio=nrow(cov_img)/ncol(cov_img)) +
  theme_void()

# Panel B - global methylation violin
message("Panel B: violin...")
glob <- read.csv("results/global_methylation/global_methylation_per_sample.csv")
glob$condition <- factor(glob$condition, levels=names(COND_COLS))
pB <- ggplot(glob, aes(x=condition, y=mean_methylation*100, fill=condition)) +
  geom_violin(alpha=0.7) +
  geom_boxplot(width=0.15, fill="white", outlier.size=1) +
  scale_fill_manual(values=COND_COLS) +
  labs(x=NULL, y="CpG methylation (%)", tag="B") +
  theme_classic(base_size=12) +
  theme(legend.position="none",
        axis.text.x=element_text(angle=30, hjust=1),
        plot.tag=element_text(face="bold", size=14))

# Panel C - correlation heatmap
message("Panel C: correlation heatmap...")
cor_mat <- read.csv("results/qc/methylation/methylation_correlation.csv",
                    row.names=1, check.names=FALSE)
df_cor <- melt(as.matrix(cor_mat), varnames=c("Sample1","Sample2"), value.name="r")
sample_order <- c("Scramble_CTRL_1","Scramble_CTRL_2","Scramble_CTRL_3",
                  "ASO_CTRL_1","ASO_CTRL_2","ASO_CTRL_3",
                  "Scramble_VPA_1","Scramble_VPA_2","Scramble_VPA_3",
                  "ASO_VPA_1","ASO_VPA_2","ASO_VPA_3")
existing <- intersect(sample_order, unique(df_cor$Sample1))
df_cor$Sample1 <- factor(df_cor$Sample1, levels=existing)
df_cor$Sample2 <- factor(df_cor$Sample2, levels=rev(existing))

pC <- ggplot(df_cor, aes(x=Sample1, y=Sample2, fill=r)) +
  geom_tile(colour="white", linewidth=0.4) +
  geom_text(aes(label=sprintf("%.3f", r)), size=2.2) +
  scale_fill_gradient2(low="#C0392B", mid="white", high="#1F3A5F",
                       midpoint=0.915, limits=c(0.83, 1.00),
                       name="r", oob=scales::squish) +
  labs(x=NULL, y=NULL, tag="C") +
  theme_classic(base_size=10) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=7),
        axis.text.y=element_text(size=7),
        plot.tag=element_text(face="bold", size=14))

# Panel D - chr1 low-res profile
message("Panel D: chr1 profile...")
meth <- readRDS("results/dmr/meth_pooled_cache.rds")
WIN_SIZE <- 1000000
CHR <- "chr1"
profiles <- lapply(names(COND_COLS), function(cond) {
  m <- meth[[cond]]
  m <- m[as.character(seqnames(m)) == CHR & m$context == "CG" & m$readsN >= 3]
  pos <- start(m)
  bin <- floor(pos/WIN_SIZE)*WIN_SIZE + WIN_SIZE/2
  df <- data.frame(pos=pos, bin=bin, readsM=m$readsM, readsN=m$readsN)
  agg <- df %>% group_by(bin) %>%
    summarise(sumM=sum(readsM), sumN=sum(readsN), .groups="drop") %>%
    filter(sumN >= 10) %>%
    mutate(meth=100*sumM/sumN, condition=cond)
  agg
})
prof_df <- do.call(rbind, profiles)
prof_df$condition <- factor(prof_df$condition, levels=names(COND_COLS))

pD <- ggplot(prof_df, aes(x=bin/1e6, y=meth, colour=condition, linetype=condition)) +
  geom_line(linewidth=0.7, alpha=0.9) +
  annotate("rect", xmin=122, xmax=147, ymin=-Inf, ymax=Inf,
           fill="grey80", alpha=0.4) +
  annotate("text", x=134.5, y=95, label="centromere", size=3, colour="grey50") +
  scale_colour_manual(values=COND_COLS, name=NULL) +
  scale_linetype_manual(values=c(Scramble_CTRL="solid", ASO_CTRL="dashed",
                                 Scramble_VPA="solid", ASO_VPA="dashed"), name=NULL) +
  scale_x_continuous(breaks=seq(0,250,50), labels=function(x) paste0(x," Mb")) +
  scale_y_continuous(limits=c(0,100), breaks=seq(0,100,20),
                     labels=function(x) paste0(x,"%")) +
  labs(x="chr1 position", y="Mean CpG methylation (%)", tag="D") +
  theme_classic(base_size=12) +
  theme(legend.position="right",
        plot.tag=element_text(face="bold", size=14))

# Assemble
message("Assembling panel...")
top <- pA + pB + plot_layout(ncol=2)
bot <- pC + pD + plot_layout(ncol=2)
combined <- top / bot + plot_layout(heights=c(1,1))

ggsave(file.path(OUT, "Fig5.1_QC_panel.pdf"),
       combined, width=16, height=12, device=cairo_pdf)
message("Saved: Fig5.1_QC_panel.pdf")
