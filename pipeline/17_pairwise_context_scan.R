#!/usr/bin/env Rscript
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})
setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")

OUT <- "results/pairwise_context_scan"
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)


message("Loading DMR sets...")
aso_alone  <- readRDS("results/dmr/dmr_ASO_CTRL_vs_Scramble_CTRL.rds")
vpa_alone  <- readRDS("results/dmr/dmr_Scramble_VPA_vs_Scramble_CTRL.rds")
aso_in_vpa <- readRDS("results/dmr/dmr_ASO_VPA_vs_Scramble_VPA.rds")
vpa_in_aso <- readRDS("results/dmr/dmr_ASO_VPA_vs_ASO_CTRL.rds")

cat("ASO alone:", length(aso_alone), "\n")
cat("VPA alone:", length(vpa_alone), "\n")
cat("ASO in VPA:", length(aso_in_vpa), "\n")
cat("VPA in ASO:", length(vpa_in_aso), "\n")


# In ASO_in_VPA but NOT in ASO_alone
message("\nGroup 1: ASO context-dependent loci...")
aso_context <- aso_in_vpa[!overlapsAny(aso_in_vpa, aso_alone)]
cat("ASO context-dependent:", length(aso_context), "\n")

# Rank by absolute methylation difference
aso_context$meth_diff <- abs(aso_context$proportion1 - aso_context$proportion2)
aso_context <- aso_context[order(aso_context$meth_diff, decreasing=TRUE)]


message("Group 2: VPA context-dependent loci...")
vpa_context <- vpa_in_aso[!overlapsAny(vpa_in_aso, vpa_alone)]
cat("VPA context-dependent:", length(vpa_context), "\n")

vpa_context$meth_diff <- abs(vpa_context$proportion1 - vpa_context$proportion2)
vpa_context <- vpa_context[order(vpa_context$meth_diff, decreasing=TRUE)]


message("Annotating top hits...")
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

annotate_top <- function(gr, n=200, label) {
  top <- gr[1:min(n, length(gr))]
  anno <- annotatePeak(top, tssRegion=c(-2000,2000),
                       TxDb=txdb, annoDb="org.Hs.eg.db",
                       verbose=FALSE)
  df <- as.data.frame(anno)
  df$group <- label
  df
}

message("  Annotating ASO context-dependent top 200...")
aso_anno <- annotate_top(aso_context, 200, "ASO_context_dependent")

message("  Annotating VPA context-dependent top 200...")
vpa_anno <- annotate_top(vpa_context, 200, "VPA_context_dependent")

# Save annotated results
write.csv(aso_anno, file.path(OUT, "ASO_context_dependent_top200.csv"),
          row.names=FALSE)
write.csv(vpa_anno, file.path(OUT, "VPA_context_dependent_top200.csv"),
          row.names=FALSE)


cat("\nTop 20 ASO context-dependent genes:\n")
print(head(aso_anno[!is.na(aso_anno$SYMBOL),
  c("seqnames","start","end","meth_diff","SYMBOL","annotation")], 20))

cat("\nTop 20 VPA context-dependent genes:\n")
print(head(vpa_anno[!is.na(vpa_anno$SYMBOL),
  c("seqnames","start","end","meth_diff","SYMBOL","annotation")], 20))


p1 <- ggplot(as.data.frame(aso_context[1:min(500,length(aso_context))]),
             aes(x=meth_diff)) +
  geom_histogram(bins=50, fill="#1F3A5F", colour="white", linewidth=0.1) +
  labs(title="ASO context-dependent loci",
       subtitle=paste0("n=", length(aso_context), " total | top 500 shown"),
       x="|Methylation difference|", y="Count") +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(face="bold"))

p2 <- ggplot(as.data.frame(vpa_context[1:min(500,length(vpa_context))]),
             aes(x=meth_diff)) +
  geom_histogram(bins=50, fill="#C0392B", colour="white", linewidth=0.1) +
  labs(title="VPA context-dependent loci",
       subtitle=paste0("n=", length(vpa_context), " total | top 500 shown"),
       x="|Methylation difference|", y="Count") +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(face="bold"))

combined <- (p1 | p2) +
  plot_annotation(
    title="Context-dependent DMRs — pairwise framework",
    theme=theme(plot.title=element_text(face="bold", size=13))
  )

ggsave(file.path(OUT, "context_dependent_distribution.pdf"),
       combined, width=12, height=5, device=cairo_pdf)

message("\nDone. Outputs in: ", OUT)
