#!/usr/bin/env Rscript
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(ggplot2)
  library(patchwork)
})

setwd("/data/scratch/bt25018/sma_epigenomics_pipeline")

# H3K9me2 signal enrichment at ASO-associated DMRs.
# Data: Marasco et al. 2022 Cell (GSE167762), HEK293T cells.
# Tests whether ASO-linked DMRs overlap with H3K9me2-marked heterochromatin.

BW_FILES <- c(
  CTR_R1 = "data/external/GSE167762_H3K9me2/GSM6063702_CTRvsInp_CTR_R1.bw",
  CTR_R2 = "data/external/GSE167762_H3K9me2/GSM6063706_CTRvsInp_CTR_R2.bw",
  ASO_R1 = "data/external/GSE167762_H3K9me2/GSM6063703_ASOvsInp_ASO_R1.bw",
  ASO_R2 = "data/external/GSE167762_H3K9me2/GSM6063707_ASOvsInp_ASO_R2.bw"
)

DMR_DIR <- "results/dmr"
OUT_DIR <- "results/h3k9me2_overlap"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

KEEP_CHRS <- paste0("chr", 1:22)

# load DMRs
message("loading DMRs...")
dmr_aso      <- readRDS(file.path(DMR_DIR, "dmr_ASO_CTRL_vs_Scramble_CTRL.rds"))
dmr_specific <- readRDS(file.path(DMR_DIR, "dmr_ASO_specific.rds"))
dmr_aso      <- dmr_aso[as.character(seqnames(dmr_aso)) %in% KEEP_CHRS]
dmr_specific <- dmr_specific[as.character(seqnames(dmr_specific)) %in% KEEP_CHRS]
message("  ASO DMRs: ", length(dmr_aso))
message("  ASO-specific DMRs: ", length(dmr_specific))

# mean bigWig signal over a set of ranges
bw_mean_over_ranges <- function(bw_file, gr) {
  message("  importing: ", basename(bw_file))
  bw <- import(bw_file, which=gr, as="NumericList")
  vapply(bw, function(x) if(length(x)==0) NA_real_ else mean(x, na.rm=TRUE),
         numeric(1))
}

# matched background regions with same width/chromosome distribution
make_background <- function(query, n_bg=NULL, seed=42) {
  if (is.null(n_bg)) n_bg <- length(query)
  set.seed(seed)
  chr_sizes <- c(
    chr1=248956422, chr2=242193529, chr3=198295559, chr4=190214555,
    chr5=181538259, chr6=170805979, chr7=159345973, chr8=145138636,
    chr9=138394717, chr10=133797422, chr11=135086622, chr12=133275309,
    chr13=114364328, chr14=107043718, chr15=101991189, chr16=90338345,
    chr17=83257441, chr18=80373285, chr19=58617616, chr20=64444167,
    chr21=46709983, chr22=50818468)
  chr_counts <- table(as.character(seqnames(query)))
  chr_counts <- chr_counts[names(chr_counts) %in% names(chr_sizes)]
  chr_probs  <- chr_counts / sum(chr_counts)
  query_widths <- width(query)
  bg_list <- list()
  attempts <- 0
  while (length(bg_list) < n_bg && attempts < n_bg*10) {
    attempts <- attempts + 1
    chr  <- sample(names(chr_probs), 1, prob=chr_probs)
    w    <- sample(query_widths, 1)
    maxs <- chr_sizes[chr] - w
    if (maxs < 1) next
    s <- sample(1:maxs, 1)
    candidate <- GRanges(chr, IRanges(s, s+w-1))
    if (length(findOverlaps(candidate, query)) == 0)
      bg_list[[length(bg_list)+1]] <- candidate
  }
  do.call(c, bg_list)
}

compute_signal <- function(gr, label) {
  message("\ncomputing signal over ", label, " (n=", length(gr), ")...")
  data.frame(
    CTR_R1 = bw_mean_over_ranges(BW_FILES["CTR_R1"], gr),
    CTR_R2 = bw_mean_over_ranges(BW_FILES["CTR_R2"], gr),
    ASO_R1 = bw_mean_over_ranges(BW_FILES["ASO_R1"], gr),
    ASO_R2 = bw_mean_over_ranges(BW_FILES["ASO_R2"], gr),
    group  = label
  )
}

# subsample ASO DMRs for speed
set.seed(1)
dmr_aso_sub <- dmr_aso[sample(length(dmr_aso), min(5000, length(dmr_aso)))]

message("making background regions...")
bg_aso      <- make_background(dmr_aso_sub)
bg_specific <- make_background(dmr_specific, n_bg=length(dmr_specific)*5)

sig_aso      <- compute_signal(dmr_aso_sub,  "ASO_DMR")
sig_bg_aso   <- compute_signal(bg_aso,       "Background")
sig_specific <- compute_signal(dmr_specific, "ASO_specific_DMR")
sig_bg_spec  <- compute_signal(bg_specific,  "Background")

# add mean columns
for (df_name in c("sig_aso","sig_bg_aso","sig_specific","sig_bg_spec")) {
  df          <- get(df_name)
  df$mean_CTR <- rowMeans(df[,c("CTR_R1","CTR_R2")], na.rm=TRUE)
  df$mean_ASO <- rowMeans(df[,c("ASO_R1","ASO_R2")], na.rm=TRUE)
  assign(df_name, df)
}

# Wilcoxon tests
test1 <- wilcox.test(sig_aso$mean_CTR,      sig_bg_aso$mean_CTR,  alternative="greater")
test2 <- wilcox.test(sig_specific$mean_CTR, sig_bg_spec$mean_CTR, alternative="greater")
message("Wilcoxon p-values:")
message("  ASO DMRs vs BG: p = ", signif(test1$p.value, 3))
message("  ASO-specific vs BG: p = ", signif(test2$p.value, 3))

summary_df <- data.frame(
  contrast         = c("ASO_DMRs", "ASO_specific_DMRs"),
  n                = c(nrow(sig_aso), nrow(sig_specific)),
  median_DMR_CTR   = c(median(sig_aso$mean_CTR,      na.rm=TRUE),
                        median(sig_specific$mean_CTR, na.rm=TRUE)),
  median_BG_CTR    = c(median(sig_bg_aso$mean_CTR,   na.rm=TRUE),
                        median(sig_bg_spec$mean_CTR,  na.rm=TRUE)),
  wilcox_p         = c(signif(test1$p.value,3), signif(test2$p.value,3))
)
print(summary_df, row.names=FALSE)
write.table(summary_df, file.path(OUT_DIR, "h3k9me2_signal_summary.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# plots
plot_df <- rbind(
  data.frame(signal=sig_aso$mean_CTR,      group="ASO DMRs",   set="All ASO DMRs"),
  data.frame(signal=sig_bg_aso$mean_CTR,   group="Background", set="All ASO DMRs"),
  data.frame(signal=sig_specific$mean_CTR, group="ASO DMRs",   set="ASO-specific DMRs"),
  data.frame(signal=sig_bg_spec$mean_CTR,  group="Background", set="ASO-specific DMRs")
)
saveRDS(plot_df, file.path(OUT_DIR, "h3k9me2_plot_df.rds"))
plot_df$group <- factor(plot_df$group, levels=c("Background","ASO DMRs"))

plot_df2 <- rbind(
  data.frame(signal=sig_aso$mean_CTR,      condition="CTR", set="All ASO DMRs"),
  data.frame(signal=sig_aso$mean_ASO,      condition="ASO", set="All ASO DMRs"),
  data.frame(signal=sig_specific$mean_CTR, condition="CTR", set="ASO-specific DMRs"),
  data.frame(signal=sig_specific$mean_ASO, condition="ASO", set="ASO-specific DMRs")
)
plot_df2$condition <- factor(plot_df2$condition, levels=c("CTR","ASO"))

make_box <- function(df, x, fill_col, title, subtitle, y_lab) {
  ggplot(df, aes_string(x=x, y="signal", fill=x)) +
    geom_boxplot(outlier.size=0.3, outlier.alpha=0.3) +
    facet_wrap(~set, scales="free_y") +
    scale_fill_manual(values=fill_col) +
    theme_classic(base_size=11) +
    theme(legend.position="none", strip.text=element_text(face="bold")) +
    labs(title=title, subtitle=subtitle, x=NULL, y=y_lab)
}

p1 <- make_box(plot_df, "group",
               c("Background"="#cccccc","ASO DMRs"="#1B4F8A"),
               "H3K9me2 at ASO DMRs vs background",
               "CTR ChIP-seq (Marasco et al. 2022)", "Mean H3K9me2 signal")
p2 <- make_box(plot_df2, "condition",
               c("CTR"="#6B7280","ASO"="#1B4F8A"),
               "H3K9me2 at DMRs: CTR vs ASO",
               "tests whether ASO increases H3K9me2 at DMR loci",
               "Mean H3K9me2 signal")

combined <- p1 / p2 +
  plot_annotation(title="H3K9me2 enrichment at WGBS DMR loci",
                  subtitle="GSE167762, HEK293T",
                  theme=theme(plot.title=element_text(face="bold", size=13)))
ggsave(file.path(OUT_DIR, "h3k9me2_signal_boxplot.pdf"), combined, width=10, height=9)

# SMN2 locus spot check
message("\nSMN2 locus H3K9me2 signal...")
smn2    <- GRanges("chr5", IRanges(70049638, 70078522))
smn2_df <- data.frame(
  condition = names(BW_FILES),
  signal    = sapply(BW_FILES, function(f) bw_mean_over_ranges(f, smn2))
)
print(smn2_df, row.names=FALSE)
write.table(smn2_df, file.path(OUT_DIR, "smn2_h3k9me2_signal.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

message("\ndone. outputs in: ", OUT_DIR)
