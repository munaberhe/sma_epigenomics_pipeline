#!/usr/bin/env Rscript
.libPaths('~/R/library')
suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
  library(ggplot2)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
})

setwd('/data/home/bt25018/sma_epigenomics_pipeline')
OUT <- 'results/tss_metaplot'
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

COLS <- c(ASO_CTRL="#1B4F8A", ASO_VPA="#B2182B",
          Scramble_VPA="#F0A500", Scramble_CTRL="#6B7280")

WIN <- 3000
BIN <- 100

message("Loading cache...")
cache <- readRDS('results/dmr/meth_pooled_cache.rds')

message("Getting TSS coordinates (chr1 only)...")
txdb  <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes <- genes(txdb)
genes <- genes[seqnames(genes) == "chr1" & width(genes) > 1000]
message("  genes: ", length(genes))

# get TSS as single points
tss <- ifelse(as.character(strand(genes)) == "+",
              start(genes), end(genes))
tss_gr <- GRanges("chr1", IRanges(tss, tss), strand=strand(genes))

# create bin windows around all TSS at once
bins <- seq(-WIN, WIN, by=BIN)
bin_mids <- bins[-length(bins)] + BIN/2

message("Building TSS windows...")
all_windows <- do.call(c, lapply(seq_along(tss_gr), function(i) {
  t <- tss[i]
  st <- as.character(strand(tss_gr[i]))
  if (st == "+") {
    GRanges("chr1", IRanges(t + bins[-length(bins)], t + bins[-1] - 1))
  } else {
    GRanges("chr1", IRanges(t - bins[-1] + 1, t - bins[-length(bins)]),
            bin_idx = rev(seq_along(bin_mids)))
  }
}))
mcols(all_windows)$gene_idx <- rep(seq_along(tss_gr), each=length(bin_mids))
mcols(all_windows)$bin_idx  <- rep(seq_along(bin_mids), length(tss_gr))

message("Computing profiles using findOverlaps...")
compute_profile <- function(meth_chr1) {
  ov <- findOverlaps(all_windows, meth_chr1)
  M  <- mcols(meth_chr1)$readsM[subjectHits(ov)]
  N  <- mcols(meth_chr1)$readsN[subjectHits(ov)]
  bi <- mcols(all_windows)$bin_idx[queryHits(ov)]
  keep <- N >= 3
  tapply(M[keep], bi[keep], sum) / tapply(N[keep], bi[keep], sum)
}

cache_chr1 <- lapply(cache, function(m) m[seqnames(m) == "chr1"])

profiles <- list()
for (cond in names(cache_chr1)) {
  message("  computing: ", cond)
  profiles[[cond]] <- compute_profile(cache_chr1[[cond]])
}

prof_df <- do.call(rbind, lapply(names(profiles), function(cond) {
  vals <- as.numeric(profiles[[cond]])
  data.frame(pos=bin_mids[as.integer(names(profiles[[cond]]))],
             meth=vals, condition=cond)
}))
prof_df$condition <- factor(prof_df$condition,
  levels=c("ASO_CTRL","Scramble_CTRL","ASO_VPA","Scramble_VPA"))

p <- ggplot(prof_df, aes(x=pos, y=meth, colour=condition)) +
  geom_line(linewidth=0.9, na.rm=TRUE) +
  geom_vline(xintercept=0, linetype="dashed",
             colour="grey40", linewidth=0.5) +
  annotate("text", x=200, y=max(prof_df$meth, na.rm=TRUE)*0.97,
           label="TSS", size=3.5, colour="grey40", hjust=0) +
  scale_colour_manual(values=COLS) +
  scale_x_continuous(breaks=seq(-3000,3000,1000),
                     labels=function(x) paste0(x/1000,"kb")) +
  scale_y_continuous(labels=scales::percent_format(1)) +
  labs(title="Average CpG methylation around TSS (chr1 genes)",
       x="Distance from TSS", y="CpG methylation", colour=NULL) +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(face="bold"), legend.position="right")

ggsave(file.path(OUT, "TSS_metaplot.png"), p,
       width=9, height=5, dpi=150, bg="white")
message("saved as PNG")
message("saved: TSS_metaplot.pdf")
