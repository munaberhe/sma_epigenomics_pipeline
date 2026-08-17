#!/usr/bin/env Rscript
# 25_metagene_profile.R
# Scaled metagene CpG methylation profile per condition.
# TSS-2kb to TES+2kb, gene body scaled to a common length (60 bins).
# Flanks: 20 bins each at fixed width. Uses chr1 protein-coding genes,
# pooled per-condition methylation (same cache as other QC figures).
# Vectorized: builds all bins for all genes as one GRanges, does one
# findOverlaps call per condition instead of per-gene looping.

.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(ggplot2)
  library(dplyr)
})

setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")
OUT <- "results/thesis_figures"

COND_COLS <- c(
  Scramble_CTRL="#6B7280", ASO_CTRL="#1F3A5F",
  Scramble_VPA="#F0A500",  ASO_VPA="#C0392B"
)

FLANK <- 2000
N_FLANK_BINS <- 20
N_BODY_BINS  <- 60
MIN_GENE_LEN <- 1000
N_TOTAL_BINS <- 2*N_FLANK_BINS + N_BODY_BINS

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_gr <- suppressMessages(genes(txdb, single.strand.genes.only=TRUE))
genes_gr <- genes_gr[seqnames(genes_gr)=="chr1"]
genes_gr <- genes_gr[width(genes_gr) >= MIN_GENE_LEN]
message("chr1 genes used: ", length(genes_gr))

message("Loading pooled methylation cache...")
cache <- readRDS("results/dmr/meth_pooled_cache.rds")
cache <- lapply(cache, function(m) m[seqnames(m)=="chr1" & mcols(m)$context=="CG"])

# build all bins for all genes upfront, vectorized
message("Building bin coordinates for all genes...")
build_all_bins <- function(genes_gr) {
  n_genes <- length(genes_gr)
  all_starts <- vector("list", n_genes)
  all_ends   <- vector("list", n_genes)
  gene_idx   <- vector("list", n_genes)
  bin_idx    <- vector("list", n_genes)

  strand_vec <- as.character(strand(genes_gr))
  gstart_vec <- start(genes_gr)
  gend_vec   <- end(genes_gr)

  for (i in seq_len(n_genes)) {
    s <- strand_vec[i]
    gs <- gstart_vec[i]; ge <- gend_vec[i]
    if (s == "-") { tmp <- gs; gs <- ge; ge <- tmp }

    flank_up   <- seq(gs - FLANK, gs, length.out=N_FLANK_BINS+1)
    body       <- seq(gs, ge, length.out=N_BODY_BINS+1)
    flank_down <- seq(ge, ge+FLANK, length.out=N_FLANK_BINS+1)
    edges <- c(flank_up, body[-1], flank_down[-1])
    if (s == "-") edges <- rev(edges)
    edges <- sort(edges)

    all_starts[[i]] <- head(edges, -1)
    all_ends[[i]]   <- tail(edges, -1)
    gene_idx[[i]]   <- rep(i, N_TOTAL_BINS)
    bin_idx[[i]]    <- seq_len(N_TOTAL_BINS)
  }

  GRanges("chr1",
    IRanges(unlist(all_starts), unlist(all_ends)),
    gene_idx=unlist(gene_idx),
    bin_idx=unlist(bin_idx))
}

all_bins <- build_all_bins(genes_gr)
message("Total bins: ", length(all_bins))

results <- list()
for (cond in names(cache)) {
  message("Processing: ", cond)
  meth <- cache[[cond]]

  ov <- findOverlaps(all_bins, meth, ignore.strand=TRUE)
  message("  overlaps found: ", length(ov))

  # aggregate readsM and readsN per bin using data.table-style vectorized sum
  bin_ids <- queryHits(ov)
  meth_ids <- subjectHits(ov)

  sumM <- tapply(mcols(meth)$readsM[meth_ids], bin_ids, sum)
  sumN <- tapply(mcols(meth)$readsN[meth_ids], bin_ids, sum)

  # build per-bin proportion, indexed to all_bins length
  prop <- rep(NA_real_, length(all_bins))
  valid_idx <- as.integer(names(sumN))
  keep <- sumN[valid_idx - valid_idx + 1] >= 5  # min 5 reads per bin
  prop_vals <- ifelse(sumN >= 5, sumM/sumN, NA)
  prop[valid_idx] <- prop_vals

  # average across genes per bin_idx position (1..N_TOTAL_BINS)
  bin_pos <- mcols(all_bins)$bin_idx
  mean_profile <- tapply(prop, bin_pos, mean, na.rm=TRUE)
  mean_profile <- mean_profile[as.character(seq_len(N_TOTAL_BINS))]

  results[[cond]] <- as.numeric(mean_profile)
}

df <- do.call(rbind, lapply(names(results), function(cond) {
  data.frame(condition=cond, bin=seq_len(N_TOTAL_BINS), meth=results[[cond]]*100)
}))
df$condition <- factor(df$condition, levels=names(COND_COLS))

write.csv(df, file.path(OUT, "metagene_profile_data.csv"), row.names=FALSE)

# Relabel conditions for the figure (script filenames are not figure labels)
LABELS <- c(
  Scramble_CTRL = "Scramble + vehicle",
  ASO_CTRL      = "ASO1 + vehicle",
  Scramble_VPA  = "Scramble + VPA",
  ASO_VPA       = "ASO1 + VPA"
)
df$condition_label <- factor(LABELS[as.character(df$condition)],
                             levels=LABELS[names(COND_COLS)])
names(COND_COLS) <- LABELS[names(COND_COLS)]

# Quantify TSS dip depth per condition: gene-body plateau minus TSS minimum
tss_bin <- N_FLANK_BINS + 1
body_plateau_bins <- (N_FLANK_BINS + 10):(N_FLANK_BINS + N_BODY_BINS - 10)
dip_stats <- df %>%
  group_by(condition) %>%
  summarise(
    body_plateau = mean(meth[bin %in% body_plateau_bins], na.rm=TRUE),
    tss_min = min(meth[bin %in% (tss_bin-2):(tss_bin+2)], na.rm=TRUE),
    dip_depth = body_plateau - tss_min,
    .groups="drop"
  )
write.csv(dip_stats, file.path(OUT, "metagene_TSS_dip_depth.csv"), row.names=FALSE)
message("\nTSS dip depth (gene-body plateau minus TSS minimum):")
print(dip_stats)

y_min <- floor(min(df$meth, na.rm=TRUE)/5)*5
y_max <- ceiling(max(df$meth, na.rm=TRUE)/5)*5

p <- ggplot(df, aes(x=bin, y=meth, colour=condition_label)) +
  geom_line(linewidth=1) +
  geom_vline(xintercept=N_FLANK_BINS+0.5, linetype="dashed", colour="grey40") +
  geom_vline(xintercept=N_FLANK_BINS+N_BODY_BINS+0.5, linetype="dashed", colour="grey40") +
  annotate("text", x=N_FLANK_BINS+0.5, y=y_max*0.99,
           label="TSS", size=3.5, hjust=1.1) +
  annotate("text", x=N_FLANK_BINS+N_BODY_BINS+0.5, y=y_max*0.99,
           label="TES", size=3.5, hjust=-0.1) +
  scale_colour_manual(values=COND_COLS, name=NULL) +
  scale_x_continuous(
    breaks=c(1, N_FLANK_BINS+1, N_FLANK_BINS+N_BODY_BINS+1, N_TOTAL_BINS),
    labels=c("-2kb", "TSS", "TES", "+2kb")
  ) +
  scale_y_continuous(breaks=seq(y_min, y_max, 5),
                     expand=expansion(mult=c(0.02, 0.05))) +
  labs(x="Gene body is length-scaled between TSS and TES; flanks are fixed-width (2kb)",
       y="Mean CpG methylation (%)") +
  theme_classic(base_size=13) +
  theme(legend.position="top",
        axis.title.x=element_text(size=9, colour="grey40"))

ggsave(file.path(OUT, "Fig_metagene_profile.pdf"),
       p, width=8, height=5.5, device=cairo_pdf)
message("Saved: Fig_metagene_profile.pdf")
