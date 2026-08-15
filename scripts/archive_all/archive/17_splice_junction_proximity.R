#!/usr/bin/env Rscript
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
})

# Tests whether ASO-specific DMRs are closer to exon-intron splice junctions
# than a matched random background.
# Splice sites from Ensembl GRCh38.109 GTF (protein-coding only):
# data/reference/splice_sites_hg38_protein_coding.bed

DMR_FILE   <- "results/dmr/dmr_ASO_specific.rds"
SPLICE_BED <- "data/reference/splice_sites_hg38_protein_coding.bed"
OUT_DIR    <- "results/splice_junction"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

KEEP_CHRS <- paste0("chr", 1:22)
N_BG      <- 1000
set.seed(42)

message("loading ASO-specific DMRs...")
dmrs <- readRDS(DMR_FILE)
dmrs <- dmrs[as.character(seqnames(dmrs)) %in% KEEP_CHRS]
message("  n = ", length(dmrs))

message("loading splice sites...")
bed <- read.table(SPLICE_BED, header=FALSE, sep="\t",
                  col.names=c("chr","start","end"))
bed <- bed[bed$chr %in% KEEP_CHRS, ]
splice_sites <- GRanges(bed$chr, IRanges(bed$end, bed$end))
message("  splice sites: ", length(splice_sites))

# distance from each DMR to nearest splice junction
message("calculating distances...")
dmr_dist      <- distanceToNearest(dmrs, splice_sites)
dmr_distances <- mcols(dmr_dist)$distance
message("  median: ", median(dmr_distances), " bp")
message("  within 300bp: ", sum(dmr_distances <= 300),
        " (", round(100*mean(dmr_distances <= 300), 1), "%)")

# matched background regions
message("generating background (n=", N_BG, ")...")
chr_sizes <- c(
  chr1=248956422, chr2=242193529, chr3=198295559, chr4=190214555,
  chr5=181538259, chr6=170805979, chr7=159345973, chr8=145138636,
  chr9=138394717, chr10=133797422, chr11=135086622, chr12=133275309,
  chr13=114364328, chr14=107043718, chr15=101991189, chr16=90338345,
  chr17=83257441,  chr18=80373285,  chr19=58617616,  chr20=64444167,
  chr21=46709983,  chr22=50818468)
chr_counts <- table(as.character(seqnames(dmrs)))
chr_probs  <- chr_counts / sum(chr_counts)
dmr_widths <- width(dmrs)
bg_list  <- list()
attempts <- 0
while (length(bg_list) < N_BG && attempts < N_BG*20) {
  attempts  <- attempts + 1
  chr  <- sample(names(chr_probs), 1, prob=chr_probs)
  w    <- sample(dmr_widths, 1)
  maxs <- chr_sizes[chr] - w
  if (maxs < 1) next
  s <- sample(1:maxs, 1)
  candidate <- GRanges(chr, IRanges(s, s+w-1))
  if (length(findOverlaps(candidate, dmrs)) == 0)
    bg_list[[length(bg_list)+1]] <- candidate
}
bg_regions   <- do.call(c, bg_list)
bg_dist      <- distanceToNearest(bg_regions, splice_sites)
bg_distances <- mcols(bg_dist)$distance
message("  background median: ", median(bg_distances), " bp")

# Wilcoxon test
wtest <- wilcox.test(dmr_distances, bg_distances, alternative="less")
message("Wilcoxon p = ", signif(wtest$p.value, 3))

# save results
results_df <- data.frame(
  type     = c(rep("ASO DMRs",   length(dmr_distances)),
               rep("Background", length(bg_distances))),
  distance = c(dmr_distances, bg_distances)
)
write.csv(results_df, file.path(OUT_DIR, "splice_junction_distances.csv"), row.names=FALSE)

summary_df <- data.frame(
  group            = c("ASO_specific_DMRs", "Background"),
  n                = c(length(dmr_distances), length(bg_distances)),
  median_dist_bp   = c(median(dmr_distances), median(bg_distances)),
  pct_within_300bp = c(round(100*mean(dmr_distances<=300), 1),
                        round(100*mean(bg_distances<=300),  1)),
  wilcox_p         = c(signif(wtest$p.value, 3), NA)
)
print(summary_df, row.names=FALSE)
write.csv(summary_df, file.path(OUT_DIR, "splice_junction_summary.csv"), row.names=FALSE)

# density plot
p <- ggplot(results_df, aes(x=distance, fill=type, colour=type)) +
  geom_density(alpha=0.4, linewidth=0.8) +
  scale_x_log10(labels=scales::comma) +
  scale_fill_manual(values=c("ASO DMRs"="#1B4F8A", "Background"="#cccccc")) +
  scale_colour_manual(values=c("ASO DMRs"="#1B4F8A", "Background"="#888888")) +
  geom_vline(xintercept=300, linetype="dashed", colour="grey40", linewidth=0.5) +
  annotate("text", x=350, y=Inf, label="300bp bin",
           hjust=0, vjust=1.5, size=3, colour="grey40") +
  labs(title="Distance from ASO-specific DMRs to nearest splice junction",
       subtitle=sprintf("n=%d DMRs vs %d background | Wilcoxon p=%s | DMR median=%d bp | BG median=%d bp",
                        length(dmr_distances), N_BG,
                        signif(wtest$p.value,3),
                        median(dmr_distances), median(bg_distances)),
       x="distance to nearest splice junction (bp, log scale)",
       y="density", fill=NULL, colour=NULL) +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(face="bold"), legend.position="top")
ggsave(file.path(OUT_DIR, "splice_junction_distance_density.pdf"), p, width=10, height=6)
message("saved: splice_junction_distance_density.pdf")
message("\ndone. outputs in: ", OUT_DIR)
