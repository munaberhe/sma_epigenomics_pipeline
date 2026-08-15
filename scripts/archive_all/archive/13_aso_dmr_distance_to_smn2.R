suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
})
.libPaths(c("~/R/library", .libPaths()))
setwd("/data/scratch/bt25018/sma_epigenomics_pipeline")

OUT_DIR <- "results/dmr_overlap"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

smn2 <- GRanges("chr5", IRanges(70049638, 70078522))

message("Loading ASO-specific DMRs...")
dmr_specific <- readRDS("results/dmr/dmr_ASO_specific.rds")
message("  n = ", length(dmr_specific))

chr <- as.character(seqnames(dmr_specific))
is_chr5 <- chr == "chr5"

message("  On chr5: ", sum(is_chr5))
message("  On other chromosomes: ", sum(!is_chr5))
message("  Trans proportion: ", round(100 * mean(!is_chr5), 1), "%")

chr5_dmrs <- dmr_specific[is_chr5]
if (length(chr5_dmrs) > 0) {
  dist_to_smn2 <- distance(chr5_dmrs, smn2)
  message("  Chr5 DMRs distance to SMN2 (bp):")
  message("    min: ",    min(dist_to_smn2, na.rm=TRUE))
  message("    median: ", median(dist_to_smn2, na.rm=TRUE))
  message("    max: ",    max(dist_to_smn2, na.rm=TRUE))
  message("    within 1Mb of SMN2: ", sum(dist_to_smn2 <= 1e6, na.rm=TRUE))
}

dist_df <- data.frame(
  dmr_id   = seq_along(dmr_specific),
  chr      = chr,
  start    = start(dmr_specific),
  end      = end(dmr_specific),
  width    = width(dmr_specific),
  location = ifelse(is_chr5, "chr5 (cis)", "other chromosome (trans)"),
  dist_to_smn2_bp = ifelse(
    is_chr5,
    distance(dmr_specific, smn2),
    NA_integer_
  )
)

write.table(dist_df,
            file.path(OUT_DIR, "aso_specific_dmr_distances.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

chr_counts <- as.data.frame(table(chr), stringsAsFactors=FALSE)
colnames(chr_counts) <- c("chr", "n")
chr_counts$chr <- factor(chr_counts$chr,
  levels=c(paste0("chr", 1:22), "chrX", "chrY"))
chr_counts <- chr_counts[!is.na(chr_counts$chr), ]

p <- ggplot(chr_counts, aes(x=chr, y=n,
       fill=ifelse(chr=="chr5", "chr5", "other"))) +
  geom_col() +
  scale_fill_manual(values=c("chr5"="#C0392B", "other"="#95A5A6"),
                    guide="none") +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs(title="Chromosomal distribution of ASO-specific DMRs (n=151)",
       subtitle="Red = chr5 (cis to SMN2); grey = other chromosomes (trans)",
       x=NULL, y="Number of DMRs")

ggsave(file.path(OUT_DIR, "aso_specific_dmr_chr_distribution.pdf"),
       p, width=10, height=5)

message("Saved: aso_specific_dmr_distances.tsv")
message("Saved: aso_specific_dmr_chr_distribution.pdf")
message("Done.")
