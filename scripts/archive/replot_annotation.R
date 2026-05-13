
.libPaths("~/R/library")
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(ggplot2)

OUT_DIR <- "results/dmr_annotation"
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

contrasts <- c("ASO_VPA_vs_ASO_CTRL", "ASO_VPA_vs_Scramble_CTRL", "VPA_vs_Scramble_CTRL")

for (contrast in contrasts) {
  message("Replotting: ", contrast)
  rds <- readRDS(file.path("results/dmr", contrast, paste0(contrast, "_all_chr.rds")))
  anno <- annotatePeak(rds, tssRegion = c(-3000, 3000),
                       TxDb = txdb, annoDb = "org.Hs.eg.db")

  # Bar chart using ggplot
  anno_df <- as.data.frame(anno)
  feat_counts <- table(anno_df$annotation)
  feat_df <- data.frame(Feature = names(feat_counts), Count = as.numeric(feat_counts))
  feat_df <- feat_df[order(feat_df$Count, decreasing = TRUE), ]

  # Simplify annotation labels — strip transcript/exon detail
  anno@anno$annotation <- gsub(" \\(.*\\)", "", anno@anno$annotation)
  anno_df$annotation   <- gsub(" \\(.*\\)", "", anno_df$annotation)

  feat_counts <- table(anno_df$annotation)
  feat_df <- data.frame(Feature = names(feat_counts), Count = as.numeric(feat_counts))
  feat_df <- feat_df[order(feat_df$Count, decreasing = TRUE), ]

  p_bar <- ggplot(feat_df, aes(x = reorder(Feature, Count), y = Count)) +
    geom_col(fill = "#065A82") +
    coord_flip() +
    labs(title = paste0("Genomic Feature Distribution\n", contrast),
         x = NULL, y = "Number of DMRs") +
    theme_bw(base_size = 11)
  ggsave(file.path(OUT_DIR, paste0(contrast, "_annotation_bar.pdf")),
         p_bar, width = 10, height = 6)

  # TSS distance using ggplot
  anno_df$dist_kb <- anno_df$distanceToTSS / 1000
  p_tss <- ggplot(anno_df, aes(x = dist_kb)) +
    geom_histogram(bins = 30, fill = "#1C7293", colour = "white") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
    labs(title = paste0("Distance to Nearest TSS\n", contrast),
         x = "Distance to TSS (bp)", y = "Count") +
    theme_bw(base_size = 11)
  ggsave(file.path(OUT_DIR, paste0(contrast, "_TSS_distance.pdf")),
         p_tss, width = 10, height = 6)

  message("  Done: ", contrast)
}
message("All replotted")
