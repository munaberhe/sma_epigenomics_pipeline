.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(GenomicRanges); library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene); library(org.Hs.eg.db)
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
  library(annotatr)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")
OUT <- "results/dmr_annotation"
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# 3 main contrasts only
CONTRASTS <- c("ASO_CTRL_vs_Scramble_CTRL",
               "Scramble_VPA_vs_Scramble_CTRL",
               "ASO_VPA_vs_Scramble_CTRL")
LABELS <- c(ASO_CTRL_vs_Scramble_CTRL = "ASO only",
            Scramble_VPA_vs_Scramble_CTRL = "VPA only",
            ASO_VPA_vs_Scramble_CTRL = "ASO+VPA")

# Radu palette
NAVY <- "#1F3A5F"; RED <- "#C0392B"
OBS_COLS <- c(Promoter="#D55E00", `5' UTR`="#E69F00", Exon="#009E73",
              Intron="#F0E442", `3' UTR`="#882255",
              Downstream="#CC79A7", Intergenic="#0072B2")

simplify_annot <- function(x) {
  dplyr::case_when(
    grepl("Promoter", x)         ~ "Promoter",
    grepl("5' UTR|5'UTR", x)     ~ "5' UTR",
    grepl("3' UTR|3'UTR", x)     ~ "3' UTR",
    grepl("Exon", x)             ~ "Exon",
    grepl("Intron", x)           ~ "Intron",
    grepl("Downstream", x)       ~ "Downstream",
    grepl("Distal Intergenic", x)~ "Intergenic",
    TRUE                         ~ "Other"
  )
}

FEAT_LEVELS <- c("Promoter","5' UTR","Exon","Intron","3' UTR","Downstream","Intergenic")

# genome background
message("building background...")
gr_bg <- genes(txdb)
anno_bg <- annotatePeak(gr_bg, tssRegion=c(-2000,2000),
                        TxDb=txdb, annoDb="org.Hs.eg.db", verbose=FALSE)
bg_df <- as.data.frame(anno_bg)
bg_df$feature <- simplify_annot(bg_df$annotation)
bg_tbl <- table(bg_df$feature) / nrow(bg_df)

# per-contrast
feat_mat <- matrix(0, nrow=length(CONTRASTS)+1, ncol=length(FEAT_LEVELS),
  dimnames=list(c("Genome", LABELS), FEAT_LEVELS))
for (cat in names(bg_tbl))
  if (cat %in% FEAT_LEVELS) feat_mat["Genome", cat] <- bg_tbl[cat]

for (ct in CONTRASTS) {
  message("  ", ct)
  csv <- file.path(OUT, paste0(ct, "_annotated.csv"))
  if (!file.exists(csv)) next
  df <- read.csv(csv, stringsAsFactors=FALSE)
  df$feature <- simplify_annot(df$annotation)
  tbl <- table(df$feature) / nrow(df)
  lab <- LABELS[ct]
  for (cat in names(tbl))
    if (cat %in% FEAT_LEVELS) feat_mat[lab, cat] <- tbl[cat]
}

feat_mat <- feat_mat / rowSums(feat_mat)

# log2 obs/exp
log2_mat <- log2((feat_mat[-1,] + 0.001) /
  (matrix(feat_mat["Genome",], nrow=nrow(feat_mat)-1,
          ncol=ncol(feat_mat), byrow=TRUE) + 0.001))
log2_mat[log2_mat < -2] <- -2
log2_mat[log2_mat >  2] <-  2

# Panel A: stacked bar
bar_df <- as.data.frame(feat_mat) %>%
  tibble::rownames_to_column("group") %>%
  pivot_longer(-group, names_to="feature", values_to="pct")
bar_df$group   <- factor(bar_df$group, levels=c("Genome", LABELS))
bar_df$feature <- factor(bar_df$feature, levels=FEAT_LEVELS)

pA <- ggplot(bar_df, aes(x=group, y=pct, fill=feature)) +
  geom_col(position="fill", width=0.75) +
  scale_fill_manual(values=OBS_COLS, name="Genomic feature") +
  scale_y_continuous(labels=scales::percent_format(1)) +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=35, hjust=1, size=9),
        legend.position="right",
        plot.title=element_text(face="bold")) +
  labs(title="Genomic feature distribution of DMRs",
       x=NULL, y="Proportion")

# Panel B: log2 obs/exp heatmap strip
heat_df <- as.data.frame(log2_mat) %>%
  tibble::rownames_to_column("contrast") %>%
  pivot_longer(-contrast, names_to="feature", values_to="log2oe")
heat_df$contrast <- factor(heat_df$contrast, levels=LABELS)
heat_df$feature  <- factor(heat_df$feature, levels=FEAT_LEVELS)

pB <- ggplot(heat_df, aes(x=feature, y=contrast, fill=log2oe)) +
  geom_tile(colour="white", linewidth=0.4) +
  geom_text(aes(label=sprintf("%+.2f", log2oe)), size=2.8, colour="black") +
  scale_fill_gradient2(low="#D55E00", mid="white", high="#0072B2",
                       midpoint=0, limits=c(-2,2), name="log2\nobs/exp") +
  theme_minimal(base_size=10) +
  theme(axis.text.x=element_text(angle=30, hjust=1),
        panel.grid=element_blank(),
        plot.title=element_text(face="bold")) +
  labs(title="log2(observed/expected) vs genome background",
       x=NULL, y=NULL)

combined <- pA / pB +
  plot_layout(heights=c(2,1)) +
  plot_annotation(
    title="DMR genomic feature enrichment — 3 main contrasts",
    subtitle="Top: feature proportions | Bottom: log2(obs/exp) vs genome background",
    theme=theme(plot.title=element_text(face="bold", size=13)))

ggsave(file.path(OUT, "bock_fig1d_style_v2.pdf"),
       combined, width=10, height=9, device=cairo_pdf)
message("saved: bock_fig1d_style_v2.pdf")
