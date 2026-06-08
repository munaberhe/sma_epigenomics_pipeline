#!/usr/bin/env Rscript
.libPaths('~/R/library')
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

setwd('/data/home/bt25018/sma_epigenomics_pipeline')
OUT <- 'results/dmr_annotation'

CONTRASTS <- c(
  'ASO_CTRL_vs_Scramble_CTRL',
  'Scramble_VPA_vs_Scramble_CTRL',
  'ASO_VPA_vs_Scramble_CTRL',
  'ASO_VPA_vs_ASO_CTRL',
  'ASO_VPA_vs_Scramble_VPA'
)

LABELS <- c(
  'ASO_CTRL_vs_Scramble_CTRL'    = 'ASO vs Scr_CTRL',
  'Scramble_VPA_vs_Scramble_CTRL'= 'VPA vs Scr_CTRL',
  'ASO_VPA_vs_Scramble_CTRL'     = 'ASO+VPA vs Scr_CTRL',
  'ASO_VPA_vs_ASO_CTRL'          = 'ASO+VPA vs ASO',
  'ASO_VPA_vs_Scramble_VPA'      = 'ASO vs Scr_VPA'
)

# Simplify annotation categories
simplify_annot <- function(x) {
  case_when(
    grepl("Promoter.*<=1kb", x)  ~ "Promoter (<=1kb)",
    grepl("Promoter.*1-2kb", x)  ~ "Promoter (1-2kb)",
    grepl("Promoter.*2-3kb", x)  ~ "Promoter (2-3kb)",
    grepl("5' UTR", x)           ~ "5' UTR",
    grepl("3' UTR", x)           ~ "3' UTR",
    grepl("Exon", x)             ~ "Exon",
    grepl("Intron", x)           ~ "Intron",
    grepl("Downstream", x)       ~ "Downstream",
    grepl("Distal Intergenic", x)~ "Distal Intergenic",
    TRUE                          ~ "Other"
  )
}

ANNOT_COLS <- c(
  "Promoter (<=1kb)"  = "#B2182B",
  "Promoter (1-2kb)"  = "#D6604D",
  "Promoter (2-3kb)"  = "#F4A582",
  "5' UTR"            = "#FEE090",
  "Exon"              = "#FFFFBF",
  "Intron"            = "#74ADD1",
  "3' UTR"            = "#ABD9E9",
  "Downstream"        = "#4575B4",
  "Distal Intergenic" = "#6B7280",
  "Other"             = "#CCCCCC"
)

all_data <- lapply(CONTRASTS, function(ct) {
  f <- paste0('results/dmr_annotation/', ct, '_annotated.csv')
  if (!file.exists(f)) return(NULL)
  df <- read.csv(f, stringsAsFactors=FALSE)
  df %>%
    mutate(
      annot_simple = simplify_annot(annotation),
      contrast = LABELS[ct],
      direction = ifelse(regionType == "gain", "Hypo", "Hyper")
    ) %>%
    count(contrast, direction, annot_simple)
}) %>% bind_rows()

all_data$contrast <- factor(all_data$contrast, levels=LABELS)
all_data$annot_simple <- factor(all_data$annot_simple, levels=names(ANNOT_COLS))
all_data$direction <- factor(all_data$direction, levels=c("Hyper","Hypo"))

# Plot 1 — stacked bar by count
p1 <- ggplot(all_data, aes(x=contrast, y=n, fill=annot_simple)) +
  geom_bar(stat="identity", position="stack") +
  facet_wrap(~direction, scales="free_y") +
  scale_fill_manual(values=ANNOT_COLS, name="Genomic feature") +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=9),
        legend.position="right",
        plot.title=element_text(face="bold")) +
  labs(title="Genomic distribution of DMRs across contrasts",
       x=NULL, y="Number of DMRs")

# Plot 2 — stacked bar by proportion
p2 <- ggplot(all_data, aes(x=contrast, y=n, fill=annot_simple)) +
  geom_bar(stat="identity", position="fill") +
  facet_wrap(~direction) +
  scale_fill_manual(values=ANNOT_COLS, name="Genomic feature") +
  scale_y_continuous(labels=scales::percent_format(1)) +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=9),
        legend.position="right",
        plot.title=element_text(face="bold")) +
  labs(title="Genomic distribution of DMRs (proportion)",
       x=NULL, y="Proportion of DMRs")

ggsave(file.path(OUT, "DMR_annotation_combined_count.pdf"),
       p1, width=12, height=6, device=cairo_pdf)
ggsave(file.path(OUT, "DMR_annotation_combined_proportion.pdf"),
       p2, width=12, height=6, device=cairo_pdf)

message("Saved: DMR_annotation_combined_count.pdf")
message("Saved: DMR_annotation_combined_proportion.pdf")

for (ct in CONTRASTS) {
  df_ct <- all_data[all_data$contrast == LABELS[ct], ]
  if (nrow(df_ct) == 0) next
  p_ct <- ggplot(df_ct, aes(x=direction, y=n, fill=annot_simple)) +
    geom_bar(stat="identity", position="fill") +
    scale_fill_manual(values=ANNOT_COLS, name="Genomic feature") +
    scale_y_continuous(labels=scales::percent_format(1)) +
    theme_classic(base_size=11) +
    theme(plot.title=element_text(face="bold"),
          legend.position="right") +
    labs(title=paste("Genomic distribution —", LABELS[ct]),
         x=NULL, y="Proportion of DMRs")
  ggsave(file.path(OUT, paste0(ct, "_annotation_bar.pdf")),
         p_ct, width=8, height=5, device=cairo_pdf)
  message("Saved: ", ct, "_annotation_bar.pdf")
}

genes <- unique(hypo$SYMBOL)
eg <- bitr(genes, fromType='SYMBOL', toType='ENTREZID', OrgDb=org.Hs.eg.db)

# KEGG
kegg <- enrichKEGG(gene=eg$ENTREZID, organism='hsa',
                   pvalueCutoff=0.05, qvalueCutoff=0.2)
if (!is.null(kegg) && nrow(kegg) > 0) {
  write.csv(as.data.frame(kegg),
            paste0(OUT, '/ASO_CTRL_vs_Scramble_CTRL_KEGG_hypo.csv'),
            row.names=FALSE)
  p <- dotplot(kegg, showCategory=20) +
    theme_classic(base_size=11) +
    labs(title="KEGG enrichment — ASO hypo DMRs")
  ggsave(paste0(OUT, '/ASO_CTRL_vs_Scramble_CTRL_KEGG_hypo_dotplot.pdf'),
         p, width=10, height=8, device=cairo_pdf)
  message("Saved KEGG hypo dotplot")
} else {
  message("No significant KEGG terms")
}

# Also GO for ASO_VPA vs Scramble_VPA (non-directional, all DMRs)
df2 <- read.csv(paste0(OUT, '/ASO_VPA_vs_Scramble_VPA_annotated.csv'))
all_genes <- unique(df2$SYMBOL[!is.na(df2$SYMBOL)])
eg2 <- bitr(all_genes, fromType='SYMBOL', toType='ENTREZID', OrgDb=org.Hs.eg.db)
go2 <- enrichGO(gene=eg2$ENTREZID, OrgDb=org.Hs.eg.db,
                ont='BP', pvalueCutoff=0.05, qvalueCutoff=0.2,
                readable=TRUE)
if (!is.null(go2) && nrow(go2) > 0) {
  write.csv(as.data.frame(go2),
            paste0(OUT, '/ASO_VPA_vs_Scramble_VPA_GO_BP.csv'),
            row.names=FALSE)
  p2 <- dotplot(go2, showCategory=20) +
    theme_classic(base_size=11) +
    labs(title="GO BP — ASO_VPA vs Scramble_VPA (all DMRs)")
  ggsave(paste0(OUT, '/ASO_VPA_vs_Scramble_VPA_GO_BP_dotplot.pdf'),
         p2, width=10, height=8, device=cairo_pdf)
  message("Saved GO BP dotplot ASO_VPA vs Scramble_VPA")
} else {
  message("No significant GO terms")
}
message("Done")
