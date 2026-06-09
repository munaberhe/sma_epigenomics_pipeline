.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(ggplot2)
  library(dplyr)
})
setwd("/data/home/bt25018/sma_epigenomics_pipeline")
OUT <- "results/dmr_annotation"



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
    grepl("Promoter", x)         ~ "Promoter",
    grepl("5' UTR", x)           ~ "5' UTR",
    grepl("3' UTR", x)           ~ "3' UTR",
    grepl("Exon", x)             ~ "Exon",
    grepl("Intron", x)           ~ "Intron",
    grepl("Downstream", x)       ~ "Downstream",
    grepl("Distal Intergenic", x)~ "Intergenic",
    TRUE                          ~ "Other"
  )
}

ANNOT_COLS <- c(
  "Promoter"   = "#B2182B",
  "5' UTR"     = "#F4A582",
  "Exon"       = "#FFFFBF",
  "Intron"     = "#74ADD1",
  "3' UTR"     = "#ABD9E9",
  "Downstream" = "#4575B4",
  "Intergenic" = "#6B7280",
  "Other"      = "#CCCCCC"
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

# pre-compute proportions to avoid stat="fill" issue
prop_data <- all_data %>%
  group_by(contrast, direction) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Plot 2 - stacked bar by proportion with percentage labels
p2 <- ggplot(prop_data, aes(x=contrast, y=prop, fill=annot_simple)) +
  geom_bar(stat="identity", position="stack") +
  geom_text(aes(label=ifelse(prop > 0.04, paste0(round(prop*100,1),"%"), ""),
                group=contrast),
            position=position_stack(vjust=0.5),
            size=2.5, colour="white", fontface="bold") +
  facet_wrap(~direction) +
  scale_fill_manual(values=ANNOT_COLS, name="Genomic feature") +
  scale_y_continuous(labels=scales::percent_format(1),
                     expand=expansion(mult=c(0, 0.02))) +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=9),
        legend.position="right",
        plot.title=element_text(face="bold")) +
  labs(title="Genomic distribution of DMRs (proportion)",
       x=NULL, y="Proportion of DMRs")

ggsave(file.path(OUT, "DMR_annotation_combined_count.pdf"),
       p1, width=12, height=6, device=cairo_pdf)
ggsave(file.path(OUT, "DMR_annotation_combined_proportion.pdf"),
       p2, width=16, height=6, device=cairo_pdf)

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

message("done. outputs in: ", OUT)
