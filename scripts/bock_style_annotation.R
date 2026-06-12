.libPaths(c('~/R/library', .libPaths()))
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(rtracklayer)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(lattice)
  library(AnnotationHub)
})
source("scripts/00_sma_palette.R")
setwd('/data/home/bt25018/sma_epigenomics_pipeline')
OUT <- 'results/dmr_annotation'

# Bock-style annotation panels: CpG context (Panel C) + genomic features (Panel D),
# each with a stacked-bar proportion plot above and a log2(obs/exp) heatmap below.
# reference: Bock et al. 2010 (Nat Biotech). colours follow that figure (Okabe-Ito).

CONTRASTS <- c(
  'ASO_CTRL_vs_Scramble_CTRL',
  'Scramble_VPA_vs_Scramble_CTRL',
  'ASO_VPA_vs_Scramble_CTRL',
  'ASO_VPA_vs_ASO_CTRL',
  'ASO_VPA_vs_Scramble_VPA'
)
LABELS <- c(
  'ASO_CTRL_vs_Scramble_CTRL'     = 'ASO vs Scr_CTRL',
  'Scramble_VPA_vs_Scramble_CTRL' = 'VPA vs Scr_CTRL',
  'ASO_VPA_vs_Scramble_CTRL'      = 'ASO+VPA vs Scr_CTRL',
  'ASO_VPA_vs_ASO_CTRL'           = 'ASO+VPA vs ASO',
  'ASO_VPA_vs_Scramble_VPA'       = 'ASO+VPA vs Scr_VPA'
)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# CpG islands from AnnotationHub; +/- 2kb shores, 2-4 kb shelves
message("loading CpG islands...")
ah  <- AnnotationHub()
cpg <- ah[["AH5086"]]
seqlevelsStyle(cpg) <- "UCSC"
shores <- c(
  GRanges(seqnames(cpg), IRanges(pmax(1,start(cpg)-2000), pmax(1,start(cpg)-1))),
  GRanges(seqnames(cpg), IRanges(end(cpg)+1, end(cpg)+2000))
)
shelves <- c(
  GRanges(seqnames(cpg), IRanges(pmax(1,start(cpg)-4000), pmax(1,start(cpg)-2001))),
  GRanges(seqnames(cpg), IRanges(end(cpg)+2001, end(cpg)+4000))
)
shores  <- shores[width(shores)  > 0]
shelves <- shelves[width(shelves) > 0]

# CpG context classifier: Island > Shore > Shelf > Open_Sea (priority order)
classify_cpg <- function(gr) {
  cls <- rep("Open_Sea", length(gr))
  cls[unique(queryHits(findOverlaps(gr, shelves)))]  <- "Shelf"
  cls[unique(queryHits(findOverlaps(gr, shores)))]   <- "Shore"
  cls[unique(queryHits(findOverlaps(gr, cpg)))]      <- "Island"
  cls
}

# collapse ChIPseeker annotation labels to 7 broad categories
simplify_annot <- function(x) {
  case_when(
    grepl("Promoter", x)          ~ "Promoter",
    grepl("5' UTR|5'UTR", x)      ~ "5' UTR",
    grepl("3' UTR|3'UTR", x)      ~ "3' UTR",
    grepl("Exon", x)              ~ "Exon",
    grepl("Intron", x)            ~ "Intron",
    grepl("Downstream", x)        ~ "Downstream",
    grepl("Distal Intergenic", x) ~ "Intergenic",
    TRUE                          ~ "Other"
  )
}

# whole-genome CpG distribution = denominator for log2(obs/exp)
message("loading genome-wide CpG background...")
cache   <- readRDS('results/dmr/meth_pooled_cache.rds')
all_cpg <- cache[[1]]
all_cpg <- all_cpg[all_cpg$readsN >= 1]
genome_cpg_ctx <- classify_cpg(all_cpg)
genome_cpg_tbl <- table(genome_cpg_ctx) / length(genome_cpg_ctx)

# every 100th CpG for the feature-background table (full set too slow)
message("annotating genome-wide CpGs for feature background...")
all_anno <- annotatePeak(all_cpg[seq(1, length(all_cpg), by=100)],
                         tssRegion=c(-2000,2000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
all_anno_df <- as.data.frame(all_anno)
all_anno_df$feature <- simplify_annot(all_anno_df$annotation)
genome_feat_tbl <- table(all_anno_df$feature) / nrow(all_anno_df)

# fill matrices: rows = Whole Genome + each contrast, cols = categories
message("processing contrasts...")
cpg_mat  <- matrix(0, nrow=length(CONTRASTS)+1,
                   ncol=4,
                   dimnames=list(c("Whole Genome", LABELS),
                                 c("Island","Shore","Shelf","Open_Sea")))
feat_mat <- matrix(0, nrow=length(CONTRASTS)+1,
                   ncol=7,
                   dimnames=list(c("Whole Genome", LABELS),
                                 c("Promoter","5' UTR","Exon","Intron",
                                   "3' UTR","Downstream","Intergenic")))

# row 1 = whole genome reference
for (cat in names(genome_cpg_tbl))
  if (cat %in% colnames(cpg_mat)) cpg_mat["Whole Genome", cat] <- genome_cpg_tbl[cat]
for (cat in names(genome_feat_tbl))
  if (cat %in% colnames(feat_mat)) feat_mat["Whole Genome", cat] <- genome_feat_tbl[cat]

for (ct in CONTRASTS) {
  message("  ", ct)
  lab <- LABELS[ct]
  dmr <- readRDS(paste0('results/dmr/dmr_', ct, '.rds'))
  dmr <- dmr[dmr$context == 'CG']

  # per-DMR CpG context
  ctx <- classify_cpg(dmr)
  ctx_tbl <- table(ctx) / length(ctx)
  for (cat in names(ctx_tbl))
    if (cat %in% colnames(cpg_mat)) cpg_mat[lab, cat] <- ctx_tbl[cat]

  # per-DMR feature annotation
  anno <- annotatePeak(dmr, tssRegion=c(-2000,2000),
                       TxDb=txdb, annoDb="org.Hs.eg.db")
  anno_df <- as.data.frame(anno)
  anno_df$feature <- simplify_annot(anno_df$annotation)
  feat_tbl <- table(anno_df$feature) / nrow(anno_df)
  for (cat in names(feat_tbl))
    if (cat %in% colnames(feat_mat)) feat_mat[lab, cat] <- feat_tbl[cat]
}

# rows -> proportions
cpg_mat  <- cpg_mat  / rowSums(cpg_mat)
feat_mat <- feat_mat / rowSums(feat_mat)

# log2 obs/exp vs whole-genome row (drop the WG row itself before dividing)
log2_cpg  <- log2((cpg_mat[-1,]  + 0.001) /
                  (matrix(cpg_mat["Whole Genome",],
                           nrow=nrow(cpg_mat)-1, ncol=ncol(cpg_mat), byrow=TRUE) + 0.001))
log2_feat <- log2((feat_mat[-1,] + 0.001) /
                  (matrix(feat_mat["Whole Genome",],
                           nrow=nrow(feat_mat)-1, ncol=ncol(feat_mat), byrow=TRUE) + 0.001))

# Okabe-Ito colourblind-safe palette (matches Bock et al. fig)
cbbPalette <- c("#000000","#E69F00","#56B4E9","#009E73",
                "#F0E442","#0072B2","#D55E00","#CC79A7")

# Bock Panel C mapping: Island=black, Shore=light-blue, Shelf=orange, OpenSea=green
cols_cpg <- c(
  Island   = "#000000",
  Shore    = "#56B4E9",
  Shelf    = "#E69F00",
  Open_Sea = "#009E73"
)

# Bock Panel D mapping: Promoter=vermillion, 5'UTR=orange, Exon=green,
# Intron=yellow, 3'UTR=black, Downstream=pink, Intergenic=dark blue
cols_feat <- c(
  "Promoter"   = "#D55E00",
  "5' UTR"     = "#E69F00",
  "Exon"       = "#009E73",
  "Intron"     = "#F0E442",
  "3' UTR"     = "#000000",
  "Downstream" = "#CC79A7",
  "Intergenic" = "#0072B2"
)

# diverging palette for the log2(obs/exp) heatmaps (orange/blue, white at 0)
cols_contrast <- rev(colorRampPalette(c(cbbPalette[7],"white",cbbPalette[6]))(60))

# tighten log2 scale: actual values are <= ~1 in magnitude, +/-3 wastes contrast.
# clip to +/-1.5 and break the colour ramp accordingly.
LOG2_LIM   <- 1.5
custom_at  <- seq(-LOG2_LIM, LOG2_LIM, length.out=61)

log2_cpg[log2_cpg < -LOG2_LIM]   <- -LOG2_LIM
log2_cpg[log2_cpg >  LOG2_LIM]   <-  LOG2_LIM
log2_feat[log2_feat < -LOG2_LIM] <- -LOG2_LIM
log2_feat[log2_feat >  LOG2_LIM] <-  LOG2_LIM

# Panel C: CpG context stacked bar ----
cpg_prop <- as.data.frame(cpg_mat) %>%
  rownames_to_column("group") %>%
  pivot_longer(-group, names_to="context", values_to="pct")
cpg_prop$group   <- factor(cpg_prop$group,
                            levels=c("Whole Genome", LABELS))
cpg_prop$context <- factor(cpg_prop$context,
                            levels=c("Island","Shore","Shelf","Open_Sea"))

p_cpg_bar <- ggplot(cpg_prop, aes(x=group, y=pct, fill=context)) +
  geom_bar(stat="identity", position="fill") +
  scale_fill_manual(values=cols_cpg, name="CpG context") +
  scale_y_continuous(labels=scales::percent_format(1)) +
  theme_classic(base_size=10) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
        legend.position="right") +
  labs(x=NULL, y="Proportion")

ggsave(file.path(OUT, "bock_panelC_cpg_bar.pdf"),
       p_cpg_bar, width=8, height=4, device=cairo_pdf)

# heatmap: enlarge canvas so the cell labels and colour key are readable
log2_cpg_df <- as.data.frame(log2_cpg) %>%
  rownames_to_column("contrast") %>%
  pivot_longer(-contrast, names_to="context", values_to="log2oe")
log2_cpg_df$contrast <- factor(log2_cpg_df$contrast, levels=LABELS)
log2_cpg_df$context  <- factor(log2_cpg_df$context, levels=c("Island","Shore","Shelf","Open_Sea"))
p_cpg_heat <- ggplot(log2_cpg_df, aes(x=context, y=contrast, fill=log2oe)) +
  geom_tile(colour="white", linewidth=0.4) +
  geom_text(aes(label=sprintf("%+.2f", log2oe)), size=3, colour="black") +
  sma_heatmap_scale(limits=c(-1.5,1.5)) +
  theme_minimal(base_size=10) +
  theme(axis.text.x=element_text(angle=30, hjust=1),
        panel.grid=element_blank(),
        plot.title=element_text(face="bold"),
        plot.margin=margin(15,20,15,15)) +
  labs(title="CpG context log2(observed/expected)", x=NULL, y=NULL)
ggsave(file.path(OUT, "bock_panelC_cpg_heatmap.pdf"),
       p_cpg_heat, width=8, height=5, device=cairo_pdf)

# Panel D: genomic feature stacked bar ----
feat_prop <- as.data.frame(feat_mat) %>%
  rownames_to_column("group") %>%
  pivot_longer(-group, names_to="feature", values_to="pct")
feat_prop$group   <- factor(feat_prop$group,
                             levels=c("Whole Genome", LABELS))
feat_prop$feature <- factor(feat_prop$feature,
                             levels=c("Promoter","5' UTR","Exon","Intron",
                                      "3' UTR","Downstream","Intergenic"))

p_feat_bar <- ggplot(feat_prop, aes(x=group, y=pct, fill=feature)) +
  geom_bar(stat="identity", position="fill") +
  scale_fill_manual(values=cols_feat, name="Genomic feature") +
  scale_y_continuous(labels=scales::percent_format(1)) +
  theme_classic(base_size=10) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
        legend.position="right") +
  labs(x=NULL, y="Proportion")

ggsave(file.path(OUT, "bock_panelD_feat_bar.pdf"),
       p_feat_bar, width=8, height=4, device=cairo_pdf)

# wider heatmap (7 categories needs more horizontal space)
log2_feat_df <- as.data.frame(log2_feat) %>%
  rownames_to_column("contrast") %>%
  pivot_longer(-contrast, names_to="feature", values_to="log2oe")
log2_feat_df$contrast <- factor(log2_feat_df$contrast, levels=LABELS)
log2_feat_df$feature  <- factor(log2_feat_df$feature,
  levels=c("Promoter","5' UTR","Exon","Intron","3' UTR","Downstream","Intergenic"))
p_feat_heat <- ggplot(log2_feat_df, aes(x=feature, y=contrast, fill=log2oe)) +
  geom_tile(colour="white", linewidth=0.4) +
  geom_text(aes(label=sprintf("%+.2f", log2oe)), size=3, colour="black") +
  sma_heatmap_scale(limits=c(-1.5,1.5)) +
  theme_minimal(base_size=10) +
  theme(axis.text.x=element_text(angle=30, hjust=1),
        panel.grid=element_blank(),
        plot.title=element_text(face="bold"),
        plot.margin=margin(15,20,15,15)) +
  labs(title="Genomic feature log2(observed/expected)", x=NULL, y=NULL)
ggsave(file.path(OUT, "bock_panelD_feat_heatmap.pdf"),
       p_feat_heat, width=10, height=5, device=cairo_pdf)

message("All Bock-style plots saved to: ", OUT)
