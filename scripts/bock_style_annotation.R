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
  'ASO_CTRL_vs_Scramble_CTRL'     = 'ASO vs Scr_CTRL',
  'Scramble_VPA_vs_Scramble_CTRL' = 'VPA vs Scr_CTRL',
  'ASO_VPA_vs_Scramble_CTRL'      = 'ASO+VPA vs Scr_CTRL',
  'ASO_VPA_vs_ASO_CTRL'           = 'ASO+VPA vs ASO',
  'ASO_VPA_vs_Scramble_VPA'       = 'ASO+VPA vs Scr_VPA'
)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Load CpG islands
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

classify_cpg <- function(gr) {
  cls <- rep("Open_Sea", length(gr))
  cls[unique(queryHits(findOverlaps(gr, shelves)))]  <- "Shelf"
  cls[unique(queryHits(findOverlaps(gr, shores)))]   <- "Shore"
  cls[unique(queryHits(findOverlaps(gr, cpg)))]      <- "Island"
  cls
}

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

# Load genome-wide CpGs for background
message("loading genome-wide CpG background...")
cache   <- readRDS('results/dmr/meth_pooled_cache.rds')
all_cpg <- cache[[1]]
all_cpg <- all_cpg[all_cpg$readsN >= 1]
genome_cpg_ctx <- classify_cpg(all_cpg)
genome_cpg_tbl <- table(genome_cpg_ctx) / length(genome_cpg_ctx)

# Load genome-wide annotation for background
message("annotating genome-wide CpGs for feature background...")
all_anno <- annotatePeak(all_cpg[seq(1, length(all_cpg), by=100)],
                         tssRegion=c(-2000,2000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
all_anno_df <- as.data.frame(all_anno)
all_anno_df$feature <- simplify_annot(all_anno_df$annotation)
genome_feat_tbl <- table(all_anno_df$feature) / nrow(all_anno_df)

# Process each contrast
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

# Fill genome row
for (cat in names(genome_cpg_tbl))
  if (cat %in% colnames(cpg_mat)) cpg_mat["Whole Genome", cat] <- genome_cpg_tbl[cat]
for (cat in names(genome_feat_tbl))
  if (cat %in% colnames(feat_mat)) feat_mat["Whole Genome", cat] <- genome_feat_tbl[cat]

for (ct in CONTRASTS) {
  message("  ", ct)
  lab <- LABELS[ct]
  dmr <- readRDS(paste0('results/dmr/dmr_', ct, '.rds'))
  dmr <- dmr[dmr$context == 'CG']

  # CpG context
  ctx <- classify_cpg(dmr)
  ctx_tbl <- table(ctx) / length(ctx)
  for (cat in names(ctx_tbl))
    if (cat %in% colnames(cpg_mat)) cpg_mat[lab, cat] <- ctx_tbl[cat]

  # Genomic features
  anno <- annotatePeak(dmr, tssRegion=c(-2000,2000),
                       TxDb=txdb, annoDb="org.Hs.eg.db")
  anno_df <- as.data.frame(anno)
  anno_df$feature <- simplify_annot(anno_df$annotation)
  feat_tbl <- table(anno_df$feature) / nrow(anno_df)
  for (cat in names(feat_tbl))
    if (cat %in% colnames(feat_mat)) feat_mat[lab, cat] <- feat_tbl[cat]
}

# Normalise rows to sum to 1
cpg_mat  <- cpg_mat  / rowSums(cpg_mat)
feat_mat <- feat_mat / rowSums(feat_mat)

# log2 obs/exp (vs Whole Genome row)
log2_cpg  <- log2((cpg_mat[-1,]  + 0.001) /
                  (matrix(cpg_mat["Whole Genome",],
                           nrow=nrow(cpg_mat)-1, ncol=ncol(cpg_mat), byrow=TRUE) + 0.001))
log2_feat <- log2((feat_mat[-1,] + 0.001) /
                  (matrix(feat_mat["Whole Genome",],
                           nrow=nrow(feat_mat)-1, ncol=ncol(feat_mat), byrow=TRUE) + 0.001))

# Colour palettes
cbbPalette <- c("#000000","#E69F00","#56B4E9","#009E73",
                "#F0E442","#0072B2","#D55E00","#CC79A7")
cols_cpg   <- cbbPalette[c(6,4,5,1)]
cols_feat  <- cbbPalette[c(7,2,3,4,5,6,1)]
cols_contrast <- rev(colorRampPalette(c(cbbPalette[7],"white",cbbPalette[6]))(60))
custom_at <- seq(-3, 3, by=(8/60))

# Clamp
log2_cpg[log2_cpg < -3]   <- -3
log2_cpg[log2_cpg > 3]    <- 3
log2_feat[log2_feat < -3] <- -3
log2_feat[log2_feat > 3]  <- 3

# Plot Panel C — CpG island
cpg_prop <- as.data.frame(cpg_mat) %>%
  rownames_to_column("group") %>%
  pivot_longer(-group, names_to="context", values_to="pct")
cpg_prop$group   <- factor(cpg_prop$group,
                            levels=c("Whole Genome", LABELS))
cpg_prop$context <- factor(cpg_prop$context,
                            levels=c("Island","Shore","Shelf","Open_Sea"))

p_cpg_bar <- ggplot(cpg_prop, aes(x=group, y=pct, fill=context)) +
  geom_bar(stat="identity", position="fill") +
  scale_fill_manual(values=setNames(cols_cpg, c("Island","Shore","Shelf","Open_Sea")),
                    name="CpG context") +
  scale_y_continuous(labels=scales::percent_format(1)) +
  theme_classic(base_size=10) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
        legend.position="right") +
  labs(x=NULL, y="Proportion")

ggsave(file.path(OUT, "bock_panelC_cpg_bar.pdf"),
       p_cpg_bar, width=8, height=4, device=cairo_pdf)

pdf(file.path(OUT, "bock_panelC_cpg_heatmap.pdf"), width=6, height=3)
print(levelplot(t(log2_cpg),
                at=custom_at, col.regions=cols_contrast,
                main="CpG context log2(obs/exp)",
                xlab="", ylab="",
                scales=list(x=list(rot=45))))
dev.off()

# Plot Panel D — genomic features
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
  scale_fill_manual(values=setNames(cols_feat,
                                    c("Promoter","5' UTR","Exon","Intron",
                                      "3' UTR","Downstream","Intergenic")),
                    name="Genomic feature") +
  scale_y_continuous(labels=scales::percent_format(1)) +
  theme_classic(base_size=10) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
        legend.position="right") +
  labs(x=NULL, y="Proportion")

ggsave(file.path(OUT, "bock_panelD_feat_bar.pdf"),
       p_feat_bar, width=8, height=4, device=cairo_pdf)

pdf(file.path(OUT, "bock_panelD_feat_heatmap.pdf"), width=7, height=3)
print(levelplot(t(log2_feat),
                at=custom_at, col.regions=cols_contrast,
                main="Genomic feature log2(obs/exp)",
                xlab="", ylab="",
                scales=list(x=list(rot=45))))
dev.off()

message("All Bock-style plots saved to: ", OUT)
