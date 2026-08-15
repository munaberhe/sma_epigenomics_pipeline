.libPaths(c('~/R/library', .libPaths()))
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(annotatr)
  library(tibble)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
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

# Load CpG islands from annotatr (hg38) - fixes AH5086 hg19 bug
library(annotatr)
  library(tibble)
cpg_anns <- build_annotations(genome="hg38",
  annotations=c("hg38_cpg_islands","hg38_cpg_shores","hg38_cpg_shelves"))
cpg_islands <- cpg_anns[cpg_anns$type=="hg38_cpg_islands"]
shores  <- cpg_anns[cpg_anns$type=="hg38_cpg_shores"]
shelves <- cpg_anns[cpg_anns$type=="hg38_cpg_shelves"]
shore_size <- 2000
shelf_size <- 2000

classify_cpg <- function(gr) {
  is_island <- length(findOverlaps(gr, cpg_islands)) > 0
  is_shore  <- length(findOverlaps(gr, shores))      > 0
  is_shelf  <- length(findOverlaps(gr, shelves))     > 0
  ifelse(is_island, "Island",
  ifelse(is_shore,  "Shore",
  ifelse(is_shelf,  "Shelf", "Open_Sea")))
}

# Genome-wide CpG context counts (expected)
message("computing genome-wide CpG context...")
genome_cpgs <- readRDS('results/dmr/meth_pooled_cache.rds')[[1]]
genome_cpgs <- genome_cpgs[genome_cpgs$readsN >= 1]

cpg_context <- data.frame(
  island = length(subsetByOverlaps(genome_cpgs, cpg_islands)),
  shore  = length(subsetByOverlaps(genome_cpgs, shores)),
  shelf  = length(subsetByOverlaps(genome_cpgs, shelves))
)
cpg_context$open_sea <- length(genome_cpgs) - rowSums(cpg_context)
expected_pct <- cpg_context / sum(cpg_context)
message("Genome-wide CpG context: ", paste(round(expected_pct*100,1), collapse=", "))

# Per-contrast CpG context
results <- lapply(CONTRASTS, function(ct) {
  message("processing: ", ct)
  dmr <- readRDS(paste0('results/dmr/dmr_', ct, '.rds'))
  dmr <- dmr[dmr$context == 'CG']

  classify_batch <- function(gr) {
    ov_island <- findOverlaps(gr, cpg_islands)
    ov_shore  <- findOverlaps(gr, shores)
    ov_shelf  <- findOverlaps(gr, shelves)
    cls <- rep("Open_Sea", length(gr))
    cls[unique(queryHits(ov_shelf))]  <- "Shelf"
    cls[unique(queryHits(ov_shore))]  <- "Shore"
    cls[unique(queryHits(ov_island))] <- "Island"
    cls
  }

  cls <- classify_batch(dmr)
  data.frame(
    contrast  = LABELS[ct],
    context   = cls,
    direction = ifelse(dmr$regionType == 'gain', 'Hypo', 'Hyper')
  )
})

all_df <- do.call(rbind, results)
all_df$contrast <- factor(all_df$contrast, levels=LABELS)
all_df$context  <- factor(all_df$context,
                           levels=c("Island","Shore","Shelf","Open_Sea"))

write.csv(all_df, file.path(OUT, 'DMR_cpg_context.csv'), row.names=FALSE)

# Panel A: 100% stacked bar
p_bar <- ggplot(all_df, aes(x=contrast, fill=context)) +
  geom_bar(position='fill') +
  facet_wrap(~direction) +
  scale_fill_manual(values=c(Island="#B2182B", Shore="#EF8A62",
                              Shelf="#FDDBC7", Open_Sea="#D1E5F0"),
                    name="CpG context") +
  scale_y_continuous(labels=scales::percent_format(1)) +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=9),
        plot.title=element_text(face='bold'),
        legend.position='right') +
  labs(x=NULL, y='Proportion of DMRs')

ggsave(file.path(OUT, 'DMR_cpg_context_bar.pdf'),
       p_bar, width=12, height=5, device=cairo_pdf)
message("saved: DMR_cpg_context_bar.pdf")

# Panel B: log2 obs/exp heatmap
obs_mat <- all_df %>%
  group_by(contrast, context) %>%
  summarise(n=n(), .groups='drop') %>%
  group_by(contrast) %>%
  mutate(pct=n/sum(n)) %>%
  select(contrast, context, pct) %>%
  pivot_wider(names_from=context, values_from=pct, values_fill=0) %>%
  column_to_rownames('contrast') %>%
  as.matrix()

exp_vec <- c(Island=as.numeric(expected_pct[1,'island']),
             Shore =as.numeric(expected_pct[1,'shore']),
             Shelf =as.numeric(expected_pct[1,'shelf']),
             Open_Sea=as.numeric(expected_pct[1,'open_sea']))

log2_mat <- log2((obs_mat + 0.001) /
                 (matrix(exp_vec[colnames(obs_mat)],
                         nrow=nrow(obs_mat),
                         ncol=ncol(obs_mat), byrow=TRUE) + 0.001))

col_fun <- colorRamp2(c(-3, 0, 3),
                      c("#2166AC", "white", "#B2182B"))

pdf(file.path(OUT, 'DMR_cpg_context_heatmap.pdf'), width=6, height=4)
draw(Heatmap(log2_mat,
             name="log2\nobs/exp",
             col=col_fun,
             cluster_rows=FALSE,
             cluster_columns=FALSE,
             cell_fun=function(j, i, x, y, width, height, fill) {
               grid.text(sprintf("%.2f", log2_mat[i,j]),
                         x, y, gp=gpar(fontsize=9))
             },
             row_names_side="left",
             column_names_rot=45,
             heatmap_legend_param=list(title="log2\nobs/exp")))
dev.off()
message("saved: DMR_cpg_context_heatmap.pdf")
message("done.")
