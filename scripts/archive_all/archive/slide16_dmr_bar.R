.libPaths('~/R/library')
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
})
setwd('/data/home/bt25018/sma_epigenomics_pipeline')

CONTRASTS <- c(
  'ASO_CTRL_vs_Scramble_CTRL',
  'Scramble_VPA_vs_Scramble_CTRL',
  'ASO_VPA_vs_Scramble_CTRL',
  'ASO_VPA_vs_ASO_CTRL',
  'ASO_VPA_vs_Scramble_VPA'
)

LABELS <- c(
  'ASO_CTRL_vs_Scramble_CTRL'     = 'ASO vs\nScr_CTRL',
  'Scramble_VPA_vs_Scramble_CTRL' = 'VPA vs\nScr_CTRL',
  'ASO_VPA_vs_Scramble_CTRL'      = 'ASO+VPA vs\nScr_CTRL',
  'ASO_VPA_vs_ASO_CTRL'           = 'ASO+VPA vs\nASO',
  'ASO_VPA_vs_Scramble_VPA'       = 'ASO+VPA vs\nScr_VPA'
)

df <- do.call(rbind, lapply(CONTRASTS, function(ct) {
  dmr <- readRDS(paste0('results/dmr/dmr_', ct, '.rds'))
  dmr <- dmr[dmr$context == 'CG']
  data.frame(
    contrast   = LABELS[ct],
    direction  = c('Hypomethylated', 'Hypermethylated'),
    n          = c(sum(dmr$regionType == 'gain'),
                   sum(dmr$regionType == 'loss'))
  )
}))

df$contrast  <- factor(df$contrast, levels=LABELS)
df$direction <- factor(df$direction, levels=c('Hypermethylated','Hypomethylated'))

p <- ggplot(df, aes(x=contrast, y=n, fill=direction)) +
  geom_bar(stat='identity', position='stack') +
  geom_text(aes(label=formatC(n, format='d', big.mark=',')),
            position=position_stack(vjust=0.5),
            size=3, colour='white', fontface='bold') +
  scale_fill_manual(values=c(Hypomethylated='#2166AC', Hypermethylated='#B2182B')) +
  scale_y_continuous(labels=scales::comma) +
  theme_classic(base_size=12) +
  theme(axis.text.x=element_text(size=10),
        plot.title=element_text(face='bold'),
        legend.position='top') +
  labs(title='How many DMRs per contrast and in which direction?',
       x=NULL, y='Number of DMRs', fill=NULL)

ggsave('results/thesis_figures/slide16_dmr_counts_bar.pdf',
       p, width=10, height=6, device=cairo_pdf)
message('Saved: slide16_dmr_counts_bar.pdf')
