.libPaths('~/R/library')
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(dplyr)
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
  'ASO_CTRL_vs_Scramble_CTRL'     = 'ASO vs Scr_CTRL',
  'Scramble_VPA_vs_Scramble_CTRL' = 'VPA vs Scr_CTRL',
  'ASO_VPA_vs_Scramble_CTRL'      = 'ASO+VPA vs Scr_CTRL',
  'ASO_VPA_vs_ASO_CTRL'           = 'ASO+VPA vs ASO',
  'ASO_VPA_vs_Scramble_VPA'       = 'ASO+VPA vs Scr_VPA'
)

df <- do.call(rbind, lapply(CONTRASTS, function(ct) {
  dmr <- readRDS(paste0('results/dmr/dmr_', ct, '.rds'))
  dmr <- dmr[dmr$context == 'CG']
  data.frame(
    contrast  = LABELS[ct],
    direction = c('Hypomethylated', 'Hypermethylated'),
    n         = c(sum(dmr$regionType == 'gain'),
                  sum(dmr$regionType == 'loss'))
  )
}))

df$contrast  <- factor(df$contrast, levels=LABELS)
df$direction <- factor(df$direction, levels=c('Hypermethylated','Hypomethylated'))

# total per contrast for top label
totals <- df %>% group_by(contrast) %>% summarise(total=sum(n))

p <- ggplot(df, aes(x=contrast, y=n, fill=direction)) +
  geom_bar(stat='identity', position='stack', width=0.6) +
  geom_text(aes(label=formatC(n, format='d', big.mark=',')),
            position=position_stack(vjust=0.5),
            size=2.8, colour='white', fontface='bold') +
  geom_text(data=totals,
            aes(x=contrast, y=total, label=formatC(total, format='d', big.mark=',')),
            inherit.aes=FALSE,
            vjust=-0.4, size=3, fontface='bold', colour='grey30') +
  scale_fill_manual(values=c(Hypomethylated='#2166AC',
                              Hypermethylated='#B2182B')) +
  scale_y_log10(labels=scales::comma,
                breaks=c(100, 1000, 10000, 100000, 1000000)) +
  annotation_logticks(sides='l') +
  theme_classic(base_size=12) +
  theme(axis.text.x=element_text(angle=30, hjust=1, size=10),
        plot.title=element_text(face='bold'),
        legend.position='top') +
  labs(title='How many DMRs per contrast and in which direction?',
       subtitle='y-axis is log10 scale',
       x=NULL, y='Number of DMRs (log10)', fill=NULL)

ggsave('results/thesis_figures/slide16_dmr_counts_bar_v2.pdf',
       p, width=10, height=6, device=cairo_pdf)
message('Saved: slide16_dmr_counts_bar_v2.pdf')
