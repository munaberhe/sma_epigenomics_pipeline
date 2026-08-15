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
    contrast = LABELS[ct],
    hypo     = sum(dmr$regionType == 'gain'),
    hyper    = -sum(dmr$regionType == 'loss')
  )
}))

df$contrast <- factor(df$contrast, levels=LABELS)

df_long <- rbind(
  data.frame(contrast=df$contrast, direction='Hypomethylated', n=df$hypo),
  data.frame(contrast=df$contrast, direction='Hypermethylated', n=df$hyper)
)
df_long$direction <- factor(df_long$direction,
                             levels=c('Hypomethylated','Hypermethylated'))

p <- ggplot(df_long, aes(x=contrast, y=n, fill=direction)) +
  geom_bar(stat='identity', width=0.6) +
  geom_hline(yintercept=0, linewidth=0.5, colour='grey30') +
  geom_text(aes(label=formatC(abs(n), format='d', big.mark=','),
                vjust=ifelse(n >= 0, -0.3, 1.2)),
            size=3, fontface='bold', colour='grey20') +
  scale_fill_manual(values=c(Hypomethylated='#2166AC',
                              Hypermethylated='#B2182B')) +
  scale_y_continuous(labels=function(x) formatC(abs(x), format='d', big.mark=',')) +
  facet_wrap(~contrast, scales='free_y', nrow=1) +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text=element_text(face='bold', size=9),
        legend.position='top',
        plot.title=element_text(face='bold')) +
  labs(x=NULL, y='Number of DMRs\n(hypo above, hyper below)',
       fill=NULL)

ggsave('results/thesis_figures/slide16_dmr_diverging.pdf',
       p, width=14, height=6, device=cairo_pdf)
message('Saved: slide16_dmr_diverging.pdf')
