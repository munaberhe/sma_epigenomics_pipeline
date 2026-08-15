.libPaths('~/R/library')
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(patchwork)
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
  dmr <- dmr[dmr$context == 'CG' & dmr$cytosinesCount >= 4]
  rbind(
    data.frame(contrast=LABELS[ct], direction='Hypomethylated',
               meth=dmr$proportion1[dmr$regionType=='gain']),
    data.frame(contrast=LABELS[ct], direction='Hypomethylated',
               meth=dmr$proportion2[dmr$regionType=='gain']),
    data.frame(contrast=LABELS[ct], direction='Hypermethylated',
               meth=dmr$proportion1[dmr$regionType=='loss']),
    data.frame(contrast=LABELS[ct], direction='Hypermethylated',
               meth=dmr$proportion2[dmr$regionType=='loss'])
  )
}))

df$contrast  <- factor(df$contrast, levels=LABELS)
df$direction <- factor(df$direction, levels=c('Hypomethylated','Hypermethylated'))

p <- ggplot(df, aes(x=contrast, y=meth, fill=direction)) +
  geom_violin(trim=FALSE, alpha=0.8, linewidth=0.3) +
  geom_boxplot(width=0.07, fill='white',
               outlier.size=0.2, outlier.alpha=0.2) +
  facet_wrap(~direction, ncol=2) +
  scale_fill_manual(values=c(Hypomethylated='#2166AC',
                              Hypermethylated='#B2182B')) +
  scale_y_continuous(labels=scales::percent_format(1), limits=c(0,1)) +
  theme_classic(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=9),
        plot.title=element_text(face='bold'),
        legend.position='none',
        strip.text=element_text(face='bold', size=11)) +
  labs(title='CpG methylation at DMR loci — hypo vs hyper separated',
       x=NULL, y='CpG methylation proportion')

ggsave('results/thesis_figures/slide18_violin_hypo_hyper.pdf',
       p, width=14, height=6, device=cairo_pdf)
message('Saved: slide18_violin_hypo_hyper.pdf')
