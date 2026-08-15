.libPaths('~/R/library')
library(ggplot2)

setwd('/data/home/bt25018/sma_epigenomics_pipeline')

go <- read.csv('results/dss_replicate/DSS_GO_BP.csv')
cat('Terms loaded:', nrow(go), '\n')

# Top 20 neural terms by raw p-value
neural <- go[grep('synap|axon|neuro|cognit|learn|memory|postsynap|presynap',
  go$Description, ignore.case=TRUE), ]
neural <- neural[neural$pvalue < 0.05, ]
neural <- head(neural[order(neural$pvalue), ], 20)
cat('Neural terms p<0.05:', nrow(neural), '\n')

neural$GeneRatio_num <- sapply(neural$GeneRatio, function(x) {
  p <- strsplit(x, '/')[[1]]
  as.numeric(p[1]) / as.numeric(p[2])
})
neural$Description <- factor(neural$Description,
  levels=rev(neural$Description))

p <- ggplot(neural, aes(x=GeneRatio_num, y=Description,
    size=Count, colour=pvalue)) +
  geom_point() +
  scale_colour_gradient(low='#B2182B', high='#4393C3',
    name='p-value') +
  scale_size_continuous(name='Gene count', range=c(2,8)) +
  labs(title='GO BP — DSS replicate-level DMRs',
       subtitle='ASO_CTRL vs Scramble_CTRL  |  nominal p<0.05  |  n=688 DMRs, n=3 replicates',
       x='Gene ratio', y=NULL) +
  theme_bw(base_size=10) +
  theme(plot.title=element_text(face='bold', size=10),
        plot.subtitle=element_text(size=8, colour='grey40'),
        axis.text.y=element_text(size=8))

ggsave('results/dss_replicate/DSS_GO_dotplot.pdf', p, width=10, height=8)
ggsave('results/meeting_plots_radu/DSS_GO_dotplot.pdf', p, width=10, height=8)
cat('Saved\n')
