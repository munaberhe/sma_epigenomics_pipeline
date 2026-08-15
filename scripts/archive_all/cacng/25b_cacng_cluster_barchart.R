.libPaths('/data/home/bt25018/R/library')
suppressPackageStartupMessages({ library(ggplot2) })
setwd('/data/home/bt25018/sma_epigenomics_pipeline')

OUT_DIR <- 'results/cacng_cluster'
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

# Values confirmed directly from annotated DMR CSVs on 2026-06-19:
#   ASO_VPA_vs_Scramble_CTRL_annotated.csv  (combination)
#   ASO_CTRL_vs_Scramble_CTRL_annotated.csv (ASO alone)
df <- data.frame(
  gene      = c('CACNG1','CACNG4','CACNG5','PRKCA'),
  position  = c('chr17:67.08 Mb','chr17:67.00 Mb','chr17:66.91 Mb','chr17:66.31 Mb'),
  aso_alone = c(NA, NA, 24.4, NA),
  aso_alone_p = c(NA, NA, 0.0008558138, NA),
  combo     = c(27.8, 26.1, -27.6, 27.0),
  combo_p   = c(5.821549e-17, 6.292656e-10, 2.246661e-17, 3.309966e-14),
  category  = c('combo-emergent','combo-emergent','reversal','combo-emergent'),
  stringsAsFactors = FALSE
)

# Order genes by chr17 position, bottom to top in barh
gene_order <- c('PRKCA','CACNG5','CACNG4','CACNG1')
df$gene <- factor(df$gene, levels=gene_order)
df$label <- paste0(df$gene, '\n', df$position)
df$label <- factor(df$label, levels=paste0(gene_order, '\n',
  df$position[match(gene_order, df$gene)]))

# Long format: one row per bar
long <- rbind(
  data.frame(label=df$label, value=df$aso_alone, pval=df$aso_alone_p,
             series='ASO_CTRL vs Scramble_CTRL', category=df$category),
  data.frame(label=df$label, value=df$combo, pval=df$combo_p,
             series='ASO_VPA vs Scramble_CTRL (combination)', category=df$category)
)
long <- long[!is.na(long$value) | long$series=='ASO_CTRL vs Scramble_CTRL', ]
# Keep the "no DMR" placeholder rows for ASO-alone where value is NA -> set near-zero marker
long$is_no_dmr <- is.na(long$value)
long$value[long$is_no_dmr] <- 0.6  # thin placeholder bar

long$plabel <- ifelse(long$is_no_dmr, 'no DMR in ASO alone',
                       sprintf('%+.1f%%  p = %.2e', long$value, long$pval))
long$plabel[!long$is_no_dmr & long$value<0] <-
  sprintf('%.1f%%', long$value[!long$is_no_dmr & long$value<0])

# Background category shading (per gene row)
band_df <- df[, c('label','category')]

p <- ggplot() +
  geom_rect(data=band_df,
            aes(ymin=as.numeric(label)-0.5, ymax=as.numeric(label)+0.5,
                xmin=-Inf, xmax=Inf, fill=category), alpha=0.18) +
  geom_col(data=long,
           aes(x=label, y=value,
               fill=ifelse(is_no_dmr, 'no_dmr',
                           ifelse(series=='ASO_CTRL vs Scramble_CTRL','aso_alone','combo'))),
           width=0.32, position=position_nudge(x=ifelse(long$series=='ASO_CTRL vs Scramble_CTRL', 0.18, -0.18)),
           colour='grey20', orientation='x') +
  coord_flip() +
  geom_hline(yintercept=0, colour='grey30', linewidth=0.5) +
  geom_text(data=long,
            aes(x=label, y=value, label=plabel,
                hjust=ifelse(value>=0, -0.05, 1.05)),
            position=position_nudge(x=ifelse(long$series=='ASO_CTRL vs Scramble_CTRL', 0.18, -0.18)),
            size=3.4, colour='grey20') +
  scale_fill_manual(values=c(
    aso_alone='#5B6B7F', combo='#B0392B', no_dmr='#C9CDD3',
    'combo-emergent'='#F5E6CC', 'reversal'='#F2D2CE'
  ), guide='none') +
  scale_y_continuous(limits=c(-32, 48), breaks=seq(-30,40,10),
                     expand=expansion(mult=c(0.02,0.02))) +
  labs(
    title='CACNG cluster (chr17): combination therapy produces a locus-specific methylation effect',
    subtitle='ASO alone produces no DMR at three of four cluster genes; CACNG5 reverses direction under combination',
    x=NULL, y='Methylation difference vs Scramble_CTRL (%)',
    caption='Source: DMRcaller-B 300 bp strict mode, label-swap null. ASO_VPA vs Scramble_CTRL and ASO_CTRL vs Scramble_CTRL contrasts. n=3 per condition.'
  ) +
  theme_minimal(base_size=13) +
  theme(
    plot.title    = element_text(face='bold', size=15, margin=margin(b=4)),
    plot.subtitle = element_text(size=11, colour='grey30', margin=margin(b=14), face='italic'),
    plot.caption  = element_text(size=8.5, colour='grey45', hjust=0, margin=margin(t=14)),
    axis.text.y   = element_text(size=11, face='bold'),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.margin   = margin(t=20, r=30, b=10, l=10)
  )

ggsave(file.path(OUT_DIR, 'cacng_cluster_combination_effect.pdf'),
       p, width=13, height=7.5, device=cairo_pdf)
message('saved: cacng_cluster_combination_effect.pdf')
