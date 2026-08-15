.libPaths('/data/home/bt25018/R/library')
suppressPackageStartupMessages({ library(ggplot2) })
setwd('/data/home/bt25018/sma_epigenomics_pipeline')

OUT_DIR <- 'results/genomewide_enhancer'

sum_df <- read.csv(file.path(OUT_DIR, 'genomewide_enhancer_enrichment_summary.csv'))

CONTRASTS <- c('ASO_CTRL_vs_Scramble_CTRL', 'ASO_VPA_vs_Scramble_CTRL',
               'ASO_VPA_vs_ASO_CTRL', 'ASO_VPA_vs_Scramble_VPA',
               'Scramble_VPA_vs_Scramble_CTRL')
sum_df$contrast <- factor(sum_df$contrast, levels=CONTRASTS)
sum_df$sig <- ifelse(sum_df$padj < 0.05, 'p.adj < 0.05', 'ns')

# n enhancers for subtitle -- pull from existing summary if present, else hardcode known value
n_enh <- 555880
N_PERM <- 1000

p_fold <- ggplot(sum_df, aes(x=contrast, y=fold, fill=sig)) +
  geom_col(width=0.65, colour='grey20') +
  geom_hline(yintercept=1, linetype='dashed', colour='grey40') +
  geom_text(aes(label=sprintf('%.2fx\np=%.2g', fold, padj)),
            vjust=ifelse(sum_df$fold > 1, -0.35, 1.35),
            size=3.2) +
  scale_fill_manual(values=c('ns'='grey80', 'p.adj < 0.05'='#C0392B')) +
  scale_y_continuous(expand=expansion(mult=c(0.05, 0.18))) +
  scale_x_discrete(labels=function(x) gsub('_', '\n', x)) +
  theme_classic(base_size=12) +
  theme(
    axis.text.x     = element_text(size=9.5, lineheight=0.85),
    axis.title.y    = element_text(margin=margin(r=10)),
    plot.title      = element_text(face='bold', size=14, margin=margin(b=4)),
    plot.subtitle   = element_text(size=10, colour='grey30', margin=margin(b=14)),
    legend.position = 'top',
    legend.title    = element_blank(),
    plot.margin     = margin(t=20, r=20, b=10, l=10)
  ) +
  labs(title='DMR-enhancer overlap, fold vs chr-matched random',
       subtitle=sprintf(
         'H9 ESC predicted non-promoter enhancers (n=%d) | %d permutations | BH FDR',
         n_enh, N_PERM),
       x=NULL, y='Fold enrichment (observed / expected)')

ggsave(file.path(OUT_DIR, 'genomewide_enhancer_fold_enrichment_v2.pdf'),
       p_fold, width=11, height=7, device=cairo_pdf)
message('saved: genomewide_enhancer_fold_enrichment_v2.pdf')
