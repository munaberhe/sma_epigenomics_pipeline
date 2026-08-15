.libPaths('~/R/library')
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

OUT <- 'results/dmr_annotation'

contrasts <- c('ASO_CTRL_vs_Scramble_CTRL',
               'ASO_VPA_vs_Scramble_CTRL',
               'Scramble_VPA_vs_Scramble_CTRL')

feature_cols <- c(
  'Promoter (<=1kb)'     = '#B2182B',
  'Promoter (1-2kb)'     = '#D6604D',
  'Promoter (2-3kb)'     = '#F4A582',
  "5' UTR"               = '#92C5DE',
  "3' UTR"               = '#4393C3',
  '1st Exon'             = '#2166AC',
  'Other Exon'           = '#053061',
  '1st Intron'           = '#4DAC26',
  'Other Intron'         = '#B8E186',
  'Downstream (<=300)'   = '#762A83',
  'Distal Intergenic'    = '#C2A5CF'
)

for (contrast in contrasts) {
  cat('Processing:', contrast, '\n')
  df <- read.csv(file.path(OUT, paste0(contrast, '_annotated.csv')))

  # Clean annotation to broad category
  df$feature <- gsub(' \\(ENST.*', '', df$annotation)
  df$feature <- gsub(' \\(.*', '', df$feature)
  df$feature <- trimws(df$feature)

  # Count and compute percentages
  counts <- df %>%
    group_by(feature) %>%
    summarise(n=n(), .groups='drop') %>%
    mutate(pct = round(n/sum(n)*100, 1)) %>%
    arrange(desc(pct))

  # Order factor
  counts$feature <- factor(counts$feature, levels=counts$feature)

  cairo_pdf(file.path(OUT, paste0(contrast, '_annotation_bar.pdf')),
            width=10, height=6)
  p <- ggplot(counts, aes(x=reorder(feature, pct), y=pct,
                           fill=feature)) +
    geom_bar(stat='identity', width=0.7) +
    geom_text(aes(label=paste0(pct,'%')),
              hjust=-0.1, size=3.5) +
    scale_fill_manual(values=feature_cols, na.value='#999999') +
    coord_flip() +
    theme_minimal(base_size=12) +
    theme(legend.position='none',
          panel.grid.major.y=element_blank()) +
    expand_limits(y=max(counts$pct)*1.15) +
    labs(title=paste0('Genomic Feature Distribution — ', contrast),
         x='', y='Percentage of DMRs (%)')
  print(p)
  dev.off()
  cat(' done\n')
}
cat('Done.\n')
