df <- read.csv('results/dmr_annotation/ASO_VPA_vs_Scramble_CTRL_annotated.csv')
df <- df[!is.na(df$SYMBOL) & df$SYMBOL != '',]
df$mb_bin <- paste0(df$seqnames, ':', floor(df$start/1e6))
tab <- table(df$mb_bin)
tab <- sort(tab, decreasing=TRUE)
cat('Top 25 densest 1Mb bins by DMR count (combination, ASO_VPA vs Scramble_CTRL):\n')
print(head(tab, 25))

cat('\n\nFor each top bin, which genes are present:\n')
top_bins <- names(head(tab, 10))
for (b in top_bins) {
  parts <- strsplit(b, ':')[[1]]
  chr <- parts[1]; mb <- as.numeric(parts[2])
  sub <- df[df$seqnames == chr & df$start >= mb*1e6 & df$start < (mb+1)*1e6, ]
  genes <- sort(unique(sub$SYMBOL))
  cat(sprintf('%s (%d DMRs): %s\n', b, nrow(sub), paste(genes, collapse=', ')))
}
