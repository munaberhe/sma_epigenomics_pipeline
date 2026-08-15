df <- read.csv('results/dmr_annotation/ASO_VPA_vs_Scramble_CTRL_annotated.csv')
df <- df[!is.na(df$SYMBOL) & df$SYMBOL != "", ]

# hypo only: regionType "gain" means cond1 (ASO_VPA) is lower than cond2 (Scramble_CTRL)
df <- df[df$regionType == "gain", ]

# p value threshold
df <- df[df$pValue < 0.001, ]

# methylation difference (proportion2 - proportion1), keep best (most extreme) row per gene
df$methDiff <- df$proportion2 - df$proportion1
df <- df[order(-df$methDiff), ]
df_unique <- df[!duplicated(df$SYMBOL), ]

cat("Top 10 hypomethylated genes, ASO_VPA vs Scramble_CTRL, sorted by methylation difference, p < 0.001:\n")
print(head(df_unique[, c("SYMBOL","seqnames","start","end","methDiff","pValue","annotation")], 10),
      row.names=FALSE)
