df <- read.csv('results/dmr_annotation/ASO_VPA_vs_Scramble_CTRL_annotated.csv')
hits <- df[df$SYMBOL=='CACNG5' & !is.na(df$SYMBOL),]
hits <- hits[order(hits$start),]
cat('All CACNG5 DMR windows (combination contrast):\n')
print(hits[,c('seqnames','start','end','proportion1','proportion2','regionType','pValue')], row.names=FALSE)
cat('\nTotal windows:', nrow(hits), '\n')

# Also pull ASO-alone windows for the same gene
df2 <- read.csv('results/dmr_annotation/ASO_CTRL_vs_Scramble_CTRL_annotated.csv')
hits2 <- df2[df2$SYMBOL=='CACNG5' & !is.na(df2$SYMBOL),]
hits2 <- hits2[order(hits2$start),]
cat('\nAll CACNG5 DMR windows (ASO alone contrast):\n')
print(hits2[,c('seqnames','start','end','proportion1','proportion2','regionType','pValue')], row.names=FALSE)
cat('\nTotal windows:', nrow(hits2), '\n')
