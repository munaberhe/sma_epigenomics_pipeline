combo <- read.csv('results/dmr_annotation/ASO_VPA_vs_Scramble_CTRL_annotated.csv')
aso   <- read.csv('results/dmr_annotation/ASO_CTRL_vs_Scramble_CTRL_annotated.csv')

best_per_gene <- function(df) {
  df <- df[!is.na(df$SYMBOL) & df$SYMBOL != "", ]
  df <- df[order(df$pValue), ]
  df[!duplicated(df$SYMBOL),
     c("SYMBOL","seqnames","start","end","proportion1","proportion2",
       "pValue","annotation","distanceToTSS")]
}
c2 <- best_per_gene(combo); names(c2)[5:9] <- paste0(names(c2)[5:9], "_combo")
a2 <- best_per_gene(aso);   names(a2)[5:9] <- paste0(names(a2)[5:9], "_aso")

m <- merge(c2, a2[, c("SYMBOL","proportion1_aso","proportion2_aso","pValue_aso")], by="SYMBOL")
m <- m[!duplicated(m$SYMBOL), ]
m$diff_combo <- m$proportion2_combo - m$proportion1_combo
m$diff_aso   <- m$proportion2_aso   - m$proportion1_aso

rev <- m[sign(m$diff_combo) != sign(m$diff_aso) &
         abs(m$diff_combo) > 0.05 & abs(m$diff_aso) > 0.05, ]
# filter out uncharacterized loci, pseudogenes, linc/mir/antisense/divergent transcripts
rev <- rev[!grepl("^LOC|^MIR|^LINC|P[0-9]*$|-AS[0-9]*$|-DT$", rev$SYMBOL), ]
rev$flip_size <- abs(rev$diff_combo) + abs(rev$diff_aso)
rev <- rev[order(-rev$flip_size), ]

cat("Named-gene reversals (excluding LOC/MIR/LINC/pseudogenes):\n")
print(head(rev[,c("SYMBOL","seqnames","start","diff_aso","diff_combo",
                  "annotation_combo","distanceToTSS_combo")], 20), row.names=FALSE)
cat("\nTotal named-gene reversals:", nrow(rev), "\n")
