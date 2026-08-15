combo <- read.csv('results/dmr_annotation/ASO_VPA_vs_Scramble_CTRL_annotated.csv')
aso   <- read.csv('results/dmr_annotation/ASO_CTRL_vs_Scramble_CTRL_annotated.csv')

best_per_gene <- function(df) {
  df <- df[!is.na(df$SYMBOL) & df$SYMBOL != "", ]
  df <- df[order(df$pValue), ]
  df[!duplicated(df$SYMBOL),
     c("SYMBOL","seqnames","start","end","proportion1","proportion2",
       "regionType","pValue","annotation","distanceToTSS")]
}

c2 <- best_per_gene(combo); names(c2)[5:10] <- paste0(names(c2)[5:10], "_combo")
a2 <- best_per_gene(aso);   names(a2)[5:10] <- paste0(names(a2)[5:10], "_aso")

m <- merge(c2[, c("SYMBOL","seqnames","start","proportion1_combo","proportion2_combo","regionType_combo",
                   "pValue_combo","annotation_combo","distanceToTSS_combo")],
           a2[, c("SYMBOL","proportion1_aso","proportion2_aso","regionType_aso",
                  "pValue_aso","annotation_aso","distanceToTSS_aso")],
           by="SYMBOL")
m <- m[!duplicated(m$SYMBOL), ]

m$diff_combo <- m$proportion2_combo - m$proportion1_combo
m$diff_aso   <- m$proportion2_aso   - m$proportion1_aso

# reversal: opposite sign, both real DMRs (not NA), reasonably sized effect on both sides
rev <- m[sign(m$diff_combo) != sign(m$diff_aso) &
         abs(m$diff_combo) > 0.05 & abs(m$diff_aso) > 0.05, ]
rev$flip_size <- abs(m$diff_combo[sign(m$diff_combo) != sign(m$diff_aso) &
                                   abs(m$diff_combo) > 0.05 & abs(m$diff_aso) > 0.05]) +
                 abs(m$diff_aso[sign(m$diff_combo) != sign(m$diff_aso) &
                                 abs(m$diff_combo) > 0.05 & abs(m$diff_aso) > 0.05])
rev <- rev[order(-rev$flip_size), ]

cat("=== All genes showing a direction reversal between ASO alone and combination ===\n")
cat("(opposite sign, |diff| > 0.05 in both contrasts)\n\n")
print(head(rev[, c("SYMBOL","seqnames","start","diff_aso","diff_combo",
                   "pValue_aso","pValue_combo","annotation_combo","distanceToTSS_combo")], 25),
      row.names=FALSE)
cat("\nTotal genes with a reversal pattern:", nrow(rev), "\n")
