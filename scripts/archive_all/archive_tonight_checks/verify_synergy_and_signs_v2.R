combo <- read.csv('results/dmr_annotation/ASO_VPA_vs_Scramble_CTRL_annotated.csv')
aso   <- read.csv('results/dmr_annotation/ASO_CTRL_vs_Scramble_CTRL_annotated.csv')
vpa   <- read.csv('results/dmr_annotation/Scramble_VPA_vs_Scramble_CTRL_annotated.csv')

best_per_gene <- function(df) {
  df <- df[!is.na(df$SYMBOL) & df$SYMBOL != "",]
  df <- df[order(df$pValue),]
  df[!duplicated(df$SYMBOL),
     c("SYMBOL","proportion1","proportion2","regionType","pValue")]
}
c2 <- best_per_gene(combo); names(c2)[2:5] <- paste0(names(c2)[2:5], "_combo")
a2 <- best_per_gene(aso);   names(a2)[2:5] <- paste0(names(a2)[2:5], "_aso")

# enforce one row per gene (best_per_gene should already do this, but guard anyway)
c2 <- c2[!duplicated(c2$SYMBOL), ]
a2 <- a2[!duplicated(a2$SYMBOL), ]

m <- merge(c2, a2, by="SYMBOL")
m <- m[!duplicated(m$SYMBOL), ]
m$diff_combo <- m$proportion2_combo - m$proportion1_combo
m$diff_aso   <- m$proportion2_aso   - m$proportion1_aso

amp <- m[sign(m$diff_combo) == sign(m$diff_aso) &
         abs(m$diff_combo) > abs(m$diff_aso),]
amp$amplification <- abs(m$diff_combo) - abs(m$diff_aso)
amp <- amp[order(-amp$amplification),]
cat("=== Top 20 combination-amplified genes ===\n")
print(head(amp[, c("SYMBOL","diff_aso","diff_combo","amplification",
                   "pValue_aso","pValue_combo")], 20), row.names=FALSE)

blunt <- m[sign(m$diff_combo) == sign(m$diff_aso) &
           abs(m$diff_combo) < abs(m$diff_aso),]
blunt$blunting <- abs(m$diff_aso) - abs(m$diff_combo)
blunt <- blunt[order(-blunt$blunting),]
cat("\n=== Top 20 combination-blunted genes ===\n")
print(head(blunt[, c("SYMBOL","diff_aso","diff_combo","blunting",
                     "pValue_aso","pValue_combo")], 20), row.names=FALSE)

cat("\n=== Locked synergy gene check (amplified list) ===\n")
locked_amp <- c("CD38","TSHZ1","HUWE1","ROCK1")
for (g in locked_amp) {
  rank <- which(amp$SYMBOL == g)
  cat(g, ": rank", if(length(rank)>0) rank else "NOT FOUND", "out of", nrow(amp), "\n")
}

cat("\n=== Locked synergy gene check (blunted list) ===\n")
locked_blunt <- c("GFRA2","SLC32A1")
for (g in locked_blunt) {
  rank <- which(blunt$SYMBOL == g)
  cat(g, ": rank", if(length(rank)>0) rank else "NOT FOUND", "out of", nrow(blunt), "\n")
}

cat("\n=== VPA-alone DMR check (is CACNG combination-specific, or just VPA?) ===\n")
v <- vpa[!is.na(vpa$SYMBOL) & vpa$SYMBOL != "",]
for (g in c("CACNG1","CACNG4","CACNG5","PRKCA",
            "ARID1B","LMO7","SEMA3C","LINC00391","POU2F1")) {
  hits <- v[v$SYMBOL == g, c("start","end","proportion1","proportion2","pValue")]
  cat(g, ": ", nrow(hits), " VPA-alone DMRs",
      if (nrow(hits) > 0) paste0(" (best p=", signif(min(hits$pValue),2), ")") else "",
      "\n", sep="")
}

cat("\n=== CACNG5 sign convention double-check ===\n")
cacng5_aso <- aso[aso$SYMBOL == "CACNG5" & !is.na(aso$SYMBOL),]
cacng5_aso <- cacng5_aso[order(cacng5_aso$pValue),][1,]
cat("ASO_CTRL_vs_Scramble_CTRL best CACNG5 row:\n")
cat("  proportion1 (cond_a=ASO_CTRL):", cacng5_aso$proportion1, "\n")
cat("  proportion2 (cond_b=Scramble_CTRL):", cacng5_aso$proportion2, "\n")
cat("  proportion2 - proportion1 =", cacng5_aso$proportion2 - cacng5_aso$proportion1, "\n")
cat("  regionType:", cacng5_aso$regionType, "\n")
