setwd("/data/home/bt25018/sma_epigenomics_pipeline")
DMR_FILES <- list(
  locked    = "results/dmr_annotation/ASO_VPA_vs_Scramble_CTRL_annotated.csv",
  sensitive = "results/dmr_annotation/ASO_VPA_vs_Scramble_CTRL_annotated.csv"
)
LOCUS_CHR <- "chr5"; QUERY_START <- 70080000; QUERY_END <- 70110000
OUT <- "results/smn2_enhancer/smn2_dmr_diagnostic.txt"
dir.create(dirname(OUT), recursive=TRUE, showWarnings=FALSE)
sink(OUT, split=TRUE)
cat("SMN2 3' DMR diagnostic\n======================\n")
cat(sprintf("Search window: %s:%d-%d\n\n", LOCUS_CHR, QUERY_START, QUERY_END))
for (label in names(DMR_FILES)) {
  f <- DMR_FILES[[label]]
  cat(sprintf("--- %s threshold (%s) ---\n", label, f))
  if (!file.exists(f)) { cat(sprintf("FILE NOT FOUND: %s\n\n", f)); next }
  d <- read.csv(f, stringsAsFactors=FALSE)
  keep <- d$seqnames==LOCUS_CHR & d$end>=QUERY_START & d$start<=QUERY_END
  sub <- d[keep,,drop=FALSE]
  cat(sprintf("DMRs in window: %d\n", nrow(sub)))
  if (nrow(sub)>0) {
    cols <- intersect(c("seqnames","start","end","regionType","proportion1",
                        "proportion2","cytosinesCount","pValue","SYMBOL","annotation"),
                      colnames(sub))
    print(sub[,cols])
    for (i in seq_len(nrow(sub))) {
      r <- sub[i,]
      delta <- (r$proportion2-r$proportion1)*100
      cat(sprintf("  row %d: %s:%d-%d delta=%+.1f%% gene=%s\n",
                  i,r$seqnames,r$start,r$end,delta,r$SYMBOL))
    }
  }
  cat("\n")
}
sink()
