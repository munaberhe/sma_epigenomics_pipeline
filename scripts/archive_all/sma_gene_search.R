# Curated list of genes with established SMA / motor neuron / neuromuscular
# junction relevance, compiled from literature checked tonight plus core
# SMA pathway genes from the broader literature.
SMA_RELEVANT_GENES <- c(
  # Core SMA pathway / SMN-related
  "SMN1", "SMN2", "GEMIN2", "GEMIN3", "GEMIN4", "GEMIN5", "GEMIN6", "GEMIN7", "GEMIN8",
  "UNRIP", "STRAP",
  # Published SMA modifiers / disease-modifying genes
  "NCALD", "PLS3", "CHP1",
  # Motor neuron disease genes (ALS/SMA overlap, distal HMN)
  "SIGMAR1", "VAPB", "DCTN1", "BSCL2", "HSPB1", "HSPB8", "GARS1", "REEP1",
  "TRPV4", "DNAJB2", "FBXO38", "ATP7A", "SLC5A7", "IGHMBP2", "AR",
  # Neuromuscular junction / AMPA-glutamate signalling (your CACNG cluster family)
  "CACNG1", "CACNG2", "CACNG3", "CACNG4", "CACNG5", "CACNG6", "CACNG7", "CACNG8",
  "PRKCA", "CHRNA1", "CHRNB1", "CHRND", "CHRNE", "CHRNG", "CHRNA9", "CHRNB3",
  "AGRN", "MUSK", "DOK7", "RAPSN", "COLQ", "LRP4", "GRIN2C", "GRIN3A",
  # GABA / inhibitory motor circuit
  "SLC32A1", "GAD1", "GAD2", "GPHN",
  # GDNF / motor neuron survival signalling
  "GFRA1", "GFRA2", "GFRA3", "GFRA4", "RET", "GDNF",
  # Cytoskeleton / axon transport relevant to motor neuron disease
  "ROCK1", "ROCK2", "DYNC1H1", "BICD2", "KIF5A", "TUBB3", "TUBB4A",
  # Splicing machinery relevant to SMN2 exon 7 regulation
  "SRSF1", "HNRNPA1", "TIA1", "HTRA2", "ELAVL1"
)

for (cn in c("ASO_VPA_vs_Scramble_CTRL", "ASO_CTRL_vs_Scramble_CTRL",
             "Scramble_VPA_vs_Scramble_CTRL", "ASO_VPA_vs_ASO_CTRL",
             "ASO_VPA_vs_Scramble_VPA")) {
  csv <- file.path("results/dmr_annotation", paste0(cn, "_annotated.csv"))
  if (!file.exists(csv)) next
  df <- read.csv(csv)
  df <- df[!is.na(df$SYMBOL) & df$SYMBOL %in% SMA_RELEVANT_GENES, ]
  if (nrow(df) == 0) next
  df <- df[order(df$pValue), ]
  df <- df[!duplicated(df$SYMBOL), ]
  cat("\n==========  ", cn, "  ==========\n")
  print(df[, c("SYMBOL","seqnames","start","pValue","regionType","annotation","distanceToTSS")],
        row.names=FALSE)
}
