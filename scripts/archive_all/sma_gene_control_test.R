SMA_RELEVANT_GENES <- c(
  "SMN1", "SMN2", "GEMIN2", "GEMIN3", "GEMIN4", "GEMIN5", "GEMIN6", "GEMIN7", "GEMIN8",
  "UNRIP", "STRAP",
  "NCALD", "PLS3", "CHP1",
  "SIGMAR1", "VAPB", "DCTN1", "BSCL2", "HSPB1", "HSPB8", "GARS1", "REEP1",
  "TRPV4", "DNAJB2", "FBXO38", "ATP7A", "SLC5A7", "IGHMBP2", "AR",
  "CACNG1", "CACNG2", "CACNG3", "CACNG4", "CACNG5", "CACNG6", "CACNG7", "CACNG8",
  "PRKCA", "CHRNA1", "CHRNB1", "CHRND", "CHRNE", "CHRNG", "CHRNA9", "CHRNB3",
  "AGRN", "MUSK", "DOK7", "RAPSN", "COLQ", "LRP4", "GRIN2C", "GRIN3A",
  "SLC32A1", "GAD1", "GAD2", "GPHN",
  "GFRA1", "GFRA2", "GFRA3", "GFRA4", "RET", "GDNF",
  "ROCK1", "ROCK2", "DYNC1H1", "BICD2", "KIF5A", "TUBB3", "TUBB4A",
  "SRSF1", "HNRNPA1", "TIA1", "HTRA2", "ELAVL1"
)
N_SMA_GENES <- length(SMA_RELEVANT_GENES)
N_PERM <- 1000
set.seed(42)

CONTRAST <- "ASO_VPA_vs_Scramble_CTRL"
csv <- file.path("results/dmr_annotation", paste0(CONTRAST, "_annotated.csv"))
df <- read.csv(csv)
df <- df[!is.na(df$SYMBOL) & df$SYMBOL != "", ]

# Universe of all genes with a DMR call available to be tested (use the
# unique gene symbols seen anywhere in the annotated table, mirroring the
# space the SMA panel was drawn from -- not literally every gene in the
# genome, just genes that show up in this annotation universe at all)
all_genes_with_dmr <- unique(df$SYMBOL)
message("Total unique genes with at least one DMR call in this contrast: ",
        length(all_genes_with_dmr))

sma_hit <- sum(SMA_RELEVANT_GENES %in% all_genes_with_dmr)
sma_hit_rate <- sma_hit / N_SMA_GENES
message("SMA panel: ", sma_hit, " of ", N_SMA_GENES, " genes hit (",
        round(sma_hit_rate*100,1), "%)")

# Universe to sample random panels from: genes that appear ANYWHERE in the
# genome annotation (i.e. genes DMRcaller's annotation pipeline could in
# principle have called). Use the gene universe from ChIPseeker's TxDb
# annotation -- approximate with all unique SYMBOLs across ALL FIVE
# contrasts combined, which gives a broad, fair sampling universe.
all_contrasts <- c("ASO_VPA_vs_Scramble_CTRL","ASO_CTRL_vs_Scramble_CTRL",
                   "Scramble_VPA_vs_Scramble_CTRL","ASO_VPA_vs_ASO_CTRL",
                   "ASO_VPA_vs_Scramble_VPA")
full_gene_universe <- character(0)
for (cn in all_contrasts) {
  f <- file.path("results/dmr_annotation", paste0(cn, "_annotated.csv"))
  if (file.exists(f)) {
    d <- read.csv(f)
    full_gene_universe <- union(full_gene_universe, unique(d$SYMBOL[!is.na(d$SYMBOL)]))
  }
}
message("Full gene universe across all 5 contrasts: ", length(full_gene_universe), " genes")

# Permutation: draw N_SMA_GENES random genes from the full universe,
# check how many show up as hit in the target contrast, repeat N_PERM times
message("\nRunning ", N_PERM, " random panel draws (n=", N_SMA_GENES, " genes each)...")
perm_hits <- integer(N_PERM)
for (i in seq_len(N_PERM)) {
  random_panel <- sample(full_gene_universe, N_SMA_GENES)
  perm_hits[i] <- sum(random_panel %in% all_genes_with_dmr)
}

perm_pval <- mean(perm_hits >= sma_hit)

cat("\n=== Control test result ===\n")
cat(sprintf("SMA panel hit rate:    %d / %d (%.1f%%)\n", sma_hit, N_SMA_GENES, sma_hit_rate*100))
cat(sprintf("Random panel mean:     %.1f / %d (%.1f%%)\n",
            mean(perm_hits), N_SMA_GENES, mean(perm_hits)/N_SMA_GENES*100))
cat(sprintf("Random panel range:    [%d, %d]\n", min(perm_hits), max(perm_hits)))
cat(sprintf("Permutation p-value (SMA panel hit rate >= random): %.4f\n", perm_pval))

if (perm_pval < 0.05) {
  cat("\nRESULT: SMA panel is hit MORE OFTEN than random gene panels of the same size (p < 0.05)\n")
} else {
  cat("\nRESULT: SMA panel hit rate is consistent with what ANY random gene panel of this size would show (p >= 0.05)\n")
  cat("This means the high hit rate reflects VPA's broad genome-wide effect, not SMA-specific targeting\n")
}
