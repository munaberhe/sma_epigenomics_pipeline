
.libPaths("~/R/library")
# filter_dmrs.R
# DMR analysis using filterDMRs on predefined promoter regions
# Complements computeDMRs — tests promoters directly as units
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

library(DMRcaller)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(GenomeInfoDb)

COV_DIR <- "results/alignments/bs/by_chr"
OUT_DIR <- "results/dmr_filterDMRs"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Parameters
MIN_COV                <- 4
CONTEXT                <- "CG"
P_VALUE_THRESHOLD      <- 0.01
MIN_CYTOSINES_COUNT    <- 4
MIN_PROPORTION_DIFF    <- 0.1
MIN_READS_PER_CYTOSINE <- 4
CHROMOSOMES            <- paste0("chr", c(1:22, "X"))

SAMPLES <- data.frame(
  name  = c("Scramble_CTRL_1", "Scramble_CTRL_2", "Scramble_CTRL_3",
            "Scramble_VPA_1",  "Scramble_VPA_2",  "Scramble_VPA_3",
            "ASO_CTRL_1",      "ASO_CTRL_2",      "ASO_CTRL_3",
            "ASO_VPA_1",       "ASO_VPA_2",       "ASO_VPA_3"),
  group = c(rep("Scramble_CTRL", 3), rep("Scramble_VPA", 3),
            rep("ASO_CTRL", 3),      rep("ASO_VPA", 3)),
  stringsAsFactors = FALSE
)

CONTRASTS <- list(
  ASO_VPA_vs_ASO_CTRL      = list(treatment = "ASO_VPA",      control = "ASO_CTRL"),
  VPA_vs_Scramble_CTRL     = list(treatment = "Scramble_VPA", control = "Scramble_CTRL"),
  ASO_VPA_vs_Scramble_CTRL = list(treatment = "ASO_VPA",      control = "Scramble_CTRL")
)

# Build promoter regions from TxDb — 2kb upstream, 200bp downstream
message("Building promoter regions from TxDb hg38...")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
promoter_regions <- promoters(txdb, upstream = 2000, downstream = 200)

# Keep only standard chromosomes
promoter_regions <- keepSeqlevels(promoter_regions,
                                  CHROMOSOMES,
                                  pruning.mode = "coarse")

# Add gene symbols via TxDb gene mapping
tx2gene <- suppressMessages(select(txdb,
  keys    = keys(txdb, "TXID"),
  columns = c("TXID", "GENEID"),
  keytype = "TXID"))
entrez_ids <- tx2gene$GENEID[match(mcols(promoter_regions)$tx_id,
                                    tx2gene$TXID)]
symbols <- suppressMessages(
  mapIds(org.Hs.eg.db, keys = as.character(entrez_ids),
         column = "SYMBOL", keytype = "ENTREZID", multiVals = "first"))
mcols(promoter_regions)$SYMBOL   <- symbols
mcols(promoter_regions)$ENTREZID <- entrez_ids

message("  Promoter regions: ", length(promoter_regions))
message("  Chromosomes: ", paste(CHROMOSOMES, collapse=", "))

# Load and pool replicates for one group on one chromosome
load_group_chromosome <- function(sample_names, chrom) {
  paths <- file.path(COV_DIR,
    paste0(sample_names, "_", chrom, ".CpG_report.txt.gz"))
  paths <- paths[file.exists(paths)]
  if (length(paths) < 2) return(NULL)
  tryCatch({
    dat <- readBismarkPool(paths)
    dat[dat$readsN >= MIN_COV]
  }, error = function(e) { message("    Load error: ", e$message); NULL })
}

for (contrast_name in names(CONTRASTS)) {
  contrast <- CONTRASTS[[contrast_name]]
  message("\n=== filterDMRs: ", contrast_name,
          " (", contrast$treatment, " vs ", contrast$control, ") ===")

  contrast_out <- file.path(OUT_DIR, contrast_name)
  dir.create(contrast_out, recursive = TRUE, showWarnings = FALSE)

  all_filtered <- list()

  for (chrom in CHROMOSOMES) {
    message("  Chromosome: ", chrom)

    treatment_samples <- SAMPLES$name[SAMPLES$group == contrast$treatment]
    control_samples   <- SAMPLES$name[SAMPLES$group == contrast$control]

    t_dat <- load_group_chromosome(treatment_samples, chrom)
    c_dat <- load_group_chromosome(control_samples,   chrom)

    if (is.null(t_dat) || is.null(c_dat)) {
      message("    Skipping — data load failed")
      next
    }

    # Subset promoters to this chromosome
    chrom_promoters <- promoter_regions[seqnames(promoter_regions) == chrom]
    if (length(chrom_promoters) == 0) next
    message("    Promoters on ", chrom, ": ", length(chrom_promoters))

    filtered <- tryCatch({
      filterDMRs(
        t_dat,
        c_dat,
        potentialDMRs        = chrom_promoters,
        context              = CONTEXT,
        test                 = "score",
        pValueThreshold      = P_VALUE_THRESHOLD,
        minCytosinesCount    = MIN_CYTOSINES_COUNT,
        minProportionDifference = MIN_PROPORTION_DIFF,
        minReadsPerCytosine  = MIN_READS_PER_CYTOSINE,
        cores                = 4
      )
    }, error = function(e) {
      message("    filterDMRs error: ", e$message)
      NULL
    })

    if (!is.null(filtered) && length(filtered) > 0) {
      message("    Differentially methylated promoters: ", length(filtered))
      all_filtered[[chrom]] <- filtered
    }
  }

  if (length(all_filtered) == 0) {
    message("  No differentially methylated promoters for ", contrast_name)
    next
  }

  # Combine across chromosomes
  message("  Combining ", length(all_filtered), " chromosomes...")
  combined <- do.call(c, unname(all_filtered))

  saveRDS(combined, file.path(contrast_out,
                              paste0(contrast_name, "_filterDMRs.rds")))

  # Save as CSV with gene symbols
  df <- as.data.frame(combined)
  # Flatten any list columns to character
  for (col in names(df)) {
    if (is.list(df[[col]])) {
      df[[col]] <- sapply(df[[col]], function(x)
        if (length(x) == 0) NA else paste(x, collapse=";"))
    }
  }
  df$SYMBOL <- as.character(mcols(combined)$SYMBOL)
  df$proportion_diff <- df$proportion2 - df$proportion1
  df <- df[order(abs(df$proportion_diff), decreasing = TRUE), ]
  write.csv(df, file.path(contrast_out,
                          paste0(contrast_name, "_filterDMRs.csv")),
            row.names = FALSE)

  message("  Total differentially methylated promoters: ", length(combined))
  message("  Hypermethylated (gain): ", sum(df$regionType == "gain", na.rm=TRUE))
  message("  Hypomethylated (loss): ", sum(df$regionType == "loss", na.rm=TRUE))

  # Top 20 by methylation difference
  message("\n  Top 20 by methylation difference:")
  top20 <- head(df[!is.na(df$SYMBOL), ], 20)
  for (i in seq_len(nrow(top20))) {
    message(sprintf("    %2d. %-12s  treat=%.1f%%  ctrl=%.1f%%  diff=%.1f%%",
                    i, top20$SYMBOL[i],
                    top20$proportion1[i]*100,
                    top20$proportion2[i]*100,
                    abs(top20$proportion_diff[i])*100))
  }
}

message("\nDone. Outputs in: ", OUT_DIR)
