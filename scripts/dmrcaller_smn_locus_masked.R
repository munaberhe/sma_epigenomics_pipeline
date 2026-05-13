#!/usr/bin/env Rscript
# =============================================================================
# dmrcaller_smn_locus_masked.R
# SMN2 locus methylation profile — SMN1-masked alignment
# Muna Berhe · bt25018 · QMUL MSc Bioinformatics
# Supervisor: Dr Radu Zabet
#
# Uses chr5 CpG reports from the SMN1-masked re-alignment to produce
# clean, unambiguous methylation profiles at the SMN2 locus.
#
# Produces:
#   - DMRcaller plotLocalMethylationProfile (3 comparisons, exon track)
#   - Exon 7 zoom plot (70,070,697–70,070,817)
#   - Updated methylation summary TSV
#
# Usage:
#   Rscript scripts/dmrcaller_smn_locus_masked.R
#   (run from project root on Apocrita after all 12 masked alignments done)
# =============================================================================

suppressPackageStartupMessages({
  library(DMRcaller)
  library(GenomicRanges)
})

.libPaths(c("~/R/library", .libPaths()))

# ── Configuration ─────────────────────────────────────────────────────────────
# Masked data — chr5 CpG reports from SMN1-masked re-alignment
BY_CHR_DIR_MASKED <- "results/alignments_smn1_masked/by_chr"
OUT_DIR           <- "results/qc/smn_locus"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

REPS <- 1:3

# hg38 SMN2 gene body coordinates (masked reads map here unambiguously)
SMN2 <- list(chr = "chr5", start = 70049638, end = 70078522)
FLANK <- 5000

# SMN2 exon coordinates (RefSeq hg38)
SMN2_EXONS <- data.frame(
  exon  = 1:9,
  start = c(70049638, 70063488, 70065918, 70066986, 70069242,
            70070513, 70070697, 70074672, 70076978),
  end   = c(70055558, 70063600, 70066017, 70067065, 70069324,
            70070628, 70070817, 70074808, 70078522)
)
# Exon 7: 70,070,697 – 70,070,817 (the therapeutic target window)
EXON7 <- list(start = 70070697, end = 70070817)

# GFF for plotLocalMethylationProfile
GFF_FILE <- "data/reference/hg38.ensGene.gtf"

# Three contrasts (same as genome-wide)
CONTRASTS <- list(
  list(name = "ASO_CTRL_vs_Scramble_CTRL",
       cond_a = "ASO_CTRL", cond_b = "Scramble_CTRL",
       label_a = "ASO_CTRL", label_b = "Scramble_CTRL",
       col_a = "#1B5478", col_b = "#6B7280"),
  list(name = "ASO_VPA_vs_Scramble_CTRL",
       cond_a = "ASO_VPA", cond_b = "Scramble_CTRL",
       label_a = "ASO_VPA", label_b = "Scramble_CTRL",
       col_a = "#C8820A", col_b = "#6B7280"),
  list(name = "Scramble_VPA_vs_Scramble_CTRL",
       cond_a = "Scramble_VPA", cond_b = "Scramble_CTRL",
       label_a = "Scramble_VPA", label_b = "Scramble_CTRL",
       col_a = "#6020A0", col_b = "#6B7280")
)

CONDITIONS_ALL <- c("ASO_CTRL", "ASO_VPA", "Scramble_CTRL", "Scramble_VPA")

# ── Step 1: Load masked chr5 CpG reports ──────────────────────────────────────
message("=== Loading masked chr5 CpG reports ===")
meth_raw <- list()

for (cond in CONDITIONS_ALL) {
  for (rep in REPS) {
    sample <- paste0(cond, "_", rep)
    path   <- file.path(BY_CHR_DIR_MASKED,
                        paste0(sample, "_chr5.CpG_report.txt.gz"))
    if (!file.exists(path)) {
      stop("Missing masked chr5 report: ", path,
           "\nHas the masked re-alignment finished for all 12 samples?")
    }
    message("  Loading: ", sample)
    meth_raw[[sample]] <- readBismark(path)
  }
}
message("All chr5 files loaded.")

# ── Step 2: Pool replicates per condition ─────────────────────────────────────
message("\n=== Pooling replicates ===")
meth_pooled <- list()
for (cond in CONDITIONS_ALL) {
  samples <- paste0(cond, "_", REPS)
  glist   <- GRangesList(lapply(samples, function(s) meth_raw[[s]]))
  meth_pooled[[cond]] <- poolMethylationDatasets(glist)
  n_cpg <- length(meth_pooled[[cond]])
  message("  ", cond, ": ", format(n_cpg, big.mark = ","), " CpGs on chr5")
}

# ── Step 3: Load GFF annotation ───────────────────────────────────────────────
message("\n=== Loading GFF annotation ===")
gff_all <- NULL
if (file.exists(GFF_FILE)) {
  library(rtracklayer)
  gff_all <- import(GFF_FILE, format = "GTF")
  gff_all <- gff_all[mcols(gff_all)$type %in% c("transcript", "exon")]
  mcols(gff_all)$ID <- mcols(gff_all)$gene_name
  message("  GFF loaded: ", length(gff_all), " features")
} else {
  message("  WARNING: GFF not found at ", GFF_FILE, " — plots will lack exon track")
}

# ── Step 4: DMRcaller local profiles — 3 comparisons ─────────────────────────
message("\n=== Plotting SMN2 local methylation profiles (masked) ===")

region_gr <- GRanges(
  seqnames = SMN2$chr,
  ranges   = IRanges(SMN2$start - FLANK, SMN2$end + FLANK)
)

plot_list <- list()

for (ct in CONTRASTS) {
  message("  Contrast: ", ct$name)
  pdf_path <- file.path(OUT_DIR,
    paste0("smn2_masked_", ct$name, ".pdf"))

  pdf(pdf_path, width = 12, height = 5)
  plotLocalMethylationProfile(
    methylationData1 = meth_pooled[[ct$cond_a]],
    methylationData2 = meth_pooled[[ct$cond_b]],
    region           = region_gr,
    context          = "CG",
    labels           = c(ct$label_a, ct$label_b),
    col              = c(ct$col_a, ct$col_b),
    GFF              = gff_all,
    windowSize       = 150
  )
  title(main = paste0("SMN2 locus (masked) — ", ct$label_a, " vs ", ct$label_b),
        sub  = "SMN1 masked at chr5:70,924,941-70,953,015 · Reads map unambiguously to SMN2")
  dev.off()
  message("    Saved: ", pdf_path)
  plot_list[[ct$name]] <- pdf_path
}

# Combined PDF — all 3 comparisons
message("\n  Generating combined PDF...")
combined_path <- file.path(OUT_DIR, "smn2_masked_all_comparisons.pdf")
pdf(combined_path, width = 12, height = 5)
for (ct in CONTRASTS) {
  plotLocalMethylationProfile(
    methylationData1 = meth_pooled[[ct$cond_a]],
    methylationData2 = meth_pooled[[ct$cond_b]],
    region           = region_gr,
    context          = "CG",
    labels           = c(ct$label_a, ct$label_b),
    col              = c(ct$col_a, ct$col_b),
    GFF              = gff_all,
    windowSize       = 150
  )
  title(main = paste0("SMN2 (masked) — ", ct$label_a, " vs ", ct$label_b))
}
dev.off()
message("  Combined PDF saved: ", combined_path)

# ── Step 5: Exon 7 zoom ───────────────────────────────────────────────────────
message("\n=== Exon 7 zoom (70,070,697–70,070,817) ===")

exon7_gr <- GRanges(
  seqnames = "chr5",
  ranges   = IRanges(EXON7$start - 500, EXON7$end + 500)
)

exon7_pdf <- file.path(OUT_DIR, "smn2_exon7_zoom_masked.pdf")
pdf(exon7_pdf, width = 10, height = 5)
plotLocalMethylationProfile(
  methylationData1 = meth_pooled[["ASO_VPA"]],
  methylationData2 = meth_pooled[["Scramble_CTRL"]],
  region           = exon7_gr,
  context          = "CG",
  labels           = c("ASO_VPA", "Scramble_CTRL"),
  col              = c("#C8820A", "#6B7280"),
  GFF              = gff_all,
  windowSize       = 50
)
# Shade exon 7 window
rect(EXON7$start, par("usr")[3], EXON7$end, par("usr")[4],
     col = adjustcolor("red", 0.08), border = NA)
abline(v = c(EXON7$start, EXON7$end), col = "red", lty = 2, lwd = 1)
title(main = "SMN2 exon 7 — CpG methylation (masked)",
      sub  = paste0("Exon 7 window: ", EXON7$start, "–", EXON7$end,
                    " · nusinersen target · red shading = exon 7"))
dev.off()
message("Exon 7 zoom saved: ", exon7_pdf)

# ── Step 6: Methylation summary table ────────────────────────────────────────
message("\n=== Methylation summary (gene body, masked) ===")

compute_weighted_mean <- function(gr, chr, start, end) {
  region <- gr[seqnames(gr) == chr &
               start(gr) >= start &
               end(gr)   <= end]
  if (length(region) == 0) return(NA)
  total_reads <- sum(region$readsN, na.rm = TRUE)
  if (total_reads == 0) return(NA)
  sum(region$readsM, na.rm = TRUE) / total_reads
}

summary_rows <- list()
for (cond in CONDITIONS_ALL) {
  meth  <- meth_pooled[[cond]]
  n_cpg <- sum(seqnames(meth) == "chr5" &
               start(meth) >= SMN2$start &
               end(meth)   <= SMN2$end,
               na.rm = TRUE)
  wmean <- compute_weighted_mean(meth, "chr5", SMN2$start, SMN2$end)
  exon7_mean <- compute_weighted_mean(meth, "chr5", EXON7$start, EXON7$end)

  summary_rows[[cond]] <- data.frame(
    condition       = cond,
    SMN2_weighted_mean = round(wmean, 4),
    SMN2_n_CpGs     = n_cpg,
    exon7_mean      = round(exon7_mean, 4)
  )
  message(sprintf("  %-20s SMN2=%.3f  exon7=%.3f  n=%d",
                  cond, wmean, exon7_mean, n_cpg))
}

summary_df  <- do.call(rbind, summary_rows)
summary_tsv <- file.path(OUT_DIR, "smn2_masked_methylation_summary.tsv")
write.table(summary_df, summary_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
message("\nSummary saved: ", summary_tsv)
print(summary_df)

message("\n=== SMN2 masked analysis complete ===")
message("Outputs in: ", OUT_DIR)
message("Key files:")
message("  ", combined_path,     "  ← show to Radu")
message("  ", exon7_pdf,         "  ← exon 7 therapeutic target")
message("  ", summary_tsv,       "  ← update thesis table")
