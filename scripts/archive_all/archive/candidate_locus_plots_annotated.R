.libPaths("~/R/library")
# candidate_locus_plots_annotated.R
# Locus-level methylation profiles for candidate DMR genes
# Monkey-patches .plotGeneticElements for staggered gene labels
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

library(DMRcaller)
library(rtracklayer)

# ── Patch .plotGeneticElements to stagger overlapping gene labels ─────────────
# Original places all labels at start(gene) at a fixed y — causes collision
# when genes are close. This version alternates y between two levels.
.plotGeneticElements_patched <- function(gff, region, col) {
  seqname <- seqnames(region)
  minPos  <- start(region)
  maxPos  <- end(region)
  gff     <- gff[queryHits(findOverlaps(gff, region))]
  start(gff) <- pmax(start(gff), minPos)
  end(gff)   <- pmin(end(gff),   maxPos)

  genes    <- gff[gff$type == "gene"]
  genesPos <- genes[strand(genes) == "+" | strand(genes) == "*"]
  genesNeg <- genes[strand(genes) == "-" | strand(genes) == "*"]
  exons    <- gff[gff$type == "exon"]
  exons    <- exons[overlapsAny(exons, genes)]
  exonsPos <- exons[strand(exons) == "+" | strand(exons) == "*"]
  exonsNeg <- exons[strand(exons) == "-" | strand(exons) == "*"]
  transposons    <- gff[gff$type == "transposable_element"]
  transposonsPos <- transposons[strand(transposons) == "+" | strand(transposons) == "*"]
  transposonsNeg <- transposons[strand(transposons) == "-" | strand(transposons) == "*"]

  negativeStrandPosition <- -0.175
  positiveStrandPosition <- -0.075

  text(maxPos + (maxPos - minPos)/100, positiveStrandPosition, "+")
  text(maxPos + (maxPos - minPos)/100, negativeStrandPosition, "-")
  lines(c(minPos, maxPos), c(-0.14, -0.14), lty = 1, lwd = 0.75, col = "black")

  # stagger_labels: alternate y position for successive genes to avoid collision
  stagger_labels <- function(genes_gr, y_levels) {
    if (length(genes_gr) == 0) return(invisible(NULL))
    xs <- start(genes_gr)
    ids <- genes_gr$ID
    # sort by position so staggering is applied left-to-right
    ord <- order(xs)
    xs  <- xs[ord]
    ids <- ids[ord]
    for (i in seq_along(xs)) {
      y <- y_levels[((i - 1) %% length(y_levels)) + 1]
      text(xs[i], y, ids[i], pos = 4, cex = 0.5)
    }
  }

  if (length(genesPos) > 0) {
    segments(start(genesPos), positiveStrandPosition, end(genesPos), positiveStrandPosition)
    # stagger between two y levels above the + strand line
    stagger_labels(genesPos, c(-0.10, -0.125))
  }
  if (length(genesNeg) > 0) {
    segments(start(genesNeg), -0.175, end(genesNeg), negativeStrandPosition)
    # stagger between two y levels below the - strand line
    stagger_labels(genesNeg, c(-0.22, -0.245))
  }
  if (length(exonsPos) > 0)
    rect(start(exonsPos), -0.05, end(exonsPos), -0.09, col = col[1], border = NA)
  if (length(exonsNeg) > 0)
    rect(start(exonsNeg), -0.16, end(exonsNeg), -0.2,  col = col[1], border = NA)
  if (length(transposonsPos) > 0) {
    rect(start(transposonsPos), -0.05, end(transposonsPos), -0.09,
         col = col[2], border = col[2], density = 30, angle = 30)
    stagger_labels(transposonsPos, c(-0.10, -0.125))
  }
  if (length(transposonsNeg) > 0) {
    rect(start(transposonsNeg), -0.16, end(transposonsNeg), -0.2,
         col = col[2], border = col[2], density = 30, angle = 30)
    stagger_labels(transposonsNeg, c(-0.22, -0.245))
  }
}

# Inject patched function into DMRcaller's private namespace
environment(.plotGeneticElements_patched) <- asNamespace("DMRcaller")
assignInNamespace(".plotGeneticElements", .plotGeneticElements_patched, ns = "DMRcaller")
message("Patched .plotGeneticElements with staggered label version")
# ─────────────────────────────────────────────────────────────────────────────

COV_DIR <- "results/alignments/bs/by_chr"
OUT_DIR <- "results/dmr/candidate_locus_plots_annotated"
GTF <- "/data/home/bt25018/sma_epigenomics_pipeline/data/reference/Homo_sapiens.GRCh38.110.gtf.gz"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SAMPLES <- data.frame(
  name = c("Scramble_CTRL_1", "Scramble_CTRL_2", "Scramble_CTRL_3",
           "Scramble_VPA_1", "Scramble_VPA_2", "Scramble_VPA_3",
           "ASO_CTRL_1", "ASO_CTRL_2", "ASO_CTRL_3",
           "ASO_VPA_1", "ASO_VPA_2", "ASO_VPA_3"),
  group = c(rep("Scramble_CTRL", 3), rep("Scramble_VPA", 3),
            rep("ASO_CTRL", 3), rep("ASO_VPA", 3)),
  stringsAsFactors = FALSE
)

load_group <- function(group_name, chrom) {
  samples <- SAMPLES$name[SAMPLES$group == group_name]
  paths <- file.path(COV_DIR, paste0(samples, "_", chrom, ".CpG_report.txt.gz"))
  paths <- paths[file.exists(paths)]
  if (length(paths) < 2) return(NULL)
  tryCatch({
    dat <- readBismarkPool(paths)
    dat[dat$readsN >= 4]
  }, error = function(e) NULL)
}

message("Loading GTF annotation...")
ge_all <- import(GTF, format = "GTF")
seqlevels(ge_all) <- paste0("chr", seqlevels(ge_all))
ge_all <- ge_all[ge_all$type %in% c("gene", "exon")]
ge_all$ID <- ifelse(!is.na(ge_all$gene_name) & ge_all$gene_name != "",
                    ge_all$gene_name, ge_all$gene_id)
message("  Features loaded: ", length(ge_all))

make_locus_plot <- function(gene, chrom, dmr_start, dmr_end,
                            treatment_group, control_group,
                            contrast_name, window = 15000) {
  message("Plotting: ", gene, " (", chrom, ":", dmr_start, "-", dmr_end, ")")

  t_dat <- load_group(treatment_group, chrom)
  c_dat <- load_group(control_group, chrom)
  if (is.null(t_dat) || is.null(c_dat)) {
    message("  Skipping — data load failed")
    return(NULL)
  }

  region <- GRanges(seqnames = Rle(chrom),
                    ranges = IRanges(dmr_start - window, dmr_end + window))
  ge_local <- ge_all[seqnames(ge_all) == chrom]

  ge_window <- subsetByOverlaps(ge_local[ge_local$type == "gene"], region)
  message("  Genes in window: ", paste(ge_window$ID, collapse=", "))

  rds_path <- file.path("results/dmr", contrast_name,
                        paste0(contrast_name, "_all_chr.rds"))
  dmr_list <- if (file.exists(rds_path)) {
    GRangesList(DMRs = readRDS(rds_path))
  } else {
    GRangesList()
  }

  outfile <- file.path(OUT_DIR, paste0(gene, "_", contrast_name, "_annotated_locus.pdf"))

  tryCatch({
    cairo_pdf(outfile, width = 16, height = 9)
    par(mar = c(6, 4, 3, 1) + 0.1)
    plotLocalMethylationProfile(
      t_dat, c_dat,
      region,
      dmr_list,
      conditionsNames = c(treatment_group, control_group),
      gff = ge_local,
      windowSize = 500,
      main = paste0(gene, " ... CG methylation ... ", treatment_group, " vs ", control_group))
    dev.off()
    message("  Saved: ", outfile)
  }, error = function(e) {
    message("  Plot error: ", e$message)
    try(dev.off(), silent = TRUE)
  })
}

# ASO_VPA vs ASO_CTRL
make_locus_plot("EZH1",   "chr17", 42742132,  42742344,  "ASO_VPA", "ASO_CTRL", "ASO_VPA_vs_ASO_CTRL")
make_locus_plot("TBP",    "chr6",  170558739, 170558882, "ASO_VPA", "ASO_CTRL", "ASO_VPA_vs_ASO_CTRL")
make_locus_plot("TUBB2B", "chr6",  3224451,   3224556,   "ASO_VPA", "ASO_CTRL", "ASO_VPA_vs_ASO_CTRL")
make_locus_plot("ABCB1",  "chr7",  87599820,  87599893,  "ASO_VPA", "ASO_CTRL", "ASO_VPA_vs_ASO_CTRL")
make_locus_plot("CAST",   "chr5",  95963581,  95963699,  "ASO_VPA", "ASO_CTRL", "ASO_VPA_vs_ASO_CTRL")
make_locus_plot("DDAH1",  "chr1",  85466937,  85467044,  "ASO_VPA", "ASO_CTRL", "ASO_VPA_vs_ASO_CTRL")

# VPA effect — Scramble_VPA vs Scramble_CTRL
make_locus_plot("FYN",    "chr6",  111778547, 111778734, "Scramble_VPA", "Scramble_CTRL", "VPA_vs_Scramble_CTRL")
make_locus_plot("SRPK1",  "chr6",  35893849,  35893980,  "Scramble_VPA", "Scramble_CTRL", "VPA_vs_Scramble_CTRL")
make_locus_plot("CACNB2", "chr10", 18397445,  18397587,  "Scramble_VPA", "Scramble_CTRL", "VPA_vs_Scramble_CTRL")
make_locus_plot("SATB1",  "chr3",  18428899,  18428989,  "Scramble_VPA", "Scramble_CTRL", "VPA_vs_Scramble_CTRL")
make_locus_plot("GTF2B",  "chr1",  88886665,  88886846,  "Scramble_VPA", "Scramble_CTRL", "VPA_vs_Scramble_CTRL")

# ASO + VPA combined — ASO_VPA vs Scramble_CTRL
make_locus_plot("FOXP1",   "chr3",  71689445,  71689619,  "ASO_VPA", "Scramble_CTRL", "ASO_VPA_vs_Scramble_CTRL")
make_locus_plot("PRPF38B", "chr1",  108700037, 108700113, "ASO_VPA", "Scramble_CTRL", "ASO_VPA_vs_Scramble_CTRL")
make_locus_plot("DSC2",    "chr18", 31099409,  31099652,  "ASO_VPA", "Scramble_CTRL", "ASO_VPA_vs_Scramble_CTRL")
make_locus_plot("ACBD6",   "chr1",  180413842, 180413989, "ASO_VPA", "Scramble_CTRL", "ASO_VPA_vs_Scramble_CTRL")

message("\nAll annotated locus plots saved to: ", OUT_DIR)
