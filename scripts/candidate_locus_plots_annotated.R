.libPaths("~/R/library")
# candidate_locus_plots_annotated.R
# Locus-level methylation profiles for candidate DMR genes
# Uses Ensembl GRCh38 GTF for gene/exon annotation on each plot
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

library(DMRcaller)
library(rtracklayer)

COV_DIR <- "results/alignments/bs/by_chr"
OUT_DIR <- "results/dmr/candidate_locus_plots_annotated"
GTF <- "/data/home/bt25018/sma_epigenomics_pipeline/data/reference/Homo_sapiens.GRCh38.109.chr.gtf.gz"
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

# load and pool replicates for a group, filter to min 4x coverage
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

# load GTF once — takes a while but only need to do it once
# originally used hg38.ensGene.gtf but gene_name was just Ensembl IDs
# switched to Ensembl release 109 which has proper gene symbols
message("Loading GTF annotation...")
ge_all <- import(GTF, format = "GTF")
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

  # tried 10000 first but 15000 shows more exon structure for larger genes like EZH1
  region <- GRanges(seqnames = Rle(chrom),
                    ranges = IRanges(dmr_start - window, dmr_end + window))

  # subset to chromosome so DMRcaller doesn't have to search the whole GTF
  ge_local <- ge_all[seqnames(ge_all) == chrom]
  message("  Gene features on ", chrom, ": ", length(ge_local))

  # load the DMR calls for this contrast if they exist
  rds_path <- file.path("results/dmr", contrast_name,
                        paste0(contrast_name, "_all_chr.rds"))
  dmr_list <- if (file.exists(rds_path)) {
    GRangesList(DMRs = readRDS(rds_path))
  } else {
    GRangesList()
  }

  outfile <- file.path(OUT_DIR, paste0(gene, "_", contrast_name, "_annotated_locus.pdf"))

  tryCatch({
    pdf(outfile, width = 14, height = 7)
    par(mar = c(4, 4, 3, 1) + 0.1)
    plotLocalMethylationProfile(
      t_dat, c_dat,
      region,
      dmr_list,
      conditionsNames = c(treatment_group, control_group),
      gff = ge_local,
      windowSize = 500,  # Radu suggested up to 500bp based on correlation plots
      main = paste0(gene, " — CG methylation — ", treatment_group, " vs ", control_group))
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

# ASO + VPA combined effect — ASO_VPA vs Scramble_CTRL
make_locus_plot("FOXP1",   "chr3",  71689445,  71689619,  "ASO_VPA", "Scramble_CTRL", "ASO_VPA_vs_Scramble_CTRL")
make_locus_plot("PRPF38B", "chr1",  108700037, 108700113, "ASO_VPA", "Scramble_CTRL", "ASO_VPA_vs_Scramble_CTRL")
make_locus_plot("DSC2",    "chr18", 31099409,  31099652,  "ASO_VPA", "Scramble_CTRL", "ASO_VPA_vs_Scramble_CTRL")
make_locus_plot("ACBD6",   "chr1",  180413842, 180413989, "ASO_VPA", "Scramble_CTRL", "ASO_VPA_vs_Scramble_CTRL")

message("\nAll annotated locus plots saved to: ", OUT_DIR)
