suppressPackageStartupMessages({
    library(DMRcaller)
    library(GenomicRanges)
})

setwd("/data/scratch/bt25018/sma_epigenomics_pipeline")

BY_CHR_DIR <- "results/alignments/bs/by_chr"
REPS       <- 1:3

LOCI <- list(
    SMN1 = list(chr = "chr5", start = 70924941, end = 70953015, strand = "+"),
    SMN2 = list(chr = "chr5", start = 70049638, end = 70078522, strand = "+")
)
FLANK <- 5000

EXONS <- list(
    SMN1 = data.frame(
        exon  = 1:9,
        start = c(70924941, 70938861, 70941290, 70942358, 70944614,
                  70946066, 70946948, 70950965, 70951942),
        end   = c(70925545, 70938977, 70941389, 70942437, 70944696,
                  70946186, 70947057, 70951088, 70953015)
    ),
    SMN2 = data.frame(
        exon  = 1:9,
        start = c(70049638, 70063488, 70065918, 70066986, 70069242,
                  70070697, 70071581, 70075600, 70076577),
        end   = c(70050242, 70063604, 70066017, 70067065, 70069324,
                  70070817, 70071690, 70075723, 70078522)
    )
)

WIN_SIZE <- 300

COMPARISONS <- list(
    list(name  = "ASO_vs_Scramble_CTRL",
         cond1 = "ASO_CTRL",  cond2 = "Scramble_CTRL",
         label = "ASO effect (CTRL background)"),
    list(name  = "VPA_vs_CTRL_ASO",
         cond1 = "ASO_CTRL",  cond2 = "ASO_VPA",
         label = "VPA effect (ASO background)"),
    list(name  = "ASO_vs_Scramble_VPA",
         cond1 = "ASO_VPA",   cond2 = "Scramble_VPA",
         label = "ASO effect (VPA background)")
)

NEEDED_CONDS <- unique(unlist(lapply(COMPARISONS, function(x) c(x$cond1, x$cond2))))

read_cpg_report <- function(f) {
    d <- read.table(gzfile(f), header = FALSE, sep = "\t",
                    col.names = c("chr","pos","strand","countM","countU","context","tri"),
                    colClasses = c("character","integer","character",
                                   "integer","integer","character","character"))
    GRanges(seqnames = d$chr, ranges = IRanges(d$pos, d$pos), strand = d$strand,
            readsM = d$countM, readsN = d$countM + d$countU,
            context = d$context, trinucleotide_context = d$tri)
}

read_condition <- function(condition) {
    files <- file.path(BY_CHR_DIR,
                       sprintf("%s_%d_chr5.CpG_report.txt.gz", condition, REPS))
    files <- files[file.exists(files)]
    if (length(files) != 3) warning(sprintf("%s: found %d/3 chr5 files", condition, length(files)))
    grs <- lapply(files, read_cpg_report)
    poolMethylationDatasets(GRangesList(grs))
}

message("Reading and pooling chr5 by condition ...")
pooled_data <- lapply(NEEDED_CONDS, function(c) { message("  ", c); read_condition(c) })
names(pooled_data) <- NEEDED_CONDS

build_gff <- function() {
    rows <- list()
    for (locus_name in names(LOCI)) {
        locus <- LOCI[[locus_name]]
        rows[[length(rows) + 1]] <- data.frame(
            chr = locus$chr, start = locus$start, end = locus$end,
            strand = locus$strand, type = "gene", name = locus_name,
            stringsAsFactors = FALSE)
        ex <- EXONS[[locus_name]]
        for (i in seq_len(nrow(ex))) {
            rows[[length(rows) + 1]] <- data.frame(
                chr = locus$chr, start = ex$start[i], end = ex$end[i],
                strand = locus$strand, type = "exon",
                name = sprintf("%s_exon%d", locus_name, ex$exon[i]),
                stringsAsFactors = FALSE)
        }
    }
    df <- do.call(rbind, rows)
    GRanges(seqnames = df$chr, ranges = IRanges(df$start, df$end),
            strand = df$strand, type = df$type, ID = df$name)
}

GEs <- build_gff()
message("Built GFF: ", length(GEs), " features")

plot_one <- function(comp, locus_name) {
    locus  <- LOCI[[locus_name]]
    region <- GRanges(seqnames = locus$chr,
                      ranges   = IRanges(locus$start - FLANK, locus$end + FLANK))
    title  <- sprintf("%s - %s vs %s (%s)",
                      locus_name, comp$cond1, comp$cond2, comp$label)
    plotLocalMethylationProfile(
        methylationData1 = pooled_data[[comp$cond1]],
        methylationData2 = pooled_data[[comp$cond2]],
        region           = region,
        DMRs             = NULL,
        conditionsNames  = c(comp$cond1, comp$cond2),
        gff              = GEs,
        windowSize       = WIN_SIZE,
        context          = "CG",
        main             = title,
        plotMeanLines    = TRUE,
        plotPoints       = TRUE
    )
    ex <- EXONS[[locus_name]]
    locus <- LOCI[[locus_name]]
    for (i in seq_len(nrow(ex))) {
        text((ex$start[i] + ex$end[i]) / 2, -0.13, labels = i, cex = 0.6, font = 2, col = ifelse(i == 7, "red", "black"))
    }
}

out_dir <- "results/qc/smn_locus"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (comp in COMPARISONS) {
    pdf_path <- file.path(out_dir, sprintf("smn_locus_dmrcaller_%s.pdf", comp$name))
    message("\nPlotting: ", comp$name, " -> ", pdf_path)
    pdf(pdf_path, width = 11, height = 8.5)
    par(mfrow = c(2, 1), mar = c(4, 4, 3, 1) + 0.1, cex = 0.9)
    plot_one(comp, "SMN1")
    plot_one(comp, "SMN2")
    dev.off()
    message("  Saved.")
}

combined_pdf <- file.path(out_dir, "smn_locus_dmrcaller_all_comparisons.pdf")
message("\nPlotting combined: ", combined_pdf)
pdf(combined_pdf, width = 16, height = 14)
par(mfrow = c(3, 2), mar = c(4, 4, 3, 1) + 0.1, cex = 0.7)
for (comp in COMPARISONS) {
    plot_one(comp, "SMN1")
    plot_one(comp, "SMN2")
}
dev.off()
message("  Saved.")

summarise <- function(condition_name, locus_name) {
    locus <- LOCI[[locus_name]]
    gr    <- pooled_data[[condition_name]]
    sel   <- as.character(seqnames(gr)) == locus$chr &
             start(gr) >= locus$start & start(gr) <= locus$end &
             gr$context == "CG" & gr$readsN > 0
    g <- gr[sel]
    if (length(g) == 0) return(data.frame(condition = condition_name, locus = locus_name,
                                           n_cpg = 0L, weighted_mean = NA_real_))
    data.frame(condition = condition_name, locus = locus_name,
               n_cpg = length(g),
               weighted_mean = round(sum(g$readsM) / sum(g$readsN), 4))
}

summary_rows <- list()
for (cond in NEEDED_CONDS) {
    for (locus_name in names(LOCI)) {
        summary_rows[[length(summary_rows) + 1]] <- summarise(cond, locus_name)
    }
}
summary_df <- do.call(rbind, summary_rows)
tsv_path <- file.path(out_dir, "smn_locus_dmrcaller_comparisons_summary.tsv")
write.table(summary_df, tsv_path, sep = "\t", quote = FALSE, row.names = FALSE)
message("\nMean methylation in gene body:")
print(summary_df, row.names = FALSE)
message("\nDone. Outputs in ", out_dir)
