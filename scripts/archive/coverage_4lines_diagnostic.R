#!/usr/bin/env Rscript
.libPaths("~/R/library")

suppressPackageStartupMessages({
    library(DMRcaller)
    library(GenomicRanges)
    library(ggplot2)
})

setwd("/data/scratch/bt25018/sma_epigenomics_pipeline")

BY_CHR_DIR <- "results/alignments/bs/by_chr"
CONDITION  <- "ASO_CTRL"
REPS       <- c(1, 2, 3)
CHROMS     <- paste0("chr", c(1:22, "X", "Y"))

read_cpg_report <- function(f) {
    d <- read.table(gzfile(f), header = FALSE, sep = "\t",
                    col.names = c("chr","pos","strand","countM","countU","context","tri"),
                    colClasses = c("character","integer","character",
                                   "integer","integer","character","character"))
    GRanges(seqnames = d$chr,
            ranges   = IRanges(d$pos, d$pos),
            strand   = d$strand,
            readsM   = d$countM,
            readsN   = d$countM + d$countU,
            context  = d$context,
            trinucleotide_context = d$tri)
}

read_replicate_safe <- function(condition, rep, chroms, dir) {
    files <- file.path(dir, sprintf("%s_%d_%s.CpG_report.txt.gz",
                                    condition, rep, chroms))
    files <- files[file.exists(files)]
    if (length(files) == 0) stop("No files found for ", condition, "_", rep)
    message(sprintf("  %s rep %d: %d per-chr files", condition, rep, length(files)))
    grs <- lapply(files, read_cpg_report)
    do.call(c, grs)
}

message("Reading replicates for ", CONDITION, " ...")
rep_data <- lapply(REPS, function(r) {
    read_replicate_safe(CONDITION, r, CHROMS, BY_CHR_DIR)
})
names(rep_data) <- paste0("rep", REPS)

n_cpg <- sapply(rep_data, length)
message("CpGs per replicate: ", paste(n_cpg, collapse = ", "))
stopifnot(length(unique(n_cpg)) == 1)

message("Pooling rep1 + rep2 + rep3 ...")
pooled <- poolMethylationDatasets(GRangesList(rep_data))

make_df <- function(gr, label) {
    data.frame(
        coverage = gr$readsN,
        track    = label,
        stringsAsFactors = FALSE
    )
}

df <- rbind(
    make_df(rep_data[[1]], "rep1"),
    make_df(rep_data[[2]], "rep2"),
    make_df(rep_data[[3]], "rep3"),
    make_df(pooled,        "pooled")
)
df$track <- factor(df$track, levels = c("rep1", "rep2", "rep3", "pooled"))
df <- df[df$coverage > 0, ]

# Sanity check: total reads
total_reads <- tapply(df$coverage, df$track, sum)
message("\nTotal reads per track:")
print(total_reads)
message("rep1+rep2+rep3 sum: ", sum(total_reads[c("rep1","rep2","rep3")]))
message("pooled total:       ", total_reads["pooled"])

XMAX <- 100

# Flip y-axis: fraction of CpGs >= threshold (matches existing QC convention)
p <- ggplot(df, aes(x = coverage, colour = track)) +
    stat_ecdf(geom = "line", linewidth = 0.9, pad = FALSE,
              aes(y = 1 - after_stat(y))) +
    coord_cartesian(xlim = c(0, XMAX)) +
    scale_colour_manual(
        values = c(rep1   = "#7CB6D6",
                   rep2   = "#4A8FB8",
                   rep3   = "#1B5478",
                   pooled = "#A84B2F")
    ) +
    geom_vline(xintercept = 10, linetype = 2, colour = "grey40") +
    annotate("text", x = 11, y = 0.55, label = "10x", hjust = 0,
             colour = "grey40", size = 3.5) +
    labs(
        title = sprintf("CpG coverage: replicates vs pooled (%s)", CONDITION),
        x     = "Coverage depth",
        y     = "Proportion of CpGs",
        colour = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
        legend.position    = "right",
        panel.grid.major.y = element_line(colour = "grey92"),
        plot.title         = element_text(face = "bold")
    )

out_dir <- "results/qc/coverage_4lines"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

png_path <- file.path(out_dir, sprintf("coverage_4lines_%s.png", CONDITION))
pdf_path <- file.path(out_dir, sprintf("coverage_4lines_%s.pdf", CONDITION))

ggsave(png_path, p, width = 7, height = 4.5, dpi = 200)
ggsave(pdf_path, p, width = 7, height = 4.5)

message("Saved: ", png_path)
message("Saved: ", pdf_path)

summarise_track <- function(label) {
    cov <- df$coverage[df$track == label]
    data.frame(
        track      = label,
        n_cpg      = length(cov),
        mean       = round(mean(cov), 2),
        median     = median(cov),
        q90        = quantile(cov, 0.9, names = FALSE),
        pct_ge_10x = round(100 * mean(cov >= 10), 1)
    )
}
summary_df <- do.call(rbind, lapply(levels(df$track), summarise_track))
message("\nNumeric summary:")
print(summary_df, row.names = FALSE)
