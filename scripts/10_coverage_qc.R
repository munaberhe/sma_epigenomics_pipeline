.libPaths(c("/data/home/bt25018/R/library", .libPaths()))
suppressPackageStartupMessages({
    library(DMRcaller)
    library(GenomicRanges)
    library(ggplot2)
    library(patchwork)
})

# working directory
setwd("/data/scratch/bt25018/sma_epigenomics_pipeline")

# input / output paths
BY_CHR_DIR <- "results/alignments/bs/by_chr"
OUT_DIR    <- "results/qc/coverage_4lines"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# experiment design
CONDITIONS <- c("ASO_CTRL", "ASO_VPA", "Scramble_CTRL", "Scramble_VPA")
REPS       <- c(1, 2, 3)
CHROMS     <- paste0("chr", c(1:22, "X", "Y"))
XMAX       <- 100   # max coverage depth to plot

# colours: 3 blues for replicates, rust for pooled
TRACK_COLOURS <- c(
    rep1   = "#7CB6D6",
    rep2   = "#4A8FB8",
    rep3   = "#1B5478",
    pooled = "#A84B2F"
)

# read a single per-chr bismark CpG report into a GRanges object
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

# load all per-chr files for one replicate and combine into one GRanges
read_replicate_safe <- function(condition, rep, chroms, dir) {
    files <- file.path(dir, sprintf("%s_%d_%s.CpG_report.txt.gz", condition, rep, chroms))
    files <- files[file.exists(files)]
    if (length(files) == 0) stop("No files found for ", condition, "_", rep)
    message(sprintf("  %s rep %d: %d chr files", condition, rep, length(files)))
    grs <- lapply(files, read_cpg_report)
    do.call(c, unname(grs))
}

# for one condition: load 3 replicates, pool them, compute retention curves
build_retention_for_condition <- function(condition) {
    message("Reading replicates for ", condition, " ...")
    rep_data <- lapply(REPS, function(r) {
        read_replicate_safe(condition, r, CHROMS, BY_CHR_DIR)
    })
    names(rep_data) <- paste0("rep", REPS)

    # sanity check: all reps should have the same number of CpGs
    n_cpg <- sapply(rep_data, length)
    message(sprintf("  CpGs per replicate: %s", paste(n_cpg, collapse = ", ")))
    stopifnot(length(unique(n_cpg)) == 1)

    # pool replicates — sums readsM and readsN across reps at each CpG
    message("  Pooling replicates...")
    pooled <- poolMethylationDatasets(GRangesList(rep_data))

    # build a long data frame with coverage for all tracks
    make_df <- function(gr, label) {
        data.frame(coverage = gr$readsN, track = label, stringsAsFactors = FALSE)
    }
    df <- rbind(
        make_df(rep_data[[1]], "rep1"),
        make_df(rep_data[[2]], "rep2"),
        make_df(rep_data[[3]], "rep3"),
        make_df(pooled,        "pooled")
    )
    df$track     <- factor(df$track, levels = c("rep1","rep2","rep3","pooled"))
    df$condition <- condition
    df <- df[df$coverage > 0, ]

    # compute retention curve: fraction of CpGs >= each threshold
    thresholds <- 0:XMAX
    retention <- do.call(rbind, lapply(levels(df$track), function(t) {
        cov <- df$coverage[df$track == t]
        data.frame(
            condition = condition,
            track     = t,
            threshold = thresholds,
            retention = sapply(thresholds, function(k) mean(cov >= k))
        )
    }))
    retention$track <- factor(retention$track, levels = c("rep1","rep2","rep3","pooled"))

    # summary stats per track
    summarise_track <- function(label) {
        cov <- df$coverage[df$track == label]
        data.frame(
            condition  = condition,
            track      = label,
            n_cpg      = length(cov),
            mean       = round(mean(cov), 2),
            median     = median(cov),
            q90        = quantile(cov, 0.9, names = FALSE),
            pct_ge_10x = round(100 * mean(cov >= 10), 1)
        )
    }
    summary_df <- do.call(rbind, lapply(levels(df$track), summarise_track))

    list(retention = retention, summary = summary_df)
}

# run for all conditions
results <- lapply(CONDITIONS, build_retention_for_condition)
names(results) <- CONDITIONS

retention_all <- do.call(rbind, lapply(results, `[[`, "retention"))
summary_all   <- do.call(rbind, lapply(results, `[[`, "summary"))
retention_all$condition <- factor(retention_all$condition, levels = CONDITIONS)

# one panel per condition
make_panel <- function(cond) {
    df_c <- retention_all[retention_all$condition == cond, ]
    ggplot(df_c, aes(x = threshold, y = retention, colour = track)) +
        geom_line(linewidth = 0.9) +
        coord_cartesian(xlim = c(0, XMAX), ylim = c(0, 1)) +
        scale_colour_manual(values = TRACK_COLOURS) +
        geom_vline(xintercept = 10, linetype = 2, colour = "grey40") +
        annotate("text", x = 11, y = 0.97, label = "10x", hjust = 0,
                 colour = "grey40", size = 3.2) +
        labs(title  = cond,
             x      = "Coverage threshold (reads per CpG)",
             y      = "Fraction of CpGs >= threshold",
             colour = NULL) +
        theme_classic(base_size = 10) +
        theme(plot.title = element_text(face = "bold", size = 11),
              legend.position = "right",
              panel.grid.major.y = element_line(colour = "grey92"))
}

panels   <- lapply(CONDITIONS, make_panel)
combined <- (panels[[1]] / panels[[2]] / panels[[3]] / panels[[4]]) +
    plot_layout(guides = "collect") +
    plot_annotation(
        title    = "CpG coverage: replicates vs pooled, per condition",
        subtitle = "Pooled (rust) should sit above replicates (blues) — confirms pooling is summing not averaging",
        theme    = theme(plot.title    = element_text(face = "bold", size = 12),
                         plot.subtitle = element_text(size = 9, colour = "grey30"))
    ) & theme(legend.position = "right")

# save
ggsave(file.path(OUT_DIR, "coverage_4lines_per_condition.pdf"),
       combined, width = 8, height = 12)
ggsave(file.path(OUT_DIR, "coverage_4lines_per_condition.png"),
       combined, width = 8, height = 12, dpi = 200)
write.table(summary_all,
            file.path(OUT_DIR, "coverage_4lines_per_condition_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("\nSummary:")
print(summary_all, row.names = FALSE)
message("Done.")
