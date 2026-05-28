suppressPackageStartupMessages({
    library(DMRcaller)
    library(GenomicRanges)
    library(ggplot2)
    library(patchwork)
})

setwd("/data/scratch/bt25018/sma_epigenomics_pipeline")

BY_CHR_DIR <- "results/alignments/bs/by_chr"
CONDITIONS <- c("ASO_CTRL", "ASO_VPA", "Scramble_CTRL", "Scramble_VPA")
REPS       <- 1:3

LOCI <- list(
    SMN1 = list(chr = "chr5", start = 70924941, end = 70953015),
    SMN2 = list(chr = "chr5", start = 70049638, end = 70078522)
)

FLANK <- 5000

COND_COLOURS <- c(
    ASO_CTRL      = "#1B5478",
    ASO_VPA       = "#A84B2F",
    Scramble_CTRL = "#7CB6D6",
    Scramble_VPA  = "#D19900"
)

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
    files <- file.path(BY_CHR_DIR, sprintf("%s_%d_chr5.CpG_report.txt.gz", condition, REPS))
    files <- files[file.exists(files)]
    if (length(files) != 3) warning(sprintf("%s: found %d/3 chr5 files", condition, length(files)))
    grs <- lapply(files, read_cpg_report)
    poolMethylationDatasets(GRangesList(grs))
}

message("Reading and pooling chr5 by condition ...")
pooled_data <- lapply(CONDITIONS, read_condition)
names(pooled_data) <- CONDITIONS

extract_locus_df <- function(condition_name, locus_name) {
    locus <- LOCI[[locus_name]]
    gr  <- pooled_data[[condition_name]]
    sel <- as.character(seqnames(gr)) == locus$chr &
           start(gr) >= (locus$start - FLANK) &
           start(gr) <= (locus$end   + FLANK) &
           gr$context == "CG" &
           gr$readsN  > 0
    g <- gr[sel]
    if (length(g) == 0) return(NULL)
    data.frame(
        pos         = start(g),
        methylation = g$readsM / g$readsN,
        coverage    = g$readsN,
        condition   = condition_name,
        locus       = locus_name
    )
}

message("Extracting per-CpG data ...")
cpg_df <- do.call(rbind, lapply(CONDITIONS, function(cond) {
    do.call(rbind, lapply(names(LOCI), function(loc) {
        extract_locus_df(cond, loc)
    }))
}))
cpg_df$condition <- factor(cpg_df$condition, levels = CONDITIONS)


compute_loess_clamped <- function(df, span = 0.75) {
    do.call(rbind, lapply(split(df, df$condition), function(sub) {
        if (nrow(sub) < 4) return(NULL)
        fit <- loess(methylation ~ pos, data = sub, weights = sub$coverage,
                     span = span, na.action = na.exclude
                     )
        xs  <- seq(min(sub$pos), max(sub$pos), length.out = 400)
        ys  <- predict(fit, newdata = data.frame(pos = xs))
        ys  <- pmin(pmax(ys, 0), 0.99)
        data.frame(pos = xs, methylation = ys, condition = sub$condition[1])
    }))
}

plot_locus_scatter <- function(df_all, locus_name) {
    locus <- LOCI[[locus_name]]
    df <- df_all[df_all$locus == locus_name, ]
    ggplot(df, aes(x = pos, y = methylation, colour = condition)) +
        geom_point(aes(size = coverage), alpha = 0.45, shape = 16) +
        geom_line(data = compute_loess_clamped(df, span = 0.75),
                  aes(x = pos, y = methylation, colour = condition),
                  linewidth = 1.1, alpha = 0.85, inherit.aes = FALSE) +
        annotate("rect", xmin = locus$start, xmax = locus$end,
                 ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.20) +
        annotate("text", x = (locus$start + locus$end) / 2, y = 1.04,
                 label = sprintf("%s gene body (%.1f kb)", locus_name,
                                 (locus$end - locus$start) / 1000),
                 size = 3.5, colour = "grey30", vjust = 0) +
        scale_colour_manual(values = COND_COLOURS) +
        scale_size_continuous(range = c(0.5, 4), guide = "none") +
        scale_x_continuous(labels = function(x) sprintf("%.2f Mb", x / 1e6),
                           breaks = scales::pretty_breaks(n = 5)) +
        coord_cartesian(ylim = c(0, 1), clip = "off") +
        labs(title = sprintf("%s locus (%s:%d-%d, +/- %d bp flank)",
                             locus_name, locus$chr, locus$start, locus$end, FLANK),
             x = "Position (chr5)", y = "CpG methylation", colour = NULL) +
        theme_classic(base_size = 11) +
        theme(plot.margin = margin(t = 18, r = 5, b = 5, l = 5)) +
        theme(legend.position = "right",
              panel.grid.major.y = element_line(colour = "grey92"),
              plot.title = element_text(face = "bold", size = 11))
}

p1 <- plot_locus_scatter(cpg_df, "SMN1")
p2 <- plot_locus_scatter(cpg_df, "SMN2")

combined <- p1 / p2 +
    plot_layout(guides = "collect") +
    plot_annotation(
        title    = "SMN locus methylation profiles (UNMASKED hg38 alignment)",
        subtitle = "Points = individual CpGs (size = coverage). Lines = LOESS smoother (coverage-weighted, span=0.75).",
        theme    = theme(plot.title = element_text(face = "bold", size = 13))
    ) & theme(legend.position = "right")

out_dir <- "results/qc/smn_locus"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(out_dir, "smn_locus_profile_unmasked.png"), combined, width = 10, height = 7, dpi = 200)
ggsave(file.path(out_dir, "smn_locus_profile_unmasked.pdf"), combined, width = 10, height = 7)
message("Saved: results/qc/smn_locus/smn_locus_profile_unmasked.pdf")

summarise_locus_raw <- function(condition_name, locus_name) {
    locus <- LOCI[[locus_name]]
    gr  <- pooled_data[[condition_name]]
    sel <- as.character(seqnames(gr)) == locus$chr &
           start(gr) >= locus$start & start(gr) <= locus$end &
           gr$context == "CG" & gr$readsN > 0
    g <- gr[sel]
    if (length(g) == 0) return(data.frame(condition = condition_name, locus = locus_name,
                                           n_cpg = 0L, weighted_mean = NA_real_,
                                           unweighted_mean = NA_real_))
    data.frame(condition = condition_name, locus = locus_name,
               n_cpg = length(g),
               weighted_mean   = sum(g$readsM) / sum(g$readsN),
               unweighted_mean = mean(g$readsM / g$readsN))
}

summary_rows <- list()
for (cond in CONDITIONS) {
    for (locus_name in names(LOCI)) {
        summary_rows[[length(summary_rows) + 1]] <- summarise_locus_raw(cond, locus_name)
    }
}
summary_df <- do.call(rbind, summary_rows)
message("\nMean methylation in gene body (per-CpG, unmasked):")
print(summary_df, row.names = FALSE, digits = 4)
