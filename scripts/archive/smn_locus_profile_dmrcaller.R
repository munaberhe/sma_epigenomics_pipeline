suppressPackageStartupMessages({
    library(DMRcaller)
    library(GenomicRanges)
    library(ggplot2)
    library(patchwork)
})

.libPaths(c("/data/home/bt25018/R/library", .libPaths()))

setwd("/data/scratch/bt25018/sma_epigenomics_pipeline")

BY_CHR_DIR <- "results/alignments/bs/by_chr"
CONDITIONS <- c("ASO_CTRL", "ASO_VPA", "Scramble_CTRL", "Scramble_VPA")
REPS       <- 1:3

LOCI <- list(
    SMN1 = list(chr = "chr5", start = 70924941, end = 70953015),
    SMN2 = list(chr = "chr5", start = 70049638, end = 70078522)
)
FLANK <- 5000

WIN_SIZE <- 1500
WIN_STEP <- 200

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

profile_in_locus <- function(gr, locus, win_size = WIN_SIZE, win_step = WIN_STEP) {
    chr   <- locus$chr
    start <- locus$start - FLANK
    end   <- locus$end   + FLANK

    win_starts <- seq(start, end - win_size, by = win_step)
    do.call(rbind, lapply(win_starts, function(ws) {
        we <- ws + win_size
        sub_region <- GRanges(seqnames = chr, ranges = IRanges(ws, we))
        p <- computeMethylationProfile(
            methylationData = gr,
            region          = sub_region,
            windowSize      = win_size,
            context         = "CG"
        )
        if (length(p) == 0) {
            if (length(p) == 0 || sum(p$readsN) < 3) {
            data.frame(midpoint = ws + win_size/2, methylation = NA_real_)
        } else {
            data.frame(midpoint = ws + win_size/2, methylation = p$Proportion[1])
        }
        #REMOVE
        } else {
            data.frame(midpoint = ws + win_size/2, methylation = p$Proportion[1])
        }
    }))
}

build_profile_df <- function(locus_name) {
    locus <- LOCI[[locus_name]]
    do.call(rbind, lapply(CONDITIONS, function(cond) {
        df <- profile_in_locus(pooled_data[[cond]], locus)
        df$condition <- cond
        df$locus     <- locus_name
        df
    }))
}

message("Computing DMRcaller profiles for SMN1 ...")
smn1_df <- build_profile_df("SMN1")
message("Computing DMRcaller profiles for SMN2 ...")
smn2_df <- build_profile_df("SMN2")

plot_locus <- function(df, locus_name) {
    locus <- LOCI[[locus_name]]
    df$condition <- factor(df$condition, levels = CONDITIONS)
    ggplot(df, aes(x = midpoint, y = methylation, colour = condition)) +
        annotate("rect", xmin = locus$start, xmax = locus$end,
                 ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.35) +
        geom_line(linewidth = 0.8, na.rm = TRUE) +
        scale_colour_manual(values = COND_COLOURS) +
        scale_x_continuous(labels = function(x) sprintf("%.2f Mb", x / 1e6),
                           breaks = scales::pretty_breaks(n = 6)) +
        coord_cartesian(ylim = c(0, 1)) +
        labs(title = sprintf("%s locus (%s:%d-%d, +/- %d bp flank)",
                             locus_name, locus$chr, locus$start, locus$end, FLANK),
             x = sprintf("Position (%s)", locus$chr),
             y = "CpG methylation", colour = NULL) +
        theme_classic(base_size = 11) +
        theme(plot.title = element_text(face = "bold", size = 11),
              legend.position = "right",
              panel.grid.major.y = element_line(colour = "grey92"))
}

p1 <- plot_locus(smn1_df, "SMN1")
p2 <- plot_locus(smn2_df, "SMN2")

combined <- p1 / p2 +
    plot_layout(guides = "collect") +
    plot_annotation(
        title    = "SMN locus methylation profiles - DMRcaller sliding window (UNMASKED hg38)",
        subtitle = sprintf("Pooled across 3 replicates per condition. Window = %d bp, step = %d bp. Shaded = gene body.",
                           WIN_SIZE, WIN_STEP),
        theme = theme(plot.title = element_text(face = "bold", size = 12),
                      plot.subtitle = element_text(size = 9, colour = "grey30"))
    ) & theme(legend.position = "right")

out_dir <- "results/qc/smn_locus"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(out_dir, "smn_locus_profile_dmrcaller.pdf"), combined, width = 10, height = 7)
ggsave(file.path(out_dir, "smn_locus_profile_dmrcaller.png"), combined, width = 10, height = 7, dpi = 200)
message("Saved: results/qc/smn_locus/smn_locus_profile_dmrcaller.pdf")

summarise_locus <- function(condition_name, locus_name) {
    locus <- LOCI[[locus_name]]
    gr  <- pooled_data[[condition_name]]
    sel <- as.character(seqnames(gr)) == locus$chr &
           start(gr) >= locus$start & start(gr) <= locus$end &
           gr$context == "CG" & gr$readsN > 0
    g <- gr[sel]
    if (length(g) == 0) return(data.frame(condition = condition_name, locus = locus_name,
                                           n_cpg = 0L, weighted_mean = NA_real_))
    data.frame(condition = condition_name, locus = locus_name,
               n_cpg = length(g),
               weighted_mean = sum(g$readsM) / sum(g$readsN))
}

summary_rows <- list()
for (cond in CONDITIONS) {
    for (locus_name in names(LOCI)) {
        summary_rows[[length(summary_rows) + 1]] <- summarise_locus(cond, locus_name)
    }
}
summary_df <- do.call(rbind, summary_rows)
message("\nMean methylation in gene body (coverage-weighted, unmasked):")
print(summary_df, row.names = FALSE, digits = 4)

write.table(summary_df,
            file = file.path(out_dir, "smn_locus_profile_dmrcaller_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Saved: results/qc/smn_locus/smn_locus_profile_dmrcaller_summary.tsv")
