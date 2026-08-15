.libPaths(c("/data/home/bt25018/R/library", .libPaths()))
# ----------------------------------------------------------------------
# SMN1 / SMN2 locus methylation profiles - BSmooth Fig 1 style
#
# Companion to smn_locus_profile_bsmooth.R. Where that script produces
# the Hansen Fig 3 view (one line per condition, all conditions on the
# same axes), this script produces the Hansen Fig 1 view:
#
#   - 3 thin replicate lines + 1 thick mean line per condition
#   - 1 facet panel per condition (4 panels stacked)
#   - CpG rug at bottom showing where data exists
#   - Exon track at the bottom showing gene structure
#   - One PDF per locus
#
# Inputs:
#   results/alignments/bs/by_chr/<condition>_<rep>_chr5.CpG_report.txt.gz
# Outputs:
#   results/qc/smn_locus/smn_locus_profile_bsmooth_fig1_SMN1.pdf
#   results/qc/smn_locus/smn_locus_profile_bsmooth_fig1_SMN2.pdf
# ----------------------------------------------------------------------

suppressPackageStartupMessages({
    library(bsseq)
    library(GenomicRanges)
    library(BiocParallel)
    library(ggplot2)
    library(patchwork)
})

setwd("/data/scratch/bt25018/sma_epigenomics_pipeline")

# ---- Configuration -----------------------------------------------------
BY_CHR_DIR <- "results/alignments/bs/by_chr"

CONDITIONS <- c("ASO_CTRL", "ASO_VPA", "Scramble_CTRL", "Scramble_VPA")
REPS       <- 1:3

# hg38 SMN1 / SMN2 gene-body coordinates
LOCI <- list(
    SMN1 = list(chr = "chr5", start = 70924941, end = 70953015),
    SMN2 = list(chr = "chr5", start = 70049638, end = 70078522)
)
FLANK <- 5000

# BSmooth parameters (matched to Fig 3 script - keep them in sync)
BSMOOTH_NS <- 70
BSMOOTH_H  <- 1000

COND_COLOURS <- c(
    ASO_CTRL      = "#1B5478",
    ASO_VPA       = "#A84B2F",
    Scramble_CTRL = "#7CB6D6",
    Scramble_VPA  = "#D19900"
)

# ---- RefSeq exon coordinates for SMN1 / SMN2 (hg38) -------------------
# Hardcoded from UCSC RefSeq NM_000344 (SMN1) and NM_017411 (SMN2).
# 9 exons each. Coordinates are 1-based, inclusive.
EXONS <- list(
    SMN1 = data.frame(
        exon = 1:9,
        start = c(70924941, 70938861, 70941290, 70942358, 70944614,
                  70946066, 70946948, 70950965, 70951942),
        end   = c(70925545, 70938977, 70941389, 70942437, 70944696,
                  70946186, 70947057, 70951088, 70953015)
    ),
    SMN2 = data.frame(
        exon = 1:9,
        start = c(70049638, 70063488, 70065918, 70066986, 70069242,
                  70070697, 70071581, 70075600, 70076577),
        end   = c(70050242, 70063604, 70066017, 70067065, 70069324,
                  70070817, 70071690, 70075723, 70078522)
    )
)

# ---- Reader: CX report -> per-CpG counts -------------------------------
read_cpg_report <- function(f) {
    d <- read.table(gzfile(f), header = FALSE, sep = "\t",
                    col.names = c("chr","pos","strand","countM","countU","context","tri"),
                    colClasses = c("character","integer","character",
                                   "integer","integer","character","character"))
    d <- d[d$context == "CG", , drop = FALSE]
    d
}

# ---- Build BSseq with one column per replicate (NOT pooled) -----------
# Anchor on the union of all CpG positions across all 12 samples, then
# fill in M and Cov per sample.
message("Reading all 12 chr5 CpG reports ...")
all_files <- character(0)
sample_labels <- character(0)
sample_conditions <- character(0)
sample_reps <- integer(0)
for (cond in CONDITIONS) {
    for (rep in REPS) {
        f <- file.path(BY_CHR_DIR,
                       sprintf("%s_%d_chr5.CpG_report.txt.gz", cond, rep))
        if (!file.exists(f)) {
            warning("Missing: ", f); next
        }
        all_files <- c(all_files, f)
        sample_labels <- c(sample_labels, sprintf("%s_rep%d", cond, rep))
        sample_conditions <- c(sample_conditions, cond)
        sample_reps <- c(sample_reps, rep)
    }
}

reps_data <- lapply(all_files, read_cpg_report)
names(reps_data) <- sample_labels

# Union of CpG positions
all_keys <- sort(unique(unlist(lapply(reps_data, function(d) paste0(d$chr, ":", d$pos)))))
chrpos <- do.call(rbind, strsplit(all_keys, ":", fixed = TRUE))
pos_df <- data.frame(chr = chrpos[,1],
                     pos = as.integer(chrpos[,2]),
                     stringsAsFactors = FALSE)
pos_df <- pos_df[order(pos_df$chr, pos_df$pos), ]
key <- paste0(pos_df$chr, ":", pos_df$pos)

# Per-sample M and Cov matrices
M_mat   <- matrix(0L, nrow = length(key), ncol = length(reps_data))
Cov_mat <- matrix(0L, nrow = length(key), ncol = length(reps_data))
for (i in seq_along(reps_data)) {
    d <- reps_data[[i]]
    idx <- match(paste0(d$chr, ":", d$pos), key)
    M_mat[idx, i]   <- d$countM
    Cov_mat[idx, i] <- d$countM + d$countU
}
colnames(M_mat)   <- sample_labels
colnames(Cov_mat) <- sample_labels

bs <- BSseq(chr = pos_df$chr, pos = pos_df$pos,
            M = M_mat, Cov = Cov_mat,
            sampleNames = sample_labels)
pData(bs)$condition <- sample_conditions
pData(bs)$replicate <- sample_reps
pData(bs)$col       <- COND_COLOURS[sample_conditions]

message(sprintf("BSseq: %d CpGs x %d samples (12 = 4 conds x 3 reps)",
                nrow(bs), ncol(bs)))

# ---- Subset to the loci of interest (+ generous flank) -----------------
# Avoid OOM by only smoothing ~80 kb of chr5, not the whole 3.1M CpGs.
SMOOTH_FLANK <- 50000
keep_regions <- GRanges(
    seqnames = sapply(LOCI, function(l) l$chr),
    ranges   = IRanges(
        start = sapply(LOCI, function(l) l$start) - SMOOTH_FLANK,
        end   = sapply(LOCI, function(l) l$end)   + SMOOTH_FLANK
    )
)
bs <- subsetByOverlaps(bs, keep_regions)
message(sprintf("After subsetting to SMN loci +/- %d kb: %d CpGs",
                SMOOTH_FLANK / 1000, nrow(bs)))

# ---- BSmooth -----------------------------------------------------------
message(sprintf("Running BSmooth (ns=%d, h=%d) ...", BSMOOTH_NS, BSMOOTH_H))
bs.smooth <- BSmooth(bs,
                     ns = BSMOOTH_NS,
                     h  = BSMOOTH_H,
                     BPPARAM = SerialParam(),
                     verbose = TRUE)

# ---- Extract smoothed methylation in a region for plotting ------------
# Returns a long data.frame: position, methylation, sample, condition, replicate.
extract_smoothed <- function(bs.smooth, chr, start, end) {
    gr <- granges(bs.smooth)
    sel <- as.character(seqnames(gr)) == chr &
           start(gr) >= start &
           start(gr) <= end
    if (!any(sel)) return(NULL)

    meth <- getMeth(bs.smooth, type = "smooth")[sel, , drop = FALSE]
    pos  <- start(gr)[sel]

    do.call(rbind, lapply(seq_len(ncol(meth)), function(i) {
        data.frame(
            position    = pos,
            methylation = meth[, i],
            sample      = colnames(meth)[i],
            condition   = pData(bs.smooth)$condition[i],
            replicate   = pData(bs.smooth)$replicate[i],
            stringsAsFactors = FALSE
        )
    }))
}

# ---- Plot one locus, Fig 1 style --------------------------------------
plot_locus_fig1 <- function(locus_name, out_pdf) {
    locus  <- LOCI[[locus_name]]
    rstart <- locus$start - FLANK
    rend   <- locus$end   + FLANK

    df <- extract_smoothed(bs.smooth, locus$chr, rstart, rend)
    if (is.null(df) || nrow(df) == 0) {
        warning(locus_name, ": no smoothed data in region"); return(invisible())
    }
    df$condition <- factor(df$condition, levels = CONDITIONS)

    # Mean line per condition (mean of the 3 replicate smoothed tracks)
    mean_df <- aggregate(methylation ~ position + condition, df, FUN = mean)

    # ---- Per-condition rug ------------------------------------------------
    # For each condition, show CpG positions covered (>= 1 read) in at least
    # one of that condition's 3 replicates. Built from the raw Cov_mat (not
    # the smoothed object) so the rug reflects what was actually measured.
    raw_pos    <- start(bs)
    raw_chr    <- as.character(seqnames(bs))
    region_sel <- raw_chr == locus$chr & raw_pos >= rstart & raw_pos <= rend
    raw_pos_r  <- raw_pos[region_sel]
    raw_cov_r  <- getCoverage(bs)[region_sel, , drop = FALSE]

    rug_df <- do.call(rbind, lapply(CONDITIONS, function(cond) {
        cond_cols <- which(pData(bs)$condition == cond)
        if (length(cond_cols) == 0) return(NULL)
        any_cov   <- rowSums(raw_cov_r[, cond_cols, drop = FALSE] >= 1) >= 1
        if (!any(any_cov)) return(NULL)
        data.frame(position  = raw_pos_r[any_cov],
                   condition = cond,
                   stringsAsFactors = FALSE)
    }))
    rug_df$condition <- factor(rug_df$condition, levels = CONDITIONS)

    # ---- Methylation panel (Fig 1 main panel) -------------------------
    p_meth <- ggplot() +
        # gene body shading
        annotate("rect",
                 xmin = locus$start, xmax = locus$end,
                 ymin = -Inf, ymax = Inf,
                 fill = "grey85", alpha = 0.35) +
        # 3 thin replicate lines per condition
        geom_line(data = df,
                  aes(x = position, y = methylation,
                      group = sample, colour = condition),
                  linewidth = 0.35, alpha = 0.65) +
        # 1 thick mean line per condition
        geom_line(data = mean_df,
                  aes(x = position, y = methylation, colour = condition),
                  linewidth = 1.3) +
        # Per-condition CpG rug at bottom
        # (CpGs covered in at least one of that condition's 3 replicates)
        geom_rug(data = rug_df,
                 aes(x = position),
                 sides = "b", colour = "black", alpha = 0.4,
                 length = unit(0.02, "npc"),
                 inherit.aes = FALSE) +
        scale_colour_manual(values = COND_COLOURS, guide = "none") +
        scale_x_continuous(
            labels = function(x) sprintf("%.2f Mb", x / 1e6),
            breaks = scales::pretty_breaks(n = 6),
            limits = c(rstart, rend),
            expand = c(0, 0)
        ) +
        coord_cartesian(ylim = c(0, 1)) +
        facet_wrap(~ condition, ncol = 1, strip.position = "right") +
        labs(
            title = sprintf("%s locus (%s:%d-%d, +/- %d bp flank)",
                            locus_name, locus$chr,
                            locus$start, locus$end, FLANK),
            subtitle = "Thin lines = 3 replicates. Thick line = mean. Rug = CpGs covered in at least one of that condition's 3 replicates.",
            x = NULL,
            y = "Methylation"
        ) +
        theme_classic(base_size = 11) +
        theme(
            plot.title = element_text(face = "bold", size = 12),
            plot.subtitle = element_text(size = 9, colour = "grey30"),
            strip.text.y = element_text(angle = 0, face = "bold", size = 9),
            strip.background = element_rect(fill = "grey95", colour = NA),
            panel.spacing = unit(0.4, "lines"),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        )

    # ---- Exon track panel (bottom) ------------------------------------
    exons <- EXONS[[locus_name]]
    p_exon <- ggplot() +
        # intron line
        annotate("segment",
                 x = locus$start, xend = locus$end,
                 y = 0.5, yend = 0.5,
                 colour = "grey40", linewidth = 0.4) +
        # exon boxes
        geom_rect(data = exons,
                  aes(xmin = start, xmax = end, ymin = 0.15, ymax = 0.85),
                  fill = "steelblue4", colour = "steelblue4") +
        # exon numbers above each box
        geom_text(data = exons,
                  aes(x = (start + end) / 2, y = 1.1, label = exon),
                  size = 2.5, colour = "grey20") +
        scale_x_continuous(
            labels = function(x) sprintf("%.2f Mb", x / 1e6),
            breaks = scales::pretty_breaks(n = 6),
            limits = c(rstart, rend),
            expand = c(0, 0)
        ) +
        coord_cartesian(ylim = c(0, 1.4)) +
        labs(x = sprintf("Position (%s)", locus$chr),
             y = NULL,
             caption = sprintf("%s exon track (RefSeq, hg38). Boxes = exons, numbered 1-9. Line = introns.",
                               locus_name)) +
        theme_classic(base_size = 10) +
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            plot.caption = element_text(hjust = 0, size = 8, colour = "grey40")
        )

    combined <- p_meth / p_exon +
        plot_layout(heights = c(7, 1.2))

    ggsave(out_pdf, combined, width = 9, height = 9)
    message("Saved: ", out_pdf)
}

out_dir <- "results/qc/smn_locus"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

plot_locus_fig1("SMN1", file.path(out_dir, "smn_locus_profile_bsmooth_fig1_SMN1.pdf"))
plot_locus_fig1("SMN2", file.path(out_dir, "smn_locus_profile_bsmooth_fig1_SMN2.pdf"))

message("\nDone. Two PDFs saved to ", out_dir)
