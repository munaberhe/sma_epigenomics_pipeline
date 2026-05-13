# ----------------------------------------------------------------------
# SMN1 / SMN2 locus methylation profiles (UNMASKED reference)
# BSmooth version - uses bsseq::BSmooth + bsseq::plotRegion
# ----------------------------------------------------------------------

suppressPackageStartupMessages({
    library(bsseq)
    library(GenomicRanges)
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

BSMOOTH_NS <- 70
BSMOOTH_H  <- 1000

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
    d <- d[d$context == "CG" & d$chr == "chr5" & d$pos >= 70044638 & d$pos <= 70958015, , drop = FALSE]
    d
}

build_pooled_bsseq <- function(condition) {
    files <- file.path(BY_CHR_DIR,
                       sprintf("%s_%d_chr5.CpG_report.txt.gz", condition, REPS))
    files <- files[file.exists(files)]
    if (length(files) != 3) warning(sprintf("%s: found %d/3 chr5 files", condition, length(files)))

    reps <- lapply(files, read_cpg_report)

    all_pos <- sort(unique(unlist(lapply(reps, function(d) paste0(d$chr, ":", d$pos)))))
    chrpos  <- do.call(rbind, strsplit(all_pos, ":", fixed = TRUE))
    pos_df  <- data.frame(chr = chrpos[,1], pos = as.integer(chrpos[,2]), stringsAsFactors = FALSE)
    pos_df  <- pos_df[order(pos_df$chr, pos_df$pos), ]
    key     <- paste0(pos_df$chr, ":", pos_df$pos)

    M_total   <- integer(nrow(pos_df))
    Cov_total <- integer(nrow(pos_df))
    for (d in reps) {
        idx <- match(paste0(d$chr, ":", d$pos), key)
        M_total[idx]   <- M_total[idx]   + d$countM
        Cov_total[idx] <- Cov_total[idx] + d$countM + d$countU
    }

    BSseq(chr   = pos_df$chr,
          pos   = pos_df$pos,
          M     = matrix(M_total,   ncol = 1),
          Cov   = matrix(Cov_total, ncol = 1),
          sampleNames = condition)
}

message("Building pooled BSseq objects per condition ...")
bs_per_cond <- lapply(CONDITIONS, build_pooled_bsseq)
names(bs_per_cond) <- CONDITIONS

message("Combining conditions into one BSseq object ...")
bs <- Reduce(combine, bs_per_cond)
pData(bs)$condition <- CONDITIONS
pData(bs)$col       <- COND_COLOURS[CONDITIONS]
pData(bs)$label     <- CONDITIONS
sampleNames(bs)     <- CONDITIONS

message(sprintf("BSseq object: %d CpGs x %d conditions", nrow(bs), ncol(bs)))

message(sprintf("Running BSmooth (ns=%d, h=%d) ...", BSMOOTH_NS, BSMOOTH_H))
bs.smooth <- BSmooth(bs,
                     ns       = BSMOOTH_NS,
                     h        = BSMOOTH_H,
                     
                     verbose  = TRUE)

plot_locus_bsmooth <- function(locus_name, out_pdf) {
    locus  <- LOCI[[locus_name]]
    region <- data.frame(chr   = locus$chr,
                         start = locus$start,
                         end   = locus$end)
    pdf(out_pdf, width = 10, height = 5)
    plotRegion(bs.smooth,
               region     = region,
               extend     = FLANK,
               main       = sprintf("%s locus (%s:%d-%d, +/- %d bp flank)",
                                    locus_name, locus$chr,
                                    locus$start, locus$end, FLANK),
               addRegions = NULL,
               annoTrack  = NULL,
               col        = pData(bs.smooth)$col,
               lty        = 1,
               lwd        = 2)
    legend("bottomleft",
           legend = CONDITIONS,
           col    = COND_COLOURS[CONDITIONS],
           lty = 1, lwd = 2, bty = "n", cex = 0.85)
    dev.off()
    message("Saved: ", out_pdf)
}

out_dir <- "results/qc/smn_locus"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

plot_locus_bsmooth("SMN1", file.path(out_dir, "smn_locus_profile_bsmooth_SMN1.pdf"))
plot_locus_bsmooth("SMN2", file.path(out_dir, "smn_locus_profile_bsmooth_SMN2.pdf"))

summarise_locus_raw <- function(condition_name, locus_name) {
    locus   <- LOCI[[locus_name]]
    bs.cond <- bs[, condition_name]
    gr      <- granges(bs.cond)
    sel     <- as.character(seqnames(gr)) == locus$chr &
               start(gr) >= locus$start & start(gr) <= locus$end
    M   <- getCoverage(bs.cond, type = "M")[sel, 1]
    Cov <- getCoverage(bs.cond, type = "Cov")[sel, 1]
    keep <- Cov > 0
    if (!any(keep)) return(data.frame(condition = condition_name, locus = locus_name,
                                       n_cpg = 0L, weighted_mean = NA_real_))
    data.frame(condition = condition_name, locus = locus_name,
               n_cpg = sum(keep),
               weighted_mean = sum(M[keep]) / sum(Cov[keep]))
}

summary_rows <- list()
for (cond in CONDITIONS) {
    for (locus_name in names(LOCI)) {
        summary_rows[[length(summary_rows) + 1]] <- summarise_locus_raw(cond, locus_name)
    }
}
summary_df <- do.call(rbind, summary_rows)

message("\nMean methylation in gene body (coverage-weighted, unmasked):")
print(summary_df, row.names = FALSE, digits = 4)

write.table(summary_df,
            file = file.path(out_dir, "smn_locus_profile_bsmooth_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Saved: ", file.path(out_dir, "smn_locus_profile_bsmooth_summary.tsv"))
