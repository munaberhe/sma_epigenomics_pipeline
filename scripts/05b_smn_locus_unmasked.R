# ----------------------------------------------------------------------
# SMN1 / SMN2 local methylation profile - DMRcaller plotLocalMethylationProfile
#
# Radu's follow-up:
#   "Please use the function in DMRcaller to plot that so we can also see
#    the exons/introns to see if we get anything over the 7th exon."
#
# EXON NUMBERING NOTE (Alberto Kornblihtt, 15 May 2026):
#   What the GTF labels as exon 2 is actually two exons: E2a and E2b.
#   This shifts Alberto's numbering vs the GTF:
#     GTF exon 1  = Alberto E1
#     GTF exon 2  = Alberto E2a
#     GTF exon 3  = Alberto E2b
#     GTF exon 4  = Alberto E3
#     GTF exon 5  = Alberto E4
#     GTF exon 6  = Alberto E5
#     GTF exon 7  = Alberto E6
#     GTF exon 8  = Alberto E7  ← THE ASO TARGET (the penultimate exon)
#     GTF exon 9  = Alberto E8
#
#   The plot labels exons using Alberto's convention (E1, E2a, E2b, E3...)
#   and marks E7 in red as the therapeutic target.
#
# Three comparisons (each isolates one variable):
#   1. ASO_CTRL  vs Scramble_CTRL   (effect of ASO alone)
#   2. ASO_CTRL  vs ASO_VPA         (effect of VPA in ASO background)
#   3. ASO_VPA   vs Scramble_VPA    (effect of ASO under VPA)
# ----------------------------------------------------------------------

suppressPackageStartupMessages({
    library(DMRcaller)
    library(GenomicRanges)
})

setwd("/data/scratch/bt25018/sma_epigenomics_pipeline")

# ---- Configuration -----------------------------------------------------
BY_CHR_DIR <- "results/alignments/bs/by_chr"
REPS       <- 1:3
FLANK      <- 5000
WIN_SIZE   <- 300

# hg38 SMN1 / SMN2 gene-body coordinates
LOCI <- list(
    SMN1 = list(chr = "chr5", start = 70924941, end = 70953015, strand = "+"),
    SMN2 = list(chr = "chr5", start = 70049638, end = 70078522, strand = "+")
)

# Exon coordinates — labelled using Alberto's convention
# SMN2: GTF exon 8 (70,076,521-70,076,574) = Alberto's E7 = ASO target
EXONS <- list(
    SMN1 = data.frame(
        exon  = c("E1","E2a","E2b","E3","E4","E5","E6","E7","E8"),
        start = c(70925030, 70938807, 70941357, 70942326, 70942686,
                  70944627, 70946033, 70951913, 70952411),
        end   = c(70925158, 70938878, 70941476, 70942526, 70942838,
                  70944722, 70946143, 70951966, 70952984),
        # E7 is the ASO target in SMN1 (same penultimate exon)
        is_target = c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE)
    ),
    SMN2 = data.frame(
        exon  = c("E1","E2a","E2b","E3","E4","E5","E6","E7","E8"),
        start = c(70049638, 70063415, 70065965, 70066934, 70067294,
                  70069235, 70070641, 70076521, 70077019),
        end   = c(70049766, 70063486, 70066084, 70067134, 70067446,
                  70069330, 70070751, 70076574, 70077592),
        # E7 = GTF exon 8 = Alberto's ASO target (penultimate exon)
        is_target = c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE)
    )
)

# Three comparisons
COMPARISONS <- list(
    list(name = "ASO_vs_Scramble_CTRL",
         cond1 = "ASO_CTRL",  cond2 = "Scramble_CTRL",
         label = "ASO effect (CTRL background)"),
    list(name = "VPA_vs_CTRL_ASO",
         cond1 = "ASO_CTRL",  cond2 = "ASO_VPA",
         label = "VPA effect (ASO background)"),
    list(name = "ASO_vs_Scramble_VPA",
         cond1 = "ASO_VPA",   cond2 = "Scramble_VPA",
         label = "ASO effect (VPA background)")
)

NEEDED_CONDS <- unique(unlist(lapply(COMPARISONS, function(x) c(x$cond1, x$cond2))))

# ---- Reader ------------------------------------------------------------
read_cpg_report <- function(f) {
    d <- read.table(gzfile(f), header = FALSE, sep = "\t",
                    col.names = c("chr","pos","strand","countM","countU","context","tri"),
                    colClasses = c("character","integer","character",
                                   "integer","integer","character","character"))
    GRanges(
        seqnames = d$chr,
        ranges   = IRanges(d$pos, d$pos),
        strand   = d$strand,
        readsM   = d$countM,
        readsN   = d$countM + d$countU,
        context  = d$context,
        trinucleotide_context = d$tri
    )
}

read_condition <- function(condition) {
    files <- file.path(BY_CHR_DIR,
                       sprintf("%s_%d_chr5.CpG_report.txt.gz", condition, REPS))
    files <- files[file.exists(files)]
    if (length(files) != 3)
        warning(sprintf("%s: found %d/3 chr5 files", condition, length(files)))
    grs <- lapply(files, read_cpg_report)
    poolMethylationDatasets(GRangesList(grs))
}

message("Reading and pooling chr5 by condition ...")
pooled_data <- lapply(NEEDED_CONDS, function(c) { message("  ", c); read_condition(c) })
names(pooled_data) <- NEEDED_CONDS

# ---- Build GFF with Alberto's exon labels ------------------------------
build_gff <- function() {
    rows <- list()
    for (locus_name in names(LOCI)) {
        locus <- LOCI[[locus_name]]
        rows[[length(rows)+1]] <- data.frame(
            chr = locus$chr, start = locus$start, end = locus$end,
            strand = locus$strand, type = "gene", name = locus_name,
            stringsAsFactors = FALSE
        )
        ex <- EXONS[[locus_name]]
        for (i in seq_len(nrow(ex))) {
            rows[[length(rows)+1]] <- data.frame(
                chr = locus$chr, start = ex$start[i], end = ex$end[i],
                strand = locus$strand, type = "exon",
                name = sprintf("%s_%s", locus_name, ex$exon[i]),
                stringsAsFactors = FALSE
            )
        }
    }
    df <- do.call(rbind, rows)
    GRanges(seqnames = df$chr, ranges = IRanges(df$start, df$end),
            strand = df$strand, type = df$type, name = df$name)
}

GEs <- build_gff()
message("Built GFF: ", length(GEs), " features")

# ---- Plot function with exon labels and E7 highlighted -----------------
plot_one <- function(comp, locus_name) {
    locus  <- LOCI[[locus_name]]
    region <- GRanges(seqnames = locus$chr,
                      ranges = IRanges(locus$start - FLANK, locus$end + FLANK))

    plotLocalMethylationProfile(
        methylationData1 = pooled_data[[comp$cond1]],
        methylationData2 = pooled_data[[comp$cond2]],
        region           = region,
        DMRs             = NULL,
        conditionsNames  = c(comp$cond1, comp$cond2),
        gff              = GEs,
        windowSize       = WIN_SIZE,
        context          = "CG",
        main             = sprintf("%s — %s vs %s (%s)",
                                   locus_name, comp$cond1, comp$cond2, comp$label),
        plotMeanLines    = TRUE,
        plotPoints       = TRUE
    )

    # Add exon labels using Alberto's convention (E1, E2a, E2b, E3...E7 in red)
    # Labels placed below the exon track using text() at y just below the plot area
    ex  <- EXONS[[locus_name]]
    usr <- par("usr")
    plt <- par("plt")
    fig <- par("fig")

    for (i in seq_len(nrow(ex))) {
        x_mid  <- (ex$start[i] + ex$end[i]) / 2
        colour <- if (ex$is_target[i]) "red" else "black"
        weight <- if (ex$is_target[i]) 2 else 1
        # Place label in the lower margin using text() in user coordinates
        # y_label sits just below the plot lower boundary (usr[3])
        
        mtext(ex$exon[i], side = 1, at = x_mid, line = -1.5,
              cex = 0.45, col = colour, font = weight)
    }


    # Gene name top-right
    text(LOCI[[locus_name]]$start - FLANK * 0.8, usr[3] + (usr[4] - usr[3]) * 0.05, labels = locus_name, cex = 0.8, font = 2, col = "black", adj = c(1, 0.5))
}

# ---- Save PDFs ---------------------------------------------------------
out_dir <- "results/qc/smn_locus"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (comp in COMPARISONS) {
    pdf_path <- file.path(out_dir,
                          sprintf("smn_locus_dmrcaller_%s.pdf", comp$name))
    message("\nPlotting: ", comp$name)
    pdf(pdf_path, width = 11, height = 8.5, bg="white")
    par(mfrow = c(2,1), mar = c(5,4,3,1)+0.1, cex = 0.9, bg="white", col.axis="black", col.lab="black", col.main="black", fg="black")
    plot_one(comp, "SMN1")
    plot_one(comp, "SMN2")
    dev.off()
    message("  Saved: ", pdf_path)
}

# Combined PDF — all 3 comparisons, width=16 for Radu
combined_pdf <- file.path(out_dir, "smn_locus_dmrcaller_all_comparisons.pdf")
message("\nPlotting combined: ", combined_pdf)
pdf(combined_pdf, width = 16, height = 14, bg="white")
par(mfrow = c(3,2), mar = c(5,4,3,1)+0.1, cex = 0.7, bg="white", col.axis="black", col.lab="black", col.main="black", fg="black")
for (comp in COMPARISONS) {
    plot_one(comp, "SMN1")
    plot_one(comp, "SMN2")
}
dev.off()
message("  Saved.")

# ---- Numeric summary ---------------------------------------------------
summarise <- function(condition_name, locus_name) {
    locus <- LOCI[[locus_name]]
    gr    <- pooled_data[[condition_name]]
    sel   <- as.character(seqnames(gr)) == locus$chr &
             start(gr) >= locus$start & start(gr) <= locus$end &
             gr$context == "CpG" & gr$readsN > 0
    g <- gr[sel]
    if (length(g) == 0)
        return(data.frame(condition=condition_name, locus=locus_name,
                          n_cpg=0L, weighted_mean=NA_real_))
    data.frame(condition=condition_name, locus=locus_name,
               n_cpg=length(g),
               weighted_mean=round(sum(g$readsM)/sum(g$readsN), 4))
}

summary_rows <- list()
for (cond in NEEDED_CONDS)
    for (locus_name in names(LOCI))
        summary_rows[[length(summary_rows)+1]] <- summarise(cond, locus_name)

summary_df <- do.call(rbind, summary_rows)
tsv_path   <- file.path(out_dir, "smn_locus_dmrcaller_comparisons_summary.tsv")
write.table(summary_df, tsv_path, sep="\t", quote=FALSE, row.names=FALSE)
message("\nSaved: ", tsv_path)
print(summary_df, row.names=FALSE)
message("\nDone. Outputs in ", out_dir)
