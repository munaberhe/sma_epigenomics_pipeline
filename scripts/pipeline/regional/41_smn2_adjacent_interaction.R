## 41_smn2_adjacent_interaction.R
## Hypothesis-driven test of the SMN2-adjacent flanking region.
##
## Why this locus. It was not selected by ranking. The region was examined
## because SMN2 is the target of the antisense oligonucleotide, so a
## treatment-associated methylation change in its flanking sequence is a
## prediction rather than a discovery. That makes it the one candidate region in
## the thesis that is not subject to selection bias, which matters when
## interpreting a p-value.
##
## Region under test: chr5:70,088,223-70,088,522 (300 bp), reported at 35.4
## percentage points of methylation loss, from 94.3 percent to 58.9 percent,
## approximately 9.7 kb beyond the SMN2 3' end. The single-region permutation
## p-value of 0.43 came from a label-based null with n=3 per condition, which
## has a coarse resolution floor. This script does three things that null cannot:
##
##   1. Reports the raw methylated and unmethylated read counts behind every
##      percentage, so the reader can see how much data supports the estimate.
##   2. Fits a binomial GLM with an ASO by VPA interaction term on those counts,
##      giving a directional test of the 2x2 design rather than a rank statistic.
##   3. Refits with quasibinomial and beta-binomial dispersion, because pooled
##      read counts are overdispersed and a plain binomial will understate the
##      standard error. The honest p-value is the overdispersed one.
##
## Reporting rule. If the interaction term is not significant under the
## overdispersed model, say so. The value of this analysis is that it is a
## pre-specified test with an interpretable effect size, not that it is positive.
##
## Usage:
##   Rscript 41_smn2_adjacent_interaction.R
##   DISCOVER_ONLY=TRUE Rscript 41_smn2_adjacent_interaction.R

script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(dirname(normalizePath(sub("^--file=", "", m[1]))))
  if (!is.null(sys.frames()[[1]]$ofile))
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  getwd()
}
source(file.path(script_dir(), "00_regional_common.R"))

set.seed(SEED)
TAG <- "41_smn2_adjacent"
DISCOVER_ONLY <- toupper(Sys.getenv("DISCOVER_ONLY", "FALSE")) == "TRUE"

REGIONS <- data.table(
  name  = c("SMN2 flanking DMR",
            "relaxed scan 1", "relaxed scan 2", "relaxed scan 3", "relaxed scan 4"),
  chr   = "chr5",
  start = c(70088223L, 70046438L, 70074838L, 70079338L, 70088438L),
  end   = c(70088522L, 70046737L, 70075137L, 70079637L, 70088487L),
  primary = c(TRUE, FALSE, FALSE, FALSE, FALSE)
)
PRIMARY <- REGIONS[primary == TRUE]

cat("== 41_smn2_adjacent_interaction ==\n")
cx <- list.files(CX_DIR, full.names = TRUE, recursive = TRUE,
                 pattern = "\\.(txt|tsv|cov|CX_report|bedGraph)(\\.gz)?$")
cat("Files under ", CX_DIR, ": ", length(cx), "\n", sep = "")
if (length(cx)) cat(paste0("  ", basename(cx), "\n"), sep = "")

if (DISCOVER_ONLY) quit(save = "no")
if (!length(cx)) stop("No per-cytosine files found under ", CX_DIR,
                      ". This script needs read-level counts, not called DMRs.")

read_cx_region <- function(path, regions) {
  x <- fread(path, header = FALSE, showProgress = FALSE)
  if (ncol(x) >= 6 && all(x$V3[1:min(50, nrow(x))] >= x$V2[1:min(50, nrow(x))]) &&
      is.numeric(x$V4)) {
    setnames(x, 1:6, c("chr", "pos", "end", "pct", "nM", "nU"))
    x[, c("end", "pct") := NULL]
    x[, strand := "*"]
  } else {
    setnames(x, 1:5, c("chr", "pos", "strand", "nM", "nU"))
    if (ncol(x) >= 6) x <- x[V6 == "CG"]
  }
  x[, chr := as.character(chr)]
  if (!any(grepl("^chr", x$chr))) x[, chr := paste0("chr", chr)]
  out <- rbindlist(lapply(seq_len(nrow(regions)), function(i) {
    r <- regions[i]
    s <- x[chr == r$chr & pos >= r$start & pos <= r$end]
    if (!nrow(s)) return(NULL)
    s[, region := r$name]
    s[]
  }))
  out
}

label_condition <- function(fn) {
  for (i in seq_len(nrow(COND)))
    if (grepl(COND$pattern[i], fn, ignore.case = TRUE)) return(COND$condition[i])
  NA_character_
}

cat("\nReading per-cytosine counts in target regions:\n")
per_cyt <- rbindlist(lapply(cx, function(p) {
  cond <- label_condition(basename(p))
  if (is.na(cond)) {
    cat("  skipping (condition not matched):", basename(p), "\n"); return(NULL)
  }
  d <- read_cx_region(p, REGIONS)
  if (is.null(d) || !nrow(d)) {
    cat("  ", basename(p), " -> ", cond, ": 0 covered CpGs in regions\n", sep = "")
    return(NULL)
  }
  d[, `:=`(file = basename(p), condition = cond,
           replicate = sub(".*?([0-9]+)\\D*$", "\\1", tools::file_path_sans_ext(basename(p))))]
  cat("  ", basename(p), " -> ", cond, ": ",
      nrow(d), " covered CpGs, ", sum(d$nM + d$nU), " reads\n", sep = "")
  d[]
}), fill = TRUE)

if (!nrow(per_cyt))
  stop("No covered CpGs recovered in any target region. Check chromosome naming ",
       "and that CX_DIR holds the pre-deduplication SMN2 extraction.")

fwrite(per_cyt, file.path(OUT_DIR, "smn2_adjacent_per_cytosine.csv"))

counts <- per_cyt[, .(n_cpg = .N,
                      meth = sum(nM), unmeth = sum(nU),
                      depth = sum(nM + nU)),
                  by = .(region, condition)]
counts[, pct_meth := 100 * meth / (meth + unmeth)]
counts <- merge(counts, COND[, .(condition, aso, vpa)], by = "condition")
counts[, condition := factor(condition, levels = COND$condition)]
setorder(counts, region, condition)

cat("\n== Raw counts by region and condition ==\n")
print(counts[, .(region, condition, n_cpg, meth, unmeth, depth,
                 pct_meth = round(pct_meth, 1))])
fwrite(counts, file.path(OUT_DIR, "smn2_adjacent_counts.csv"))

cat("\nMean read depth per covered CpG, primary region:\n")
print(counts[region == PRIMARY$name,
             .(condition, depth_per_cpg = round(depth / n_cpg, 2))])

fit_region <- function(cd, region_name) {
  d <- cd[region == region_name]
  if (nrow(d) < 4) {
    cat("\n[", region_name, "] fewer than four conditions covered; ",
        "no interaction test possible.\n", sep = "")
    return(NULL)
  }
  y <- cbind(meth = d$meth, unmeth = d$unmeth)
  cat("\n=====================================================\n")
  cat("Region: ", region_name, "  (", d$chr[1] %||% "chr5", ")\n", sep = "")
  cat("=====================================================\n")

  cat("\n-- Binomial GLM, ASO * VPA (understates SE, shown for reference) --\n")
  m_bin <- glm(y ~ aso * vpa, data = d, family = binomial())
  print(summary(m_bin)$coefficients)

  disp <- sum(residuals(m_bin, type = "pearson")^2) / df.residual(m_bin)
  cat("\nResidual overdispersion (Pearson chi-square / df): ",
      round(disp, 2), "\n", sep = "")
  if (is.finite(disp) && disp > 1.5)
    cat("Overdispersed. The binomial p-values above are anti-conservative and\n",
        "must not be reported. Use the quasibinomial fit below.\n", sep = "")

  cat("\n-- Quasibinomial GLM, ASO * VPA (report this one) --\n")
  m_q <- glm(y ~ aso * vpa, data = d, family = quasibinomial())
  print(summary(m_q)$coefficients)

  ci <- try(suppressMessages(confint(m_q)), silent = TRUE)
  if (!inherits(ci, "try-error")) {
    cat("\nProfile confidence intervals (log odds):\n"); print(round(ci, 3))
  }

  g <- function(a, v) d[aso == a & vpa == v]$pct_meth
  dd <- (g(1, 1) - g(0, 1)) - (g(1, 0) - g(0, 0))
  cat("\nDifference of differences (percentage points): ", round(dd, 2), "\n", sep = "")
  cat("  ASO effect on VPA background: ", round(g(1, 1) - g(0, 1), 2), " pp\n", sep = "")
  cat("  ASO effect on control background: ", round(g(1, 0) - g(0, 0), 2), " pp\n", sep = "")
  cat("  VPA effect on ASO background: ", round(g(1, 1) - g(1, 0), 2), " pp\n", sep = "")
  cat("  VPA effect on control background: ", round(g(0, 1) - g(0, 0), 2), " pp\n", sep = "")

  if (requireNamespace("aod", quietly = TRUE)) {
    cat("\n-- Beta-binomial GLM (aod::betabin) --\n")
    d2 <- copy(d)[, n := meth + unmeth]
    bb <- try(aod::betabin(cbind(meth, unmeth) ~ aso * vpa, ~ 1, data = d2),
              silent = TRUE)
    if (!inherits(bb, "try-error")) print(bb) else
      cat("betabin did not converge with four groups; quasibinomial stands.\n")
  } else {
    cat("\n(aod not installed; skipping beta-binomial. ",
        "install.packages('aod') to add it.)\n", sep = "")
  }

  cat("\n-- Fisher exact, ASO_VPA versus Scramble_VPA read counts --\n")
  a <- d[aso == 1 & vpa == 1]; b <- d[aso == 0 & vpa == 1]
  if (nrow(a) && nrow(b)) {
    tb <- matrix(c(a$meth, a$unmeth, b$meth, b$unmeth), nrow = 2,
                 dimnames = list(c("meth", "unmeth"), c("ASO_VPA", "Scramble_VPA")))
    print(tb); print(fisher.test(tb))
  }

  invisible(list(binomial = m_bin, quasi = m_q, dispersion = disp, dd = dd))
}

`%||%` <- function(a, b) if (is.null(a)) b else a
fits <- list()
fits[[PRIMARY$name]] <- fit_region(counts, PRIMARY$name)
for (rn in REGIONS[primary == FALSE]$name) {
  if (rn %in% counts$region) fits[[rn]] <- fit_region(counts, rn)
}

rep_counts <- per_cyt[, .(meth = sum(nM), unmeth = sum(nU)),
                      by = .(region, condition, file)]
rep_counts <- merge(rep_counts, COND[, .(condition, aso, vpa)], by = "condition")
n_per_cond <- rep_counts[region == PRIMARY$name, .N, by = condition]
cat("\n== Replicate-level check, primary region ==\n")
print(n_per_cond)

if (nrow(n_per_cond) && max(n_per_cond$N) > 1) {
  d <- rep_counts[region == PRIMARY$name]
  d[, prop := meth / (meth + unmeth)]
  cat("\nPer-replicate methylation proportion:\n")
  print(d[, .(condition, file, meth, unmeth, prop = round(prop, 3))])
  m <- glm(cbind(meth, unmeth) ~ aso * vpa, data = d, family = quasibinomial())
  cat("\nQuasibinomial on replicates:\n")
  print(summary(m)$coefficients)
  cat("\nThis is the test to quote if replicates are available, because it\n")
  cat("treats the replicate rather than the read as the independent unit.\n")
} else {
  cat("Counts are pooled, one value per condition. With four observations and\n")
  cat("three parameters plus an interaction the model has zero residual degrees\n")
  cat("of freedom, so no p-value from the pooled fit is interpretable. Report\n")
  cat("the raw counts and the difference of differences as descriptive, and say\n")
  cat("that replicate-level testing was not possible at this locus depth.\n")
}

pf <- counts[region == PRIMARY$name]
if (nrow(pf)) {
  ci <- t(mapply(function(m, n) {
    if (n == 0) return(c(NA, NA))
    unlist(binom.test(m, n)$conf.int) * 100
  }, pf$meth, pf$meth + pf$unmeth))
  pf[, `:=`(lo = ci[, 1], hi = ci[, 2])]
  pf[, background := ifelse(vpa == 1, "VPA background", "Control background")]
  pf[, oligo := ifelse(aso == 1, "ASO1", "Scramble")]

  pdf(file.path(OUT_DIR, "Fig_smn2_adjacent_2x2.pdf"), width = 6.5, height = 4.5)
  print(
    ggplot(pf, aes(x = oligo, y = pct_meth, fill = background)) +
      geom_col(position = position_dodge(0.7), width = 0.6) +
      geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(0.7),
                    width = 0.15, colour = COL_NEUT) +
      geom_text(aes(label = sprintf("%.1f%%\n%d/%d", pct_meth, meth, meth + unmeth)),
                position = position_dodge(0.7), vjust = -0.4, size = 2.7) +
      scale_fill_manual(values = c("Control background" = COL_NEUT,
                                   "VPA background" = COL_HYPO)) +
      scale_y_continuous(limits = c(0, 118), breaks = seq(0, 100, 25)) +
      labs(x = NULL, y = "Methylation (%)",
           title = sprintf("%s, %s:%s-%s",
                           PRIMARY$name, PRIMARY$chr,
                           format(PRIMARY$start, big.mark = ","),
                           format(PRIMARY$end, big.mark = ",")),
           subtitle = paste("Bars show pooled read proportions with binomial",
                            "95% intervals. Labels give methylated over total reads.")) +
      theme_thesis()
  )
  dev.off()
  cat("\nWrote Fig_smn2_adjacent_2x2.pdf\n")
}

cat("\nOutputs in", OUT_DIR, "\n")
log_session(TAG)
