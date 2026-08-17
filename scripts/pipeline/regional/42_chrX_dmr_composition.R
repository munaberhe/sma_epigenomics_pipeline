## 42_chrX_dmr_composition.R
## Quantifies the chrX representation among DMRs in each pairwise contrast.
##
## Status. The chrX excess in the ASO-alone set was noted during review and never
## quantified. This script tests it properly. It has never been run under the
## pairwise framework, so all four contrasts are treated as first-pass.
##
## Three questions, in order of how defensible the answer will be:
##
##   1. Is chrX over-represented relative to its covered-CpG share? This is the
##      only one with a real null. Coverage weighting is essential: HEK293T is
##      female-derived and X coverage differs from autosomal coverage, so a raw
##      chromosome proportion would be uninterpretable.
##   2. Does the hypo to hyper ratio on chrX differ from the autosomal ratio?
##      X-inactivation maintains one silenced X, so the starting methylation
##      distribution on X is not the autosomal one, and a directional difference
##      may reflect that baseline rather than treatment.
##   3. Do the chrX DMRs fall near the X-inactivation centre or in genes with a
##      documented escape phenotype? This is descriptive positioning, not a test.
##
## Interpretation caution. Any chrX effect in this design is confounded with
## X-inactivation status and with the copy-number and karyotypic instability
## documented in HEK293-lineage cells. The script prints that caveat with the
## result so it cannot be reported without it.
##
## Usage:
##   Rscript 42_chrX_dmr_composition.R
##   DISCOVER_ONLY=TRUE Rscript 42_chrX_dmr_composition.R

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
TAG <- "42_chrX"
DISCOVER_ONLY <- toupper(Sys.getenv("DISCOVER_ONLY", "FALSE")) == "TRUE"

cat("== 42_chrX_dmr_composition ==\n")
files <- discover_dmr_files()

if (DISCOVER_ONLY) quit(save = "no")

dmr_files <- resolve_contrast_files(files)

cat("\nReading DMR tables:\n")
dmr <- lapply(dmr_files, read_dmr)
names(dmr) <- names(dmr_files)

cache_bg <- list.files(CACHE_DIR, pattern = "^binned_covered_cpg_.*\\.rds$",
                       full.names = TRUE)
if (length(cache_bg)) {
  cat("\nUsing coverage cache:", basename(cache_bg[1]), "\n")
  bgb <- readRDS(cache_bg[1])
  chr_bg <- bgb[, .(covered = sum(covered)), by = chr]
} else {
  stop("No coverage cache found in ", CACHE_DIR,
       ". Run 40_regional_hotspot_scan.R first: it builds and caches the ",
       "covered-CpG background that this script needs for its null.")
}
chr_bg[, frac_covered := covered / sum(covered)]
cat("\nchrX share of covered CpGs: ",
    round(100 * chr_bg[chr == "chrX"]$frac_covered, 3), "%\n", sep = "")

test_chrX <- function(dt, label) {
  d <- copy(dt)
  n_tot <- nrow(d)
  n_x   <- sum(d$chr == "chrX")
  p_x   <- chr_bg[chr == "chrX"]$frac_covered
  exp_x <- n_tot * p_x
  bt <- binom.test(n_x, n_tot, p = p_x)
  probs <- chr_bg[match(MAIN_CHR, chr)]$frac_covered
  probs[is.na(probs)] <- 0
  probs <- probs / sum(probs)
  ix <- which(MAIN_CHR == "chrX")
  sims <- rmultinom(N_PERM, size = n_tot, prob = probs)[ix, ]
  p_emp <- (sum(sims >= n_x) + 1) / (N_PERM + 1)
  data.table(
    contrast   = label,
    n_total    = n_tot,
    n_chrX     = n_x,
    pct_chrX   = 100 * n_x / n_tot,
    expected   = exp_x,
    pct_exp    = 100 * p_x,
    obs_exp    = n_x / exp_x,
    log2_oe    = log2(n_x / exp_x),
    p_binom    = bt$p.value,
    ci_lo      = 100 * bt$conf.int[1],
    ci_hi      = 100 * bt$conf.int[2],
    p_perm     = p_emp
  )
}

q1 <- rbindlist(lapply(names(dmr), function(l) test_chrX(dmr[[l]], l)))
q1[, p_binom_adj := p.adjust(p_binom, method = "BH")]
q1[, p_perm_adj  := p.adjust(p_perm,  method = "BH")]
q1[, contrast := factor(contrast, levels = names(CONTRASTS))]
setorder(q1, contrast)

cat("\n== chrX representation against covered-CpG expectation ==\n")
print(q1[, .(contrast, n_total, n_chrX,
             pct_chrX = round(pct_chrX, 3), pct_exp = round(pct_exp, 3),
             obs_exp = round(obs_exp, 2), log2_oe = round(log2_oe, 2),
             p_perm, p_perm_adj = signif(p_perm_adj, 3),
             p_binom = signif(p_binom, 3))])
fwrite(q1, file.path(OUT_DIR, "chrX_representation.csv"))

cat("\nNote on resolution: the permutation p cannot fall below ",
    signif(1 / (N_PERM + 1), 3), " with ", N_PERM, " iterations. A value at\n",
    "that floor means the observed count was never reached in any draw, not\n",
    "that the p-value is exactly that number.\n", sep = "")

dir_tab <- rbindlist(lapply(names(dmr), function(l) {
  d <- dmr[[l]]
  if (all(is.na(d$direction))) {
    cat("\n[", l, "] no effect-size column found; direction test skipped.\n", sep = "")
    return(NULL)
  }
  d[, arm := ifelse(chr == "chrX", "chrX", "autosome")]
  t <- d[, .(hyper = sum(direction == "hyper"),
             hypo  = sum(direction == "hypo")), by = arm]
  t[, pct_hypo := 100 * hypo / (hypo + hyper)]
  t[, contrast := l]
  t[]
}))

if (nrow(dir_tab)) {
  cat("\n== Directional composition, chrX versus autosomes ==\n")
  print(dir_tab[, .(contrast, arm, hyper, hypo, pct_hypo = round(pct_hypo, 1))])
  fwrite(dir_tab, file.path(OUT_DIR, "chrX_direction.csv"))
  cat("\nFisher exact, hypo/hyper on chrX versus autosomes:\n")
  for (l in unique(dir_tab$contrast)) {
    s <- dir_tab[contrast == l]
    if (nrow(s) < 2) next
    m <- as.matrix(s[order(arm), .(hypo, hyper)])
    rownames(m) <- s[order(arm)]$arm
    ft <- fisher.test(m)
    cat("  ", l, ": OR = ", round(unname(ft$estimate), 3),
        ", 95% CI ", round(ft$conf.int[1], 3), " to ", round(ft$conf.int[2], 3),
        ", p = ", signif(ft$p.value, 3), "\n", sep = "")
  }
  cat("\nOR above 1 means chrX regions are more often hypomethylated than\n")
  cat("autosomal regions in that contrast. Because one X is constitutively\n")
  cat("silenced and therefore methylated at baseline, a hypomethylation excess\n")
  cat("on chrX is expected from the starting distribution and is not evidence\n")
  cat("of a treatment-specific X effect.\n")
}

XIC <- GRanges("chrX", IRanges(72000000, 74000000))
XIST <- GRanges("chrX", IRanges(73820651, 73852753))
ESCAPE <- data.table(
  gene  = c("XIST", "KDM6A", "DDX3X", "ZFX", "EIF1AX", "USP9X", "KDM5C",
            "RPS4X", "SMC1A", "PUDP", "STS", "KAL1_ANOS1", "TXLNG", "CA5B"),
  chr   = "chrX",
  start = c(73820651, 44873182, 41333308, 24148933, 20125515, 41082943, 53190852,
            72272611, 53374149, 55010351, 7137498, 8595996, 16804500, 15756830),
  end   = c(73852753, 45112602, 41364472, 24216255, 20148491, 41215399, 53220514,
            72278243, 53424307, 55070595, 7288363, 8815466, 16863172, 15790756)
)

cat("\n== chrX DMR positioning ==\n")
pos_out <- rbindlist(lapply(names(dmr), function(l) {
  dx <- dmr[[l]][chr == "chrX"]
  if (!nrow(dx)) { cat("  ", l, ": no chrX regions\n", sep = ""); return(NULL) }
  gr <- as_gr(dx)
  n_xic  <- sum(overlapsAny(gr, XIC))
  n_xist <- sum(overlapsAny(gr, XIST))
  esc_gr <- GRanges(ESCAPE$chr, IRanges(ESCAPE$start, ESCAPE$end), gene = ESCAPE$gene)
  ov <- findOverlaps(gr, esc_gr)
  genes <- if (length(ov)) sort(unique(esc_gr$gene[subjectHits(ov)])) else character(0)
  cat("  ", l, ": ", nrow(dx), " chrX regions; ", n_xic,
      " within the Xq13 XIC interval; ", n_xist, " overlapping XIST",
      sep = "")
  if (length(genes)) cat("; escape-gene overlaps: ", paste(genes, collapse = ", "))
  cat("\n")
  data.table(contrast = l, n_chrX = nrow(dx), n_in_XIC = n_xic,
             n_in_XIST = n_xist,
             escape_genes = paste(genes, collapse = ";"))
}))
if (nrow(pos_out)) fwrite(pos_out, file.path(OUT_DIR, "chrX_positioning.csv"))

pdf(file.path(OUT_DIR, "Fig_chrX_representation.pdf"), width = 7, height = 4.5)
print(
  ggplot(q1, aes(x = contrast, y = log2_oe,
                 fill = p_perm_adj < 0.05)) +
    geom_col(width = 0.6) +
    geom_hline(yintercept = 0, colour = COL_NEUT, linewidth = 0.4) +
    geom_text(aes(label = sprintf("%s / %.0f", format(n_chrX, big.mark = ","), expected)),
              vjust = ifelse(q1$log2_oe >= 0, -0.4, 1.3), size = 2.8) +
    scale_fill_manual(values = c(`FALSE` = COL_GREY, `TRUE` = COL_HYPER),
                      labels = c("BH >= 0.05", "BH < 0.05")) +
    labs(x = NULL, y = "log2 observed / expected chrX regions",
         title = "chrX representation among DMRs, against a covered-CpG null",
         subtitle = sprintf("Labels give observed over expected counts. %d randomisations, seed %d.",
                            N_PERM, SEED)) +
    theme_thesis()
)
dev.off()

cat("\n-----------------------------------------------------------\n")
cat("Mandatory caveat for any chrX statement in the thesis:\n")
cat("HEK293-lineage cells are karyotypically unstable with documented\n")
cat("copy-number variation, and X-inactivation gives chrX a baseline\n")
cat("methylation distribution unlike the autosomes. An enrichment here is\n")
cat("therefore consistent with a chromosome-level baseline difference and\n")
cat("cannot be attributed to treatment without an orthogonal measurement.\n")
cat("-----------------------------------------------------------\n")

cat("\nOutputs in", OUT_DIR, "\n")
log_session(TAG)
