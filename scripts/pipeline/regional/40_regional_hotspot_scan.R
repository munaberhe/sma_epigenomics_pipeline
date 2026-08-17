## 40_regional_hotspot_scan.R
## Regional (10 Mb bin) and chromosome-level DMR density scan for all four
## pairwise contrasts, with a covered-CpG background.
##
## Rationale. Per-locus label permutation on n=3 per condition has a hard
## resolution floor, because 3 versus 3 admits only a small number of distinct
## label arrangements. Regional tests avoid that floor: the unit of inference is
## a bin containing many CpGs and many called regions, and the null is built by
## randomising region positions within the covered genome rather than by
## shuffling sample labels. This is the same class of null used for the enhancer
## enrichment, which is why that analysis reached the permutation floor while
## the single-locus tests did not.
##
## Status. No chromosome- or region-level hotspot search has been run under the
## pairwise framework. This script is that search, run de novo on the four
## current contrasts. Any earlier hotspot observation, including the chr13 one,
## came from a contrast definition that has since been retired and carries no
## standing here. chr13 is reported below as a named legacy locus being tested,
## not as a result being confirmed. If it does not reappear, that is the finding.
##
## Critical control. Raw DMR density per bin tracks sequencing coverage. Every
## expectation below is weighted by the number of CpGs actually covered in that
## bin, so a bin is only called dense if it exceeds what its own coverage
## predicts. Without this no density signal can be distinguished from a
## coverage artefact.
##
## Usage on Apocrita:
##   Rscript 40_regional_hotspot_scan.R            # full run
##   DISCOVER_ONLY=TRUE Rscript 40_regional_hotspot_scan.R   # path check only

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
DISCOVER_ONLY <- toupper(Sys.getenv("DISCOVER_ONLY", "FALSE")) == "TRUE"
TAG <- "40_regional"

cat("== 40_regional_hotspot_scan ==\n")
files <- discover_dmr_files()

if (DISCOVER_ONLY) {
  cat("\nDISCOVER_ONLY set. Also listing methylation directory:\n")
  print(head(list.files(METH_DIR, full.names = FALSE), 40))
  quit(save = "no")
}

dmr_files <- resolve_contrast_files(files)
cat("\nResolved contrast files:\n"); print(dmr_files)

cat("\nReading DMR tables:\n")
dmr <- lapply(dmr_files, read_dmr)
names(dmr) <- names(dmr_files)

bins <- make_bins()
cat("\n", nrow(bins), " bins of ", BIN_SIZE / 1e6, " Mb across ",
    length(MAIN_CHR), " chromosomes\n", sep = "")

METH_PATTERN <- Sys.getenv("METH_PATTERN", "")
MIN_COV      <- as.integer(Sys.getenv("MIN_COV", "4"))
cache_bg <- file.path(CACHE_DIR, sprintf("binned_covered_cpg_%dMb_cov%d.rds",
                                         BIN_SIZE / 1e6, MIN_COV))

bin_background <- function() {
  if (file.exists(cache_bg)) {
    cat("Using cached coverage background:", basename(cache_bg), "\n")
    return(readRDS(cache_bg))
  }
  mf <- list.files(METH_DIR, full.names = TRUE, recursive = TRUE,
                   pattern = if (nzchar(METH_PATTERN)) METH_PATTERN else
                     "\\.(txt|tsv|cov|CX_report|bedGraph)(\\.gz)?$")
  if (!length(mf)) {
    warning("No methylation files found under ", METH_DIR,
            ". Falling back to a UNIFORM background. ",
            "Results are NOT coverage-controlled and must not be reported ",
            "as enrichment. Set METH_PATTERN and rerun.")
    bg <- copy(bins)[, covered := (end - start + 1L)]
    return(bg)
  }
  # Use one representative sample (ASO_CTRL_1), but read ALL its per-chromosome
  # files, not just the first one, so every chromosome gets real coverage.
  rep_sample <- "ASO_CTRL_1"
  mf_rep <- mf[grepl(rep_sample, basename(mf))]
  if (!length(mf_rep)) mf_rep <- mf  # fallback to all files if pattern fails
  cat("Building coverage background from ", length(mf_rep),
      " per-chromosome files for representative sample: ", rep_sample, "\n", sep = "")

  agg_list <- lapply(mf_rep, function(use) {
    x <- fread(use, header = FALSE, showProgress = FALSE)
    if (ncol(x) >= 6 && is.numeric(x$V2) && is.numeric(x$V3) &&
        all(x$V3[1:min(100, nrow(x))] >= x$V2[1:min(100, nrow(x))])) {
      setnames(x, 1:6, c("chr", "start", "end", "pct", "nM", "nU"))
      x[, pos := start]
    } else {
      setnames(x, 1:5, c("chr", "pos", "strand", "nM", "nU"))
      if (ncol(x) >= 6) x <- x[V6 == "CG"]
    }
    x[, chr := as.character(chr)]
    if (!any(grepl("^chr", x$chr))) x[, chr := paste0("chr", chr)]
    x <- x[chr %in% MAIN_CHR][(nM + nU) >= MIN_COV]
    x[, bin_start := ((pos - 1L) %/% as.integer(BIN_SIZE)) * as.integer(BIN_SIZE) + 1L]
    x[, .(covered = .N), by = .(chr, start = bin_start)]
  })
  agg <- rbindlist(agg_list)
  agg <- agg[, .(covered = sum(covered)), by = .(chr, start)]
  cat("  Total covered CpGs across all chromosomes: ",
      format(sum(agg$covered), big.mark = ","), "\n", sep = "")
  bg <- merge(bins, agg, by = c("chr", "start"), all.x = TRUE)
  bg[is.na(covered), covered := 0L]
  saveRDS(bg, cache_bg)
  bg
}

bg <- bin_background()
cat("Background: ", format(sum(bg$covered), big.mark = ","),
    " covered CpGs over ", sum(bg$covered > 0), " informative bins\n", sep = "")

scan_contrast <- function(dt, label) {
  cat("\n-- ", label, " --\n", sep = "")
  d <- copy(dt)
  d[, mid := as.integer((start + end) / 2)]
  d[, bin_start := ((mid - 1L) %/% as.integer(BIN_SIZE)) * as.integer(BIN_SIZE) + 1L]
  obs <- d[, .(n_dmr = .N,
               n_hyper = sum(direction == "hyper", na.rm = TRUE),
               n_hypo  = sum(direction == "hypo",  na.rm = TRUE)),
           by = .(chr, start = bin_start)]
  tab <- merge(bg, obs, by = c("chr", "start"), all.x = TRUE)
  for (v in c("n_dmr", "n_hyper", "n_hypo")) tab[is.na(get(v)), (v) := 0L]
  informative <- tab$covered > 0
  p_bin <- tab$covered / sum(tab$covered)
  N <- sum(tab$n_dmr)
  tab[, expected := N * p_bin]
  tab[, obs_exp := ifelse(expected > 0, n_dmr / expected, NA_real_)]
  tab[, log2_oe := log2(obs_exp)]
  cat("  ", format(N, big.mark = ","), " regions placed into ",
      sum(tab$n_dmr > 0), " bins\n", sep = "")
  cat("  running ", N_PERM, " genomic randomisations\n", sep = "")
  ge <- integer(nrow(tab))
  for (i in seq_len(N_PERM)) {
    sim <- rmultinom(1, size = N, prob = p_bin)[, 1]
    ge <- ge + as.integer(sim >= tab$n_dmr)
  }
  tab[, p_emp := (ge + 1) / (N_PERM + 1)]
  tab[!informative, p_emp := NA_real_]
  tab[, p_adj := p.adjust(p_emp, method = "BH")]
  tab[, contrast := label]
  chrom <- tab[, .(n_dmr = sum(n_dmr), n_hyper = sum(n_hyper),
                   n_hypo = sum(n_hypo), covered = sum(covered)),
               by = chr]
  pc <- chrom$covered / sum(chrom$covered)
  chrom[, expected := sum(n_dmr) * pc]
  chrom[, obs_exp := n_dmr / expected]
  chrom[, log2_oe := log2(obs_exp)]
  gec <- integer(nrow(chrom))
  for (i in seq_len(N_PERM)) {
    sim <- rmultinom(1, size = sum(chrom$n_dmr), prob = pc)[, 1]
    gec <- gec + as.integer(sim >= chrom$n_dmr)
  }
  chrom[, p_emp := (gec + 1) / (N_PERM + 1)]
  chrom[, p_adj := p.adjust(p_emp, method = "BH")]
  chrom[, contrast := label]
  list(bins = tab[], chrom = chrom[])
}

res <- lapply(names(dmr), function(l) scan_contrast(dmr[[l]], l))
names(res) <- names(dmr)
bin_all   <- rbindlist(lapply(res, `[[`, "bins"))
chrom_all <- rbindlist(lapply(res, `[[`, "chrom"))
bin_all[, contrast := factor(contrast, levels = names(CONTRASTS))]
chrom_all[, contrast := factor(contrast, levels = names(CONTRASTS))]
chrom_all[, chr := factor(chr, levels = MAIN_CHR)]

bin_suffix <- sprintf("%dMb", as.integer(BIN_SIZE/1e6))
fwrite(bin_all,   file.path(OUT_DIR, paste0("regional_bins_all_contrasts_", bin_suffix, ".csv")))
fwrite(chrom_all, file.path(OUT_DIR, paste0("regional_chromosome_all_contrasts_", bin_suffix, ".csv")))

cat("\n== Bins passing BH < 0.05, ranked by log2 observed/expected ==\n")
top <- bin_all[!is.na(p_adj) & p_adj < 0.05 & obs_exp > 1][order(-log2_oe)]
if (nrow(top)) {
  print(head(top[, .(contrast, chr, start, end, n_dmr,
                     expected = round(expected, 1),
                     log2_oe = round(log2_oe, 2),
                     p_emp, p_adj = signif(p_adj, 3))], 30))
} else {
  cat("None. Report this as a null result: DMR density is consistent with\n")
  cat("covered-CpG density genome-wide at 10 Mb resolution.\n")
}
fwrite(top, file.path(OUT_DIR, paste0("regional_bins_significant_", bin_suffix, ".csv")))

LEGACY_CHR <- c("chr13", "chrX")
cat("\n== Legacy loci under the pairwise framework, chromosome level ==\n")
cat("Neither locus was defined by the four current contrasts. Treat a null\n")
cat("result here as informative, not as a failure.\n")
print(chrom_all[chr %in% LEGACY_CHR,
                .(contrast, chr, n_dmr, expected = round(expected, 1),
                  log2_oe = round(log2_oe, 2), p_emp, p_adj = signif(p_adj, 3))])

cat("\n== chr13 bins, all contrasts ==\n")
print(bin_all[chr == "chr13",
              .(contrast, start, n_dmr, expected = round(expected, 1),
                log2_oe = round(log2_oe, 2), p_emp,
                p_adj = signif(p_adj, 3))][order(contrast, start)])

cat("\n== chr13 replication verdict ==\n")
for (l in names(CONTRASTS)) {
  sub <- bin_all[chr == "chr13" & contrast == l & covered > 0]
  if (!nrow(sub)) { cat("  ", l, ": no informative chr13 bins\n", sep = ""); next }
  hit <- sub[!is.na(p_adj) & p_adj < 0.05 & obs_exp > 1]
  cat("  ", l, ": ", nrow(hit), " of ", nrow(sub),
      " chr13 bins enriched at BH < 0.05", sep = "")
  if (nrow(hit)) {
    cat(" (max log2 O/E ", round(max(hit$log2_oe), 2), " at ",
        format(hit[which.max(log2_oe)]$start, big.mark = ","), ")", sep = "")
  } else {
    cat(" (max log2 O/E ", round(max(sub$log2_oe, na.rm = TRUE), 2), ")", sep = "")
  }
  cat("\n")
}
cat("If every contrast returns zero enriched chr13 bins, the reportable\n")
cat("statement is that the legacy chr13 hotspot does not replicate under the\n")
cat("pairwise framework once coverage is controlled.\n")

pdf(file.path(OUT_DIR, "Fig_regional_bin_density.pdf"), width = 11, height = 7)
bin_plot <- bin_all[covered > 0]
bin_plot[, chr := factor(chr, levels = MAIN_CHR)]
bin_plot[, sig := !is.na(p_adj) & p_adj < 0.05]
print(
  ggplot(bin_plot, aes(x = start / 1e6, y = log2_oe)) +
    geom_hline(yintercept = 0, colour = COL_NEUT, linewidth = 0.3) +
    geom_point(aes(colour = sig), size = 0.8) +
    scale_colour_manual(values = c(`FALSE` = COL_GREY, `TRUE` = COL_HYPER),
                        labels = c("BH >= 0.05", "BH < 0.05")) +
    facet_grid(contrast ~ chr, scales = "free_x", space = "free_x") +
    labs(x = "Position (Mb)",
         y = "log2 observed / expected DMR count",
         title = "Regional DMR density against a covered-CpG background",
         subtitle = sprintf("%.0f Mb bins, %d genomic randomisations, seed %d",
                            BIN_SIZE / 1e6, N_PERM, SEED)) +
    theme_thesis() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.spacing.x = unit(0.02, "lines"),
          strip.text.x = element_text(size = 6, angle = 90))
)
dev.off()

pdf(file.path(OUT_DIR, "Fig_regional_chromosome_oe.pdf"), width = 9, height = 6)
chrom_all[, sig := !is.na(p_adj) & p_adj < 0.05]
print(
  ggplot(chrom_all, aes(x = chr, y = log2_oe, fill = sig)) +
    geom_col() +
    geom_hline(yintercept = 0, colour = COL_NEUT, linewidth = 0.3) +
    scale_fill_manual(values = c(`FALSE` = COL_GREY, `TRUE` = COL_HYPER),
                      labels = c("BH >= 0.05", "BH < 0.05")) +
    facet_wrap(~ contrast, ncol = 1, scales = "free_y") +
    labs(x = NULL, y = "log2 observed / expected DMR count",
         title = "Chromosome-level DMR density against a covered-CpG background",
         subtitle = sprintf("%d genomic randomisations, seed %d", N_PERM, SEED)) +
    theme_thesis() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 7))
)
dev.off()

cat("\nWrote figures and tables to", OUT_DIR, "\n")
log_session(TAG)
