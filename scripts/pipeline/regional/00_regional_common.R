## 00_regional_common.R
## Shared configuration and helpers for the regional analysis suite.
## Sourced by 40_, 41_ and 42_. Edit the CONFIG block, nothing else.
.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(ggplot2)
})

PROJ      <- Sys.getenv("SMA_PROJ", "/gpfs/scratch/bt25018/sma_epigenomics_pipeline")
RESULTS   <- file.path(PROJ, "results")
DMR_DIR   <- file.path(RESULTS, "dmr")
METH_DIR  <- file.path(RESULTS, "alignments", "bs", "by_chr")
CX_DIR    <- file.path(RESULTS, "smn2_predup_cx")
OUT_DIR   <- file.path(RESULTS, "regional_suite")
CACHE_DIR <- file.path(OUT_DIR, "cache")
dir.create(OUT_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

SEED     <- 42
N_PERM   <- 1000
BIN_SIZE <- 1e6
SUBSAMPLE <- 50000

CONTRASTS <- c(
  "ASO alone"     = "ASO_CTRL_vs_Scramble_CTRL",
  "VPA alone"     = "Scramble_VPA_vs_Scramble_CTRL",
  "ASO in VPA"    = "ASO_VPA_vs_Scramble_VPA",
  "VPA in ASO"    = "ASO_VPA_vs_ASO_CTRL"
)

COL_HYPER <- "#A84B2F"
COL_HYPO  <- "#D19900"
COL_NEUT  <- "#1B474D"
COL_GREY  <- "#9A9A9A"

HG38 <- data.table(
  chr = c(paste0("chr", 1:22), "chrX", "chrY"),
  len = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
          159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
          114364328, 107043718, 101991189,  90338345,  83257441,  80373285,
           58617616,  64444167,  46709983,  50818468, 156040895,  57227415)
)
MAIN_CHR <- HG38$chr[HG38$chr != "chrY"]

theme_thesis <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = base_size + 1),
      plot.subtitle = element_text(size = base_size - 1, colour = "grey30"),
      legend.position = "top",
      legend.title = element_blank()
    )
}

discover_dmr_files <- function(dir = DMR_DIR, verbose = TRUE) {
  f <- list.files(dir, pattern = "^dmr_.*\\.rds$",
                  full.names = TRUE, recursive = FALSE)
  f <- f[!grepl("tested_windows|by_chr|dmr_all_contrasts", f)]
  if (verbose) {
    cat("Files under", dir, ":\n")
    if (!length(f)) cat("  none found\n") else
      cat(paste0("  ", basename(f), "  (", round(file.size(f) / 1e6, 1), " MB)\n"),
          sep = "")
  }
  invisible(f)
}

resolve_contrast_files <- function(files, contrasts = CONTRASTS) {
  out <- character(0)
  for (lab in names(contrasts)) {
    pat <- contrasts[[lab]]
    hit <- files[grepl(pat, basename(files), ignore.case = TRUE)]
    if (length(hit) == 0)
      stop("No DMR file matched pattern '", pat, "' for contrast '", lab,
           "'. Run discover_dmr_files() and edit CONTRASTS.")
    if (length(hit) > 1)
      stop("Pattern '", pat, "' matched ", length(hit), " files: ",
           paste(basename(hit), collapse = ", "),
           ". Make CONTRASTS more specific.")
    out[lab] <- hit
  }
  out
}

read_dmr <- function(path, verbose = TRUE) {
  if (grepl("\\.rds$", path)) {
    gr <- readRDS(path)
    dt <- as.data.table(as.data.frame(gr))
    setnames(dt, "seqnames", "chr", skip_absent=TRUE)
    dt[, chr := as.character(chr)]
    if ("proportion1" %in% names(dt) && "proportion2" %in% names(dt)) {
      dt[, diff := proportion1 - proportion2]
      dt[, direction := ifelse(diff > 0, "hyper", "hypo")]
    } else if ("regionType" %in% names(dt)) {
      dt[, direction := ifelse(regionType == "gain", "hypo", "hyper")]
    } else {
      dt[, direction := NA_character_]
    }
    dt <- dt[chr %in% MAIN_CHR]
    if (verbose) cat("  ", basename(path), ": ", format(nrow(dt), big.mark=","),
                     " regions (RDS)\n", sep="")
    return(dt[])
  }

  dt <- if (grepl("\\.bed(\\.gz)?$", path)) {
    fread(path, header = FALSE, select = 1:3,
          col.names = c("chr", "start", "end"))
  } else {
    fread(path)
  }
  nm <- tolower(names(dt))
  pick <- function(cands) {
    i <- which(nm %in% cands)
    if (length(i)) names(dt)[i[1]] else NA_character_
  }
  c_chr   <- pick(c("chr", "seqnames", "chrom", "chromosome", "v1"))
  c_start <- pick(c("start", "chromstart", "v2"))
  c_end   <- pick(c("end", "chromend", "stop", "v3"))
  c_diff  <- pick(c("meth_diff", "methdiff", "diff", "proportiondifference",
                    "delta", "delta_meth", "sumreadsm1", "meandiff"))
  if (any(is.na(c(c_chr, c_start, c_end))))
    stop("Could not identify chr/start/end in ", basename(path),
         ". Columns present: ", paste(names(dt), collapse = ", "))
  keep <- c(chr = c_chr, start = c_start, end = c_end)
  if (!is.na(c_diff)) keep <- c(keep, diff = c_diff)
  dt <- dt[, ..keep]
  setnames(dt, names(keep))
  dt[, chr := as.character(chr)]
  if (!any(grepl("^chr", dt$chr))) dt[, chr := paste0("chr", chr)]
  dt <- dt[chr %in% MAIN_CHR]
  dt[, start := as.integer(start)][, end := as.integer(end)]
  if ("diff" %in% names(dt)) {
    dt[, direction := ifelse(diff > 0, "hyper", "hypo")]
  } else {
    dt[, direction := NA_character_]
  }
  if (verbose) {
    cat("  ", basename(path), ": ", format(nrow(dt), big.mark = ","),
        " regions on main chromosomes", sep = "")
    if (!is.na(c_diff)) {
      cat("; effect column '", c_diff, "'; ",
          format(sum(dt$direction == "hyper"), big.mark = ","), " hyper / ",
          format(sum(dt$direction == "hypo"),  big.mark = ","), " hypo", sep = "")
    }
    cat("\n")
  }
  dt[]
}

as_gr <- function(dt) {
  GRanges(dt$chr, IRanges(dt$start, dt$end))
}

make_bins <- function(bin = BIN_SIZE, chrs = MAIN_CHR) {
  h <- HG38[chr %in% chrs]
  rbindlist(lapply(seq_len(nrow(h)), function(i) {
    starts <- seq(1, h$len[i], by = bin)
    data.table(chr = h$chr[i], start = starts,
               end = pmin(starts + bin - 1L, h$len[i]))
  }))[, bin_id := paste0(chr, ":", start)][]
}

log_session <- function(tag) {
  f <- file.path(OUT_DIR, paste0("sessionInfo_", tag, ".txt"))
  writeLines(capture.output(sessionInfo()), f)
  cat("Session info written to", f, "\n")
}

COND <- data.table(
  condition = c("Scramble_CTRL", "ASO_CTRL", "Scramble_VPA", "ASO_VPA"),
  pattern   = c("Scramble_CTRL|Scr_CTRL", "ASO_CTRL", "Scramble_VPA|Scr_VPA", "ASO_VPA"),
  aso       = c(0L, 1L, 0L, 1L),
  vpa       = c(0L, 0L, 1L, 1L)
)
