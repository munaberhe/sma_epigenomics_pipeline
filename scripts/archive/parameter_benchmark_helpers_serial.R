# parameter_benchmark_helpers_serial.R
# PLAN B: serial real+scrambled execution with all 16 cores in ONE SnowParam.
# SMA Epigenomics Pipeline — Muna Berhe, QMUL

library(DMRcaller)
library(BiocParallel)
library(parallel)

COV_DIR    <- "results/alignments/bs/by_chr"
OUT_DIR    <- "results/dmr_benchmark_serial"
CHROM      <- "chr1"
REGION_END <- 248956422
RDS_PATH   <- "results/dmr_benchmark/chr1_data.rds"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

INNER_WORKERS <- 16
BPPARAM_FAST  <- SnowParam(workers = INNER_WORKERS, type = "SOCK")

WINDOW_SIZES <- c(100, 200, 300, 500, 1000, 2000)
METHODS      <- c("bins", "neighbourhood", "noise_filter")
MODES        <- c(FALSE, TRUE)

load_chr1_data <- function() {
  if (file.exists(RDS_PATH)) {
    message("Loading chr1 data from cache: ", RDS_PATH)
    dat <- readRDS(RDS_PATH)
    message("  ASO_VPA CpGs:  ", length(dat$aso_vpa))
    message("  ASO_CTRL CpGs: ", length(dat$aso_ctrl))
    return(dat)
  }
  message("No cache at ", RDS_PATH, " — falling back to readBismarkPool().")
  load_one <- function(samples) {
    paths <- file.path(COV_DIR, paste0(samples, "_", CHROM, ".CpG_report.txt.gz"))
    paths <- paths[file.exists(paths)]
    readBismarkPool(paths)
  }
  list(aso_vpa  = load_one(c("ASO_VPA_1", "ASO_VPA_2", "ASO_VPA_3")),
       aso_ctrl = load_one(c("ASO_CTRL_1", "ASO_CTRL_2", "ASO_CTRL_3")))
}

get_thresholds <- function(method, strict) {
  if (!strict) {
    list(pval=0.05, minCpG=1, minDiff=0.1, minSize=1,
         minGap=if (method == "noise_filter") 0 else 200, test="fisher")
  } else {
    list(pval=0.01, minCpG=4,
         minDiff=if (method == "noise_filter") 0.4 else 0.2,
         minSize=if (method == "neighbourhood") 1 else 50,
         minGap=if (method == "noise_filter") 0 else 200,
         test=if (method == "noise_filter") "score" else "fisher")
  }
}

run_dmrs_one <- function(treat, ctrl, method, ws, kernel="triangular",
                         strict=FALSE, region) {
  th <- get_thresholds(method, strict)
  tryCatch({
    if (method == "noise_filter") {
      computeDMRs(treat, ctrl, regions=region, context="CG",
        method="noise_filter", windowSize=ws, kernelFunction=kernel,
        test=th$test, pValueThreshold=th$pval, minCytosinesCount=th$minCpG,
        minProportionDifference=th$minDiff, minGap=th$minGap, minSize=th$minSize,
        minReadsPerCytosine=4, BPPARAM=BPPARAM_FAST, parallel=TRUE)
    } else if (method == "bins") {
      computeDMRs(treat, ctrl, regions=region, context="CG",
        method="bins", binSize=ws, test=th$test, pValueThreshold=th$pval,
        minCytosinesCount=th$minCpG, minProportionDifference=th$minDiff,
        minGap=th$minGap, minSize=th$minSize, minReadsPerCytosine=4,
        BPPARAM=BPPARAM_FAST, parallel=TRUE)
    } else {
      computeDMRs(treat, ctrl, regions=region, context="CG",
        method="neighbourhood", test=th$test, pValueThreshold=th$pval,
        minCytosinesCount=th$minCpG, minProportionDifference=th$minDiff,
        minGap=th$minGap, minSize=th$minSize, minReadsPerCytosine=4,
        BPPARAM=BPPARAM_FAST, parallel=TRUE)
    }
  }, error = function(e) {
    message("  Error: ", e$message)
    structure(list(), class="dmr_error", message=e$message)
  })
}

run_dmrs_pair_parallel <- function(treat_real, ctrl_real, treat_scr, ctrl_scr,
                                   method, ws, kernel, strict, region) {
  t_start <- Sys.time()
  message("  [serial] Real branch...")
  s <- Sys.time()
  out_real <- run_dmrs_one(treat_real, ctrl_real, method, ws, kernel, strict, region)
  t_real <- as.numeric(difftime(Sys.time(), s, units="secs"))
  message("  [serial] Real branch done in ", round(t_real), "s")
  message("  [serial] Scrambled branch...")
  s <- Sys.time()
  out_scr <- run_dmrs_one(treat_scr, ctrl_scr, method, ws, kernel, strict, region)
  t_scr <- as.numeric(difftime(Sys.time(), s, units="secs"))
  message("  [serial] Scrambled branch done in ", round(t_scr), "s")
  list(real=out_real, scr=out_scr, t_real=t_real, t_scr=t_scr,
       t_wall=as.numeric(difftime(Sys.time(), t_start, units="secs")))
}

is_error <- function(x) inherits(x, "dmr_error")

dmr_size_summary <- function(dmrs) {
  if (is_error(dmrs) || length(dmrs) == 0)
    return(list(median=NA_real_, mean=NA_real_, q1=NA_real_, q3=NA_real_,
                min=NA_real_, max=NA_real_))
  w <- as.numeric(width(dmrs))
  list(median=median(w), mean=round(mean(w),1),
       q1=as.numeric(quantile(w,0.25)), q3=as.numeric(quantile(w,0.75)),
       min=min(w), max=max(w))
}

make_result_row <- function(pair, method, ws, mode_lab, kernel, null_name) {
  err_real <- if (is_error(pair$real)) attr(pair$real, "message") else NA_character_
  err_scr  <- if (is_error(pair$scr))  attr(pair$scr,  "message") else NA_character_
  n_real <- if (is_error(pair$real)) NA_integer_ else length(pair$real)
  n_scr  <- if (is_error(pair$scr))  NA_integer_ else length(pair$scr)
  ratio <- if (!is.na(n_real) && !is.na(n_scr) && n_scr > 0) round(n_real/n_scr,2) else
           if (!is.na(n_real) && !is.na(n_scr) && n_scr == 0) Inf else NA_real_
  sz_real <- dmr_size_summary(pair$real)
  sz_scr  <- dmr_size_summary(pair$scr)
  data.frame(method=method, window_size=ws, mode=mode_lab, kernel=kernel,
    scramble_method=null_name, n_real=n_real, n_scrambled=n_scr, ratio=ratio,
    size_real_median=sz_real$median, size_real_mean=sz_real$mean,
    size_real_q1=sz_real$q1, size_real_q3=sz_real$q3,
    size_real_min=sz_real$min, size_real_max=sz_real$max,
    size_scr_median=sz_scr$median, size_scr_mean=sz_scr$mean,
    size_scr_q1=sz_scr$q1, size_scr_q3=sz_scr$q3,
    size_scr_min=sz_scr$min, size_scr_max=sz_scr$max,
    t_real_secs=round(pair$t_real,1), t_scr_secs=round(pair$t_scr,1),
    t_wall_secs=round(pair$t_wall,1), error_real=err_real, error_scr=err_scr,
    stringsAsFactors=FALSE)
}

run_benchmark_loop <- function(aso_vpa, aso_ctrl, aso_vpa_scr, aso_ctrl_scr, null_name) {
  region <- GRanges(CHROM, IRanges(1, REGION_END))
  ckpt_file <- file.path(OUT_DIR, paste0("checkpoint_", null_name, ".rds"))
  if (file.exists(ckpt_file)) {
    message("Resuming from checkpoint: ", ckpt_file)
    results <- readRDS(ckpt_file)
  } else {
    results <- list()
  }
  for (method in METHODS) {
    for (ws in WINDOW_SIZES) {
      if (method == "neighbourhood" && ws != WINDOW_SIZES[1]) next
      ker <- if (method == "noise_filter") "triangular" else "NA"
      for (strict in MODES) {
        mode_lab <- if (strict) "strict" else "loose"
        key <- paste(method, ws, mode_lab, ker, sep="_")
        if (!is.null(results[[key]])) {
          message("\n--- [skip cached] ", key, " ---"); next
        }
        message("\n--- ", mode_lab, " | ", method, " | ws=", ws, " | kernel=", ker, " ---")
        message("  Running real + scrambled SERIALLY (Plan B)...")
        pair <- run_dmrs_pair_parallel(aso_vpa, aso_ctrl, aso_vpa_scr, aso_ctrl_scr,
                                       method, ws, ker, strict, region)
        n_real <- if (is_error(pair$real)) NA_integer_ else length(pair$real)
        n_scr  <- if (is_error(pair$scr))  NA_integer_ else length(pair$scr)
        message("  Real DMRs: ", n_real, "  (", round(pair$t_real), "s)")
        message("  Scrambled: ", n_scr,  "  (", round(pair$t_scr),  "s)")
        message("  Cell wallclock (serial): ", round(pair$t_wall), "s")
        results[[key]] <- make_result_row(pair, method, ws, mode_lab, ker, null_name)
        saveRDS(results, ckpt_file)
      }
    }
  }
  summary_df <- do.call(rbind, results)
  write.csv(summary_df,
            file.path(OUT_DIR, paste0("parameter_benchmark_", null_name, ".csv")),
            row.names=FALSE)
  writeLines(capture.output(sessionInfo()),
             file.path(OUT_DIR, paste0("sessionInfo_", null_name, ".txt")))
  message("\nDone. Results: ", OUT_DIR, "/parameter_benchmark_", null_name, ".csv")
}
