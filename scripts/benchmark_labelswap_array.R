.libPaths("~/R/library")
source("scripts/parameter_benchmark_helpers.R")

TASK_ID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

methods      <- c("bins", "neighbourhood", "noise_filter")
window_sizes <- c(100, 200, 300, 500, 1000, 2000)
modes        <- c(FALSE, TRUE)

combos <- expand.grid(method=methods, window_size=window_sizes,
                      strict=modes, stringsAsFactors=FALSE)
combos <- combos[!(combos$method=="neighbourhood" &
                   combos$window_size != window_sizes[1]), ]
combos <- combos[!(combos$method=="noise_filter" & combos$strict==FALSE), ]

if (TASK_ID > nrow(combos)) {
  message("TASK_ID ", TASK_ID, " exceeds combo table (nrow=", nrow(combos), ") -- skipping")
  quit(status=0)
}
method <- combos$method[TASK_ID]
ws     <- combos$window_size[TASK_ID]
strict <- combos$strict[TASK_ID]
mode_lab <- if(strict) "strict" else "loose"
message("Task ", TASK_ID, ": ", method, " ws=", ws, " ", mode_lab)

dat <- load_chr1_data()
aso_vpa      <- dat$aso_vpa
aso_ctrl     <- dat$aso_ctrl
aso_vpa_scr  <- aso_ctrl
aso_ctrl_scr <- aso_vpa

region <- GRanges(CHROM, IRanges(1, REGION_END))
ker <- if (method == "noise_filter") "triangular" else "NA"

pair <- run_dmrs_pair_parallel(aso_vpa, aso_ctrl, aso_vpa_scr, aso_ctrl_scr,
                               method, ws, ker, strict, region)
row <- make_result_row(pair, method, ws, mode_lab, ker, "label_swap")

out_csv <- file.path(OUT_DIR,
  sprintf("parameter_benchmark_labelswap_%s_ws%d_%s.csv", method, ws, mode_lab))
write.csv(row, out_csv, row.names=FALSE)
message("Saved: ", out_csv)
