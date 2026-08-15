.libPaths("~/R/library")
OUT_DIR <- "results/dmr_benchmark"

combine_pattern <- function(pattern, label) {
  files <- list.files(OUT_DIR, pattern=pattern, full.names=TRUE)
  if (length(files) == 0) {
    message("No files found for: ", pattern); return(NULL)
  }
  df <- do.call(rbind, lapply(files, read.csv))
  message(label, ": ", nrow(df), " rows from ", length(files), " files")
  df
}

ls_df <- combine_pattern("parameter_benchmark_labelswap.*\\.csv", "Label swap")
rc_df <- combine_pattern("parameter_benchmark_readcount.*\\.csv", "Read count")
mg_df <- combine_pattern("parameter_benchmark_mingap.*\\.csv",    "MinGap")

if (!is.null(ls_df)) write.csv(ls_df, file.path(OUT_DIR, "parameter_benchmark_label_swap.csv"),    row.names=FALSE)
if (!is.null(rc_df)) write.csv(rc_df, file.path(OUT_DIR, "parameter_benchmark_archie_scramble.csv"), row.names=FALSE)
if (!is.null(mg_df)) write.csv(mg_df, file.path(OUT_DIR, "benchmark_nb_mingap_v2.csv"),             row.names=FALSE)

message("Done. Combined CSVs saved to: ", OUT_DIR)
