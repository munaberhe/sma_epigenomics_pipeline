#!/usr/bin/env Rscript

parse_bismark_report <- function(path) {
  lines <- readLines(path)

  get_num <- function(pattern) {
    l <- grep(pattern, lines, value = TRUE)[1]
    if (length(l) == 0 || is.na(l)) return(NA_real_)
    as.numeric(gsub("[^0-9]", "", l))
  }

  total    <- get_num("Sequence pairs analysed in total")
  unique   <- get_num("Number of paired-end alignments with a unique best hit")
  no_aln   <- get_num("Sequence pairs with no alignments under any condition")
  not_uniq <- get_num("Sequence pairs did not map uniquely")

  # ASO_CTRL_1_1_val_1_bismark_bt2_PE_report.txt -> ASO_CTRL_1
  base <- basename(path)
  sample <- sub("_\\d+_val_\\d+_bismark_bt2_PE_report.txt$", "", base)

  data.frame(
    Sample          = sample,
    Total_reads     = total,
    Uniquely_mapped = unique,
    Pct_unique      = ifelse(is.na(total) | is.na(unique),
                             NA, round(unique / total * 100, 1)),
    No_alignments   = no_aln,
    Pct_no_align    = ifelse(is.na(total) | is.na(no_aln),
                             NA, round(no_aln / total * 100, 1)),
    Not_unique      = not_uniq,
    Pct_not_unique  = ifelse(is.na(total) | is.na(not_uniq),
                             NA, round(not_uniq / total * 100, 1)),
    stringsAsFactors = FALSE
  )
}

rep_dir <- "results/alignments/bs"
reports <- list.files(rep_dir, pattern = "_bismark_bt2_PE_report.txt$", full.names = TRUE)

if (length(reports) == 0) {
  stop("No *_bismark_bt2_PE_report.txt files found in ", rep_dir)
}

rows <- lapply(reports, parse_bismark_report)
df   <- do.call(rbind, rows)
df   <- df[order(df$Sample), ]

out <- "results/qc/alignment_summary_table.csv"
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
write.table(df, file = out, sep = ",", row.names = FALSE)
message("Done: ", out)
print(df)
