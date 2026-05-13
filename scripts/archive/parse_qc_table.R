#!/usr/bin/env Rscript

read_mqc <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  x <- read.table(path, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, check.names = FALSE)

  needed <- c("Sample", "Total Sequences", "avg_sequence_length",
              "%GC", "total_deduplicated_percentage")
  missing <- setdiff(needed, names(x))
  if (length(missing) > 0) {
    stop("Missing columns in ", path, ": ",
         paste(missing, collapse = ", "),
         "\nAvailable: ", paste(names(x), collapse = ", "))
  }

  data.frame(
    sample    = x[["Sample"]],
    raw_reads = x[["Total Sequences"]],
    read_len  = x[["avg_sequence_length"]],
    pct_gc    = x[["%GC"]],
    pct_dup   = x[["total_deduplicated_percentage"]],
    stringsAsFactors = FALSE
  )
}

pre  <- read_mqc("results/qc/multiqc_data/multiqc_fastqc.txt")
pre$stage <- "pre_trim"

post <- read_mqc("results/qc/trimmed/multiqc_data/multiqc_fastqc.txt")
post$stage <- "post_trim"

combined <- rbind(pre, post)
combined <- combined[order(combined$sample, combined$stage), ]

out <- "results/qc/pre_post_trimming_table.csv"
write.table(combined, file = out, sep = ",", row.names = FALSE)
message("Done: ", out)
print(head(combined, 12))
