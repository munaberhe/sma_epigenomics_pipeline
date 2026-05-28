#!/usr/bin/env Rscript

read_mqc_full <- function(path, stage_label) {
  if (!file.exists(path)) stop("File not found: ", path)
  x <- read.table(path, header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, check.names = FALSE)
  x$stage <- stage_label
  x
}

pre  <- read_mqc_full("results/qc/multiqc_data/multiqc_fastqc.txt",
                      "pre_trim")
post <- read_mqc_full("results/qc/trimmed/multiqc_data/multiqc_fastqc.txt",
                      "post_trim")

combined <- rbind(pre, post)

# Optional: order nicely by Sample then stage
if ("Sample" %in% names(combined)) {
  combined <- combined[order(combined$Sample, combined$stage), ]
}

out <- "results/qc/pre_post_trimming_table_full.csv"
write.table(combined, file = out, sep = ",", row.names = FALSE)
message("Done: ", out)
print(head(combined, 5))
