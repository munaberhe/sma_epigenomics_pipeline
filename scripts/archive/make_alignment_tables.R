#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
})

qc_dir   <- "results/qc"
align_csv <- file.path(qc_dir, "alignment_summary_table.csv")
trim_csv  <- file.path(qc_dir, "pre_post_trimming_table_full.csv")

if (!file.exists(align_csv)) {
  stop("Missing ", align_csv, " – run parse_alignment_table.R first?")
}
if (!file.exists(trim_csv)) {
  stop("Missing ", trim_csv, " – pre_post_trimming_table_full.csv not found")
}

align <- read.csv(align_csv, stringsAsFactors = FALSE)
trim  <- read.csv(trim_csv, stringsAsFactors = FALSE)

# Expect columns: Sample, Total.Sequences, stage (pre_trim / post_trim)
if (!all(c("Sample", "Total.Sequences", "stage") %in% names(trim))) {
  stop("pre_post_trimming_table_full.csv must contain Sample, Total.Sequences, stage columns")
}

# Collapse pre/post into one row per logical sample (strip lane/val suffix)
# ASO_CTRL_1_1, ASO_CTRL_1_2 -> ASO_CTRL_1
sample_core <- sub("_(\\d+)(_val_\\d+)?$", "", trim$Sample)

trim_long <- trim %>%
  mutate(Sample_core = sample_core) %>%
  select(Sample_core, stage, Total.Sequences)

trim_wide <- trim_long %>%
  tidyr::pivot_wider(
    id_cols = Sample_core,
    names_from = stage,
    values_from = Total.Sequences,
    values_fn = list(Total.Sequences = sum)
  )

names(trim_wide) <- sub("^pre_trim$", "Raw_reads", names(trim_wide))
names(trim_wide) <- sub("^post_trim$", "Trimmed_reads", names(trim_wide))

# Now align$Sample is already like ASO_CTRL_1, ASO_VPA_1, etc.
trim_sub <- trim_wide
names(trim_sub)[names(trim_sub) == "Sample_core"] <- "Sample"

# Merge on Sample
merged <- merge(trim_sub, align, by = "Sample", all.x = TRUE, sort = TRUE)

# Derived percentages
merged$Pct_retained <- with(merged,
                            ifelse(is.na(Raw_reads) | Raw_reads == 0,
                                   NA, round(Trimmed_reads / Raw_reads * 100, 1)))

# Reorder / select key columns for talk
out_df <- merged[, c("Sample",
                     "Raw_reads",
                     "Trimmed_reads",
                     "Pct_retained",
                     "Total_reads",
                     "Uniquely_mapped",
                     "Pct_unique",
                     "No_alignments",
                     "Pct_no_align")]

out_file <- file.path(qc_dir, "pre_post_alignment_table.csv")
write.csv(out_df, out_file, row.names = FALSE)
message("Wrote ", out_file)
