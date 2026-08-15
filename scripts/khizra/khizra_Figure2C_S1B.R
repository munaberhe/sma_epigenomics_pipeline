library(GenomicRanges)
library(rtracklayer)
library(UpSetR)
library(ggplot2)
library(gridExtra)
library(grid)
library(ggrepel)
library(reshape2)

setwd('/Users/hmz453/Documents/PostDoc/Research_Articles/Pancretic_Islets/Codes/H9/Analysis/Final_Plots/')


load("/Users/hmz453/Documents/PostDoc/Research_Articles/Pancretic_Islets/Codes/H9/CNN/H9_CNN_pancretic_islet_non_promoter_enhancers_0.51(200bp).Rda")
CNN <- non_promoter_enhancers

load("/Users/hmz453/Documents/PostDoc/Research_Articles/Pancretic_Islets/Codes/H9/LR/H9_LR_pancretic_islet_non_promoter_enhancers_0.04(200bp).Rda")
LR <- non_promoter_enhancers

load("/Users/hmz453/Documents/PostDoc/Research_Articles/Pancretic_Islets/Codes/H9/XGBoost/H9_XGBoost_pancretic_islet_non_promoter_enhancers_0.04(200bp).Rda")
Xgboost <- non_promoter_enhancers

enhancers_data <- import("/Users/hmz453/Documents/PostDoc/Research_Articles/Pancretic_Islets/Codes/H9/Enhancers_data/regulome_clustering_V7_hg38_sorted.bed.gz")

Active_enhancers_1 <- enhancers_data[enhancers_data$name == 'Active_enhancers_I']
Active_enhancers_2 <- enhancers_data[enhancers_data$name == 'Active_enhancers_II']
Active_enhancers_3 <- enhancers_data[enhancers_data$name == 'Active_enhancers_III']

Active_enhancers_1_CTCF <- enhancers_data[enhancers_data$name == 'Active_enhancers_I_CTCF']
Active_enhancers_2_CTCF <- enhancers_data[enhancers_data$name == 'Active_enhancers_II_CTCF']
Active_enhancers_3_CTCF <- enhancers_data[enhancers_data$name == 'Active_enhancers_III_CTCF']

Active_enhancers <- c(Active_enhancers_1, Active_enhancers_2, Active_enhancers_3,
                      Active_enhancers_1_CTCF, Active_enhancers_2_CTCF,Active_enhancers_3_CTCF)

Active_enhancers_gr <- GRanges(seqnames = seqnames(Active_enhancers), ranges=IRanges(start = start(Active_enhancers), end= end(Active_enhancers)), name=Active_enhancers$name)


Inactive_enhancers_1 <- enhancers_data[enhancers_data$name == 'Inactive_enhancers']
Inactive_enhancers_CTCF <- enhancers_data[enhancers_data$name == 'Inactive_enhancers_CTCF']

Inactive_enhancers <- c(Inactive_enhancers_1, Inactive_enhancers_CTCF)


Inactive_enhancers_gr <- GRanges(seqnames = seqnames(Inactive_enhancers), ranges=IRanges(start = start(Inactive_enhancers), end= end(Inactive_enhancers)), name=Inactive_enhancers$name)

Active_promoters_1 <- enhancers_data[enhancers_data$name == 'Active_promoters']
Active_promoters_CTCF <- enhancers_data[enhancers_data$name == 'Active_promoters_CTCF']

Active_promoters <- c(Active_promoters_1, Active_promoters_CTCF)
Active_promoters_gr <- GRanges(seqnames = seqnames(Active_promoters), ranges=IRanges(start = start(Active_promoters), end= end(Active_promoters)), name=Active_promoters$name)


motif_list <- list(
  CNN = CNN,
  LR = LR,
  XGBoost = Xgboost,
  Active_enhancers = Active_enhancers_gr,
  Inactive_enhancers = Inactive_enhancers_gr,
  Active_promoters = Active_promoters_gr
)

options(scipen = 999)

model_names <- names(motif_list)

# Initialize overlap matrix
overlap_mat <- matrix(
  0,
  nrow = length(motif_list),
  ncol = length(motif_list),
  dimnames = list(model_names, model_names)
)

# Compute pairwise overlaps
for (i in seq_along(motif_list)) {
  for (j in seq_along(motif_list)) {
    overlap_mat[i, j] <- length(
      intersect(motif_list[[i]], motif_list[[j]])
    )
  }
}

# Total counts per set
total_counts <- sapply(motif_list, length)

# Axis labels with counts
axis_labels <- paste0(
  model_names,
  "\n(n = ", total_counts, ")"
)

library(pheatmap)

# Extract set sizes automatically from diagonal
set_sizes <- diag(overlap_mat)

# Compute percent overlap matrix automatically
overlap_pct <- outer(
  set_sizes, set_sizes,
  FUN = function(x, y) overlap_mat / pmin(x, y)
) * 100

pdf(
  "/Users/hmz453/Documents/PostDoc/Research_Articles/Pancretic_Islets/Codes/H9/Analysis/Final_Plots/Models_regulome_overlap_Heatmap_200bp.pdf",
  width = 5,
  height = 4
)

p <- pheatmap(
  overlap_pct,
  display_numbers = overlap_mat,   # keep counts inside cells
  number_color = "black",
  fontsize_number = 7,
  fontsize_row = 8,
  fontsize_col = 8,
  labels_row = axis_labels,
  labels_col = axis_labels,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  breaks = seq(0,100,length.out=101),
  legend_breaks = seq(0,100,20),
  legend_labels = paste0(seq(0,100,20), "%"),
  color = colorRampPalette(c(
    "#FFF0F6",
    "#FADADD",
    "#E8BFD8",
    "#C77DBB"
  ))(100),
  main = "Overlap between different classes",
  border_color = "grey70"
)

print(p)
dev.off()

### Upset Plot #################################################

options(scipen = 999)

enhancers_set <- list(
  CNN = CNN,
  LR = LR,
  XGBoost = Xgboost,
  Active_enhancers = Active_enhancers_gr,
  Inactive_enhancers = Inactive_enhancers_gr,
  Active_promoters = Active_promoters_gr
)

all_regions <- reduce(unlist(GRangesList(enhancers_set)))

binary_matrix <- matrix(0, nrow = length(all_regions), ncol = length(enhancers_set))
colnames(binary_matrix) <- names(enhancers_set)

for (i in 1:length(enhancers_set)) {
  # Find overlaps between all regions and each peak set
  overlaps <- findOverlaps(all_regions, enhancers_set[[i]])
  binary_matrix[queryHits(overlaps), i] <- 1
}

# Convert to data frame for UpSetR
binary_df <- as.data.frame(binary_matrix)
upset_plot <- upset(binary_df,
                    sets = names(enhancers_set),
                    main.bar.color = "black",   # deep plum (vertical bars)
                    sets.bar.color = "black",   # teal (horizontal bars)
                    order.by = "freq",
                    empty.intersections = "on",
                    mainbar.y.label = "Intersection Size",
                    sets.x.label = "Peaks per Histone Mark",
                    mb.ratio = c(0.7, 0.3),
                    point.size = 2.8,
                    line.size = 1.2,
                    text.scale = c(1.3, 1.3, 1, 1, 1.5, 1.5))

# Save as PDF
pdf("/Users/hmz453/Documents/PostDoc/Research_Articles/Pancretic_Islets/Codes/H9/Analysis/Final_Plots/ML_Models_predicted_enhancers_upset_plot(200bp).pdf", width = 20, height = 8)
print(upset_plot)
dev.off()





