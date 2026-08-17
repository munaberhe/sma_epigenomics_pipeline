.libPaths(c("~/R/library", .libPaths()))
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
setwd("/gpfs/scratch/bt25018/sma_epigenomics_pipeline")
OUT <- "results/thesis_figures"

COND_COLS <- c(Scramble_CTRL="#6B7280", ASO_CTRL="#1F3A5F",
               Scramble_VPA="#F0A500", ASO_VPA="#C0392B")
N_FLANK_BINS <- 20; N_BODY_BINS <- 60; N_TOTAL_BINS <- 100

df <- read.csv(file.path(OUT, "metagene_profile_data.csv"))
df$condition <- factor(df$condition, levels=names(COND_COLS))

LABELS <- c(Scramble_CTRL="Scramble + vehicle", ASO_CTRL="ASO1 + vehicle",
            Scramble_VPA="Scramble + VPA", ASO_VPA="ASO1 + VPA")
df$condition_label <- factor(LABELS[as.character(df$condition)], levels=LABELS[names(COND_COLS)])
names(COND_COLS) <- LABELS[names(COND_COLS)]

tss_bin <- N_FLANK_BINS + 1
body_plateau_bins <- (N_FLANK_BINS + 10):(N_FLANK_BINS + N_BODY_BINS - 10)
dip_stats <- df %>%
  group_by(condition) %>%
  summarise(body_plateau = mean(meth[bin %in% body_plateau_bins], na.rm=TRUE),
            tss_min = min(meth[bin %in% (tss_bin-2):(tss_bin+2)], na.rm=TRUE),
            dip_depth = body_plateau - tss_min, .groups="drop")
write.csv(dip_stats, file.path(OUT, "metagene_TSS_dip_depth.csv"), row.names=FALSE)
message("TSS dip depth:"); print(dip_stats)

y_min <- floor(min(df$meth, na.rm=TRUE)/5)*5
y_max <- ceiling(max(df$meth, na.rm=TRUE)/5)*5

p <- ggplot(df, aes(x=bin, y=meth, colour=condition_label)) +
  geom_line(linewidth=1) +
  geom_vline(xintercept=N_FLANK_BINS+0.5, linetype="dashed", colour="grey40") +
  geom_vline(xintercept=N_FLANK_BINS+N_BODY_BINS+0.5, linetype="dashed", colour="grey40") +
  annotate("text", x=N_FLANK_BINS+0.5, y=y_max*0.99, label="TSS", size=3.5, hjust=1.1) +
  annotate("text", x=N_FLANK_BINS+N_BODY_BINS+0.5, y=y_max*0.99, label="TES", size=3.5, hjust=-0.1) +
  scale_colour_manual(values=COND_COLS, name=NULL) +
  scale_x_continuous(breaks=c(1, N_FLANK_BINS+1, N_FLANK_BINS+N_BODY_BINS+1, N_TOTAL_BINS),
                     labels=c("-2kb", "TSS", "TES", "+2kb")) +
  scale_y_continuous(breaks=seq(y_min, y_max, 5), expand=expansion(mult=c(0.02, 0.05))) +
  labs(x="Gene body is length-scaled between TSS and TES; flanks are fixed-width (2kb)",
       y="Mean CpG methylation (%)") +
  theme_classic(base_size=13) +
  theme(legend.position="top", axis.title.x=element_text(size=9, colour="grey40"))

ggsave(file.path(OUT, "Fig_metagene_profile.pdf"), p, width=8, height=5.5, device=cairo_pdf)
message("Saved")
