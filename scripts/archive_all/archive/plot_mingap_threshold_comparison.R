.libPaths("~/R/library")
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

OUT_DIR  <- "results/dmr_benchmark"
PLOT_DIR <- "results/dmr_benchmark/plots_v2"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- minDiff=0.2 data ---
# Label swap from original CSV
ls <- read.csv(file.path(OUT_DIR, "parameter_benchmark_neighbourhood_mingap.csv"))
ls <- ls %>% filter(window_size > 0) %>%
  mutate(scramble_method="label_swap", minDiff="0.2")

# Stratified from v2 checkpoint
ckpt <- readRDS(file.path(OUT_DIR, "checkpoint_neighbourhood_mingap_v2.rds"))
strat <- do.call(rbind, ckpt) %>%
  rename(window_size=window_size) %>%
  mutate(minDiff="0.2")

df_02 <- bind_rows(
  ls %>% select(window_size, mode, n_real, n_scrambled, ratio, scramble_method, minDiff),
  strat %>% select(window_size, mode, n_real, n_scrambled, ratio, scramble_method, minDiff)
)

# --- minDiff=0.4 data (Radu's params) ---
fin <- read.csv(file.path(OUT_DIR, "parameter_benchmark_neighbourhood_mingap_final.csv"))
# Use minSize=100 only for cleaner plot, rename mingap to window_size
df_04 <- fin %>%
  filter(minsize == 100) %>%
  rename(window_size=mingap) %>%
  mutate(minDiff="0.4 (Radu)", mode="strict") %>%
  select(window_size, mode, n_real, n_scrambled, ratio, scramble_method, minDiff)

df_all <- bind_rows(df_02, df_04)

df_all$scramble_method <- recode(df_all$scramble_method,
  label_swap = "Label-swap",
  stratified = "Stratified scramble",
  archie     = "Read count permutation")

df_all$minDiff <- factor(df_all$minDiff,
  levels = c("0.2", "0.4 (Radu)"))

colours_null <- c(
  "Label-swap"             = "#1565C0",
  "Stratified scramble"    = "#2E7D32",
  "Read count permutation" = "#AD1457")

theme_pub <- theme_bw(base_size = 12) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill="white", colour="black"),
        strip.text       = element_text(face="bold"))

# Panel A — real DMR counts, strict only, facet by minDiff
df_strict <- df_all %>% filter(mode == "strict", window_size > 0)

pA <- ggplot(df_strict,
    aes(x=window_size, y=n_real, colour=scramble_method,
        shape=scramble_method, group=scramble_method)) +
  geom_line(linewidth=0.9) +
  geom_point(size=3) +
  facet_wrap(~minDiff, ncol=2, scales="free_y") +
  scale_x_continuous(breaks=c(100,200,500,1000,2000),
                     labels=c("100","200","500","1000","2000")) +
  scale_colour_manual(values=colours_null, name="Null model") +
  scale_shape_manual(values=c(15,17,18), name="Null model") +
  labs(title="A   Real DMR counts — strict threshold",
       subtitle="minDiff=0.2 shows clear signal; Radu's minDiff=0.4 gives near-zero DMRs",
       x="minGap (bp)", y="n real DMRs") +
  theme_pub

# Panel B — signal/noise ratio, strict only
df_strict_ratio <- df_strict %>%
  filter(!is.infinite(ratio), !is.na(ratio))

pB <- ggplot(df_strict_ratio,
    aes(x=window_size, y=ratio, colour=scramble_method,
        shape=scramble_method, group=scramble_method)) +
  geom_hline(yintercept=1, linetype="dashed", colour="grey50") +
  geom_line(linewidth=0.9) +
  geom_point(size=3) +
  facet_wrap(~minDiff, ncol=2, scales="free_y") +
  scale_x_continuous(breaks=c(100,200,500,1000,2000),
                     labels=c("100","200","500","1000","2000")) +
  scale_colour_manual(values=colours_null, name="Null model") +
  scale_shape_manual(values=c(15,17,18), name="Null model") +
  labs(title="B   Signal/noise ratio — strict threshold",
       subtitle="Ratio > 1 = real data calls more DMRs than scrambled",
       x="minGap (bp)", y="Signal/noise ratio (real/scrambled)") +
  theme_pub

combined <- pA / pB + plot_layout(guides="collect") &
  theme(legend.position="bottom")

ggsave(file.path(PLOT_DIR, "nb_mingap_threshold_comparison.pdf"),
       combined, width=12, height=10)
message("Done.")
