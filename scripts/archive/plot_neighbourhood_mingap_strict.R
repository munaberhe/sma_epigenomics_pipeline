.libPaths("~/R/library")
library(ggplot2)
library(dplyr)
library(tidyr)

# Copy the uploaded CSV to results dir first
csv_path <- "results/dmr_benchmark/parameter_benchmark_neighbourhood_mingap.csv"
out_dir  <- "results/dmr_benchmark/plots_v2"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

dat <- read.csv(csv_path)

strict <- dat %>%
  filter(mode == "strict", window_size > 0) %>%
  arrange(window_size)

plot_dat <- strict %>%
  select(window_size, n_real, n_scrambled) %>%
  pivot_longer(c(n_real, n_scrambled), names_to = "type", values_to = "n_dmr") %>%
  mutate(type = recode(type, n_real = "Real DMRs", n_scrambled = "Scrambled DMRs"))

p <- ggplot(plot_dat,
    aes(x = window_size, y = n_dmr, colour = type, linetype = type)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.5, shape = 15) +
  scale_x_continuous(breaks = strict$window_size, labels = strict$window_size) +
  scale_colour_manual(values = c("Real DMRs" = "#1f78b4", "Scrambled DMRs" = "#33a02c")) +
  scale_linetype_manual(values = c("Real DMRs" = "solid", "Scrambled DMRs" = "dashed")) +
  labs(title = "Neighbourhood minGap (strict)",
       subtitle = "Real vs scrambled DMR counts across window sizes",
       x = "window size (bp)", y = "Number of DMRs",
       colour = NULL, linetype = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top", legend.direction = "horizontal",
        axis.text = element_text(colour = "black"),
        plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "nb_mingap_strict_counts.pdf"), p, width = 6, height = 4)
message("Done.")
