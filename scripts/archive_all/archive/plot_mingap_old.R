.libPaths("~/R/library")
library(ggplot2)
library(dplyr)
library(tidyr)

df <- read.csv("results/dmr_benchmark/parameter_benchmark_neighbourhood_mingap_final.csv")

df_clean <- df %>%
  rename(minGap = mingap, mode = scramble_method) %>%
  mutate(mode = factor(mode, levels = c("label_swap", "archie", "stratified")))

p_ratio <- ggplot(df_clean, aes(x = minGap, y = ratio, colour = mode)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = c(100, 200, 300, 500, 1000, 2000)) +
  labs(x = "minGap (bp)", y = "Real / scrambled DMR ratio",
       colour = "Null model",
       title = "Neighbourhood (NB) - minGap sweep",
       subtitle = "ASO VPA vs ASO CTRL") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")

ggsave("results/dmr_benchmark/plots_v2/fig_nb_mingap_ratio.pdf", p_ratio, width = 6, height = 4)

df_long <- df_clean %>%
  select(minGap, mode, n_real, n_scrambled) %>%
  pivot_longer(c(n_real, n_scrambled), names_to = "type", values_to = "n_dmrs") %>%
  mutate(type = recode(type, n_real = "Real DMRs", n_scrambled = "Scrambled DMRs"))

p_counts <- ggplot(df_long, aes(x = minGap, y = n_dmrs, colour = type, linetype = mode)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = c(100, 200, 300, 500, 1000, 2000)) +
  scale_colour_manual(values = c("Real DMRs" = "#1f77b4", "Scrambled DMRs" = "#d62728")) +
  labs(x = "minGap (bp)", y = "Number of DMRs",
       colour = "Series", linetype = "Null model",
       title = "Neighbourhood (NB) - real vs scrambled counts") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")

ggsave("results/dmr_benchmark/plots_v2/fig_nb_mingap_counts.pdf", p_counts, width = 6, height = 4)
message("Done.")
