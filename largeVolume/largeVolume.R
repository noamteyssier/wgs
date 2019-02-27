library(tidyverse)
library(ggpubr)

setwd("~/bin/wgs/largeVolume")

collapsed <- read_tsv("data/original/collapsed.tab.gz")
metrics <- read_tsv("data/original/maskedMetrics.tab")

norm_collapsed <- read_tsv("data/normalized/ds_collapsed.tab.gz")
merged_metrics <- read_tsv("data/normalized/merged_masked.tab")

################################################################################
# variance comparison between large volume and small volume sWGA
################################################################################

wide_collapsed <- collapsed %>%
  spread(key = rep, val = depth)

wide_norm_collapsed <- norm_collapsed %>%
  spread(key = rep, val = depth)

correlation_plot <- ggplot(wide_collapsed, aes(x = log10(`1`), y=log10(`2`))) +
  geom_point(shape = 21, size = 0.5, alpha = 0.5, aes(fill = interaction(swga, type))) +
  facet_grid(type~swga) +
  stat_cor() +
  theme_classic()

norm_correlation_plot <- ggplot(wide_norm_collapsed, aes(x = log10(`1`), y=log10(`2`))) +
  geom_point(shape = 21, size = 0.5, alpha = 0.5, aes(fill = interaction(swga, type))) +
  facet_grid(type~swga) +
  stat_cor() +
  theme_classic()

ggsave('plots/correlation.png', correlation_plot, width = 10, height = 8)
ggsave('plots/norm_correlation.png', norm_correlation_plot, width = 10, height = 8)


################################################################################
# percentile coverage comparison between the two
################################################################################

depthPCT <- metrics %>%
  gather(
    key = 'depth', value = 'percentage',
    -type, -density, -ext, -swga, -rep
  ) %>%
  filter(grepl('PCT_[0-9]+X', depth)) %>%
  mutate(depth = gsub('PCT_', '', depth)) %>%
  mutate(depth = gsub('X', '', depth) %>% as.integer())

merge_depthPCT <- merged_metrics %>%
  gather(
    key = 'depth', value = 'percentage',
    -type, -density, -ext, -swga,
  ) %>%
  filter(grepl('PCT_[0-9]+X', depth)) %>%
  mutate(depth = gsub('PCT_', '', depth)) %>%
  mutate(depth = gsub('X', '', depth) %>% as.integer())


percentiles <- ggplot(depthPCT, aes(x = depth, y = percentage)) +
  geom_point(
    shape = 21, aes(
      fill = interaction(type, swga)
      )) +
  geom_line(
    aes(
      group = interaction(rep, swga),
      colour = interaction(type, swga))
  ) +
  geom_hline(yintercept = 0.5, linetype = 'dashed') +
  geom_vline(xintercept = 10) +
  facet_wrap(~type, ncol=1) +
  theme_classic()

merged_percentiles <- ggplot(merge_depthPCT, aes(x = depth, y = percentage)) +
  geom_point(
    shape = 21, aes(
      fill = interaction(type, swga)
      )) +
  geom_line(
    aes(
      group = interaction(ext, swga),
      colour = interaction(type, swga))
  ) +
  geom_hline(yintercept = 0.5, linetype = 'dashed') +
  geom_vline(xintercept = 10) +
  facet_wrap(~type, ncol=1) +
  theme_classic()

ggsave('plots/percentiles.pdf', percentiles, width = 10, height = 6)
ggsave('plots/merged_percentiles.pdf', merged_percentiles, width = 10, height = 6)
