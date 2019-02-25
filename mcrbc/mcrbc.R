library(tidyverse)
library(ggpubr)

setwd("~/bin/wgs/mcrbc")

collapsed <- read_tsv("data/original/collapsed.tab.gz")
metrics <- read_tsv("data/original/maskedMetrics.tab")

norm_collapsed <- read_tsv("data/normalized/ds_collapsed.tab.gz")
merged_metrics <- read_tsv("data/normalized/maskedMetrics.tab")

################################################################################
# variance comparison between large volume and small volume sWGA
################################################################################

wide_collapsed <- collapsed %>%
  spread(key = rep, val = depth)

norm_wide_collapsed <- collapsed %>%
  spread(key = rep, val = depth)

correlation_plot <- ggplot(wide_collapsed, aes(x = log10(`1`), y=log10(`2`))) +
  geom_point(shape = 21, size = 0.5, alpha = 0.5, aes(fill = interaction(amp, digest))) +
  facet_grid(digest~amp) +
  stat_cor() +
  theme_classic()

norm_correlation_plot <- ggplot(norm_wide_collapsed, aes(x = log10(`1`), y=log10(`2`))) +
  geom_point(shape = 21, size = 0.5, alpha = 0.5, aes(fill = interaction(amp, digest))) +
  facet_grid(digest~amp) +
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
    -amp, -density, -ext, -digest, -rep
  ) %>%
  filter(grepl('PCT_[0-9]+X', depth)) %>%
  mutate(depth = gsub('PCT_', '', depth)) %>%
  mutate(depth = gsub('X', '', depth) %>% as.integer())

merged_depthPCT <- merged_metrics %>%
  gather(
    key = 'depth', value = 'percentage',
    -amp, -density, -ext, -digest
  ) %>%
  filter(grepl('PCT_[0-9]+X', depth)) %>%
  mutate(depth = gsub('PCT_', '', depth)) %>%
  mutate(depth = gsub('X', '', depth) %>% as.integer())


percentiles <- ggplot(depthPCT, aes(x = depth, y = percentage)) +
  geom_point(
    shape = 21, aes(
      fill = interaction(digest, amp)
      )) +
  geom_line(
    aes(
      group = interaction(interaction(rep, amp), digest),
      colour = interaction(digest, amp))
  ) +
  geom_hline(yintercept = 0.5, linetype = 'dashed') +
  geom_vline(xintercept = 10) +
  scale_y_continuous(breaks = c(seq(0,1,0.1))) +
  scale_x_continuous(breaks = c(seq(0,100,10))) +
  theme_classic()

merged_percentiles <- ggplot(merged_depthPCT, aes(x = depth, y = percentage)) +
  geom_point(
    shape = 21, aes(
      fill = interaction(digest, amp)
      )) +
  geom_line(
    aes(
      group = interaction(amp, digest),
      colour = interaction(digest, amp))
  ) +
  geom_hline(yintercept = 0.5, linetype = 'dashed') +
  geom_vline(xintercept = 10) +
  scale_y_continuous(breaks = c(seq(0,1,0.1))) +
  scale_x_continuous(breaks = c(seq(0,100,10))) +
  theme_classic()

ggsave('plots/percentiles.png', percentiles, width = 10, height = 6)
ggsave('plots/merged_percentiles.png', merged_percentiles, width = 10, height = 6)
