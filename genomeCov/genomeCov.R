library(tidyverse)
library(ggpubr)
library(wesanderson)

setwd("~/bin/wgs/genomeCov")

collapsed <- read_tsv("data/ds_collapsed.tab.gz")
masked <- read_tsv("data/ds_merged_rep_metric.tab")
flagstat <- read_tsv("data/flagstat.tab")
ss_flagstat <- read_tsv("data/subsampled_flagstat.tab")

###########################################
# summary plots of sequencing and mapping #
###########################################

flagstat <- flagstat %>%
  gather('flagstat', 'count', -density, -extraction, -swga, -rep)
flagstat$flagstat <- factor(flagstat$flagstat, levels = c('total', 'mapped', 'paired'))

mpt <- ggplot(data = flagstat, aes(x = interaction(swga, extraction), y = count, fill=flagstat)) +
  geom_bar(stat='identity', position='dodge') +
  coord_polar() +
  facet_wrap(~density) +
  scale_fill_brewer(palette = "Reds") +
  theme_bw() +
  xlab("") +
  theme(legend.position='bottom')


ss_flagstat <- ss_flagstat %>%
  gather('flagstat', 'count', -density, -extraction, -swga, -rep)
ss_flagstat$flagstat <- factor(ss_flagstat$flagstat, levels = c('total', 'mapped', 'paired'))

ss_mpt <- ggplot(data = ss_flagstat, aes(x = interaction(swga, extraction), y = count, fill=flagstat)) +
  geom_bar(stat='identity', position='dodge') +
  coord_polar() +
  facet_wrap(~density) +
  scale_fill_brewer(palette = "Blues") +
  theme_bw()  +
  xlab("") +
  theme(legend.position='bottom')

#####################################################
# variance comparison between normalized replicates #
#####################################################

wide_collapsed <- collapsed %>%
  select(-snum) %>%
  spread(key = rep, val = depth)

correlation_plot <- ggplot(wide_collapsed, aes(x = log10(`1`), y=log10(`2`))) +
  geom_point(shape = 21, size = 0.5, alpha = 0.5, aes(fill = interaction(swga, ext))) +
  facet_grid(ext~swga~density) +
  stat_cor() +
  theme_classic() +
  xlab("Replicate 1 (log10)") +
  ylab("Replicate 2 (log10)")


##################################################
# Percentile Coverage Comparison Between Samples #
##################################################

depthPCT <- masked %>%
  gather(
    key = 'depth', value = 'percentage',
    -density, -ext, -swga
  ) %>%
  filter(grepl('PCT_[0-9]+X', depth)) %>%
  mutate(depth = gsub('PCT_', '', depth)) %>%
  mutate(depth = gsub('X', '', depth) %>% as.integer())

depth_percentage_plot <- ggplot(depthPCT, aes(x = depth, y = percentage)) +
  geom_point(
    shape = 21, aes(
      fill = interaction(swga, ext)
      )) +
  geom_line(
    aes(
      group = interaction(ext , swga),
      colour = interaction(swga, ext))
  ) +
  geom_hline(yintercept = 0.5, linetype = 'dashed') +
  geom_vline(xintercept = 10) +
  facet_wrap(~density, ncol=1) +
  scale_y_continuous(breaks = c(seq(0,1,0.1))) +
  scale_x_continuous(breaks = c(seq(0,100,10))) +
  theme_classic()


##############
# Save Plots #
##############

ggsave("plots/MappedPairedTotal.pdf", mpt, width = 10, height = 8)
ggsave("plots/SSMappedPairedTotal.pdf", ss_mpt, width = 10, height = 8)
ggsave("plots/correlation_plot.png", correlation_plot, width = 10, height = 8)
ggsave("plots/percentile_plot.pdf", depth_percentage_plot, width = 10, height = 6)
