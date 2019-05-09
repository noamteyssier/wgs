library(tidyverse)

df <- read_tsv("data/merged_metrics.tab")

df <- df %>%
  gather("depth", "percentage", -digest, -extraction, -swga, -density) %>%
  filter(grepl('PCT_[[:digit:]]+X', depth)) %>%
  mutate(
    depth = gsub('PCT_', '', depth),
    depth = as.numeric(gsub('X', '', depth)))


ggplot(df, aes(x = depth, y = percentage, color=as.factor(density), linetype=digest)) +
  geom_line(size = 1) +
  facet_wrap(extraction~swga, ncol=2) +
  theme_classic() +
  geom_hline(yintercept=c(0.5, 0.8), linetype=3, size= 0.2) +
  geom_vline(xintercept=40, linetype=3, size = 0.2)
