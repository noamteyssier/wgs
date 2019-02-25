library(tidyverse)
library(reshape2)
library(RColorBrewer)

meta <- read_tsv('/home/noam/bin/visuals/WGS/sWGA_merge/meta.tab.txt')
metrics <- read_tsv('/home/noam/bin/visuals/WGS/sWGA_merge/sampleMetrics.tab.txt')
covHist <- read_tsv('/home/noam/bin/visuals/WGS/sWGA_merge/coverageHistograms.tab.txt')

# merge meta with main
sampleMetrics <- left_join(metrics, meta, by='SAMPLE_NAME')


# add a count of the number of sWGA
sampleMetrics <- sampleMetrics %>%
  mutate(COUNT_SWGA = count.fields(textConnection(SWGA_COMBINATION), sep = '_'))

# gather Percentile information into long form
quantileInfo <- melt(sampleMetrics,
  id.vars=c('EXTRACTION_METHOD', 'SWGA_COMBINATION', 'PARASITE_DENSITY')) %>%
  filter(grepl('PCT', variable)) %>%
  filter(!grepl('EXC', variable)) %>%
  mutate(value = as.double(value))


# calculate binary percentage of coverage
depthBool <- covHist %>%
  mutate(depth_bool = case_when(
    DEPTH == 0 ~ 'no_depth',
    DEPTH > 0 ~ 'depth'
    )) %>%
    group_by(SAMPLE_NAME, depth_bool) %>%
    summarise(m = sum(COUNT))  %>%
    group_by(SAMPLE_NAME) %>%
    spread(key = depth_bool, value=m) %>%
    mutate(pcNotCovered = no_depth / (no_depth + depth)) %>%
    left_join(meta, by='SAMPLE_NAME')


# calculate binned percentage of coverage
binnedHist <- covHist %>%
  mutate(binnedDepth = case_when(
    DEPTH == 0 ~ 'zero',
    DEPTH > 0 & DEPTH < 100 ~ '1-TO-100',
    DEPTH >= 100 & DEPTH < 250 ~ '100-TO-250',
      DEPTH == 250 ~ '250',
    )) %>%
  group_by(SAMPLE_NAME, binnedDepth) %>%
  summarise(binTotal = sum(COUNT))

# calculate total base total // binFractions // and merge with meta
binnedHist <- binnedHist %>%
  group_by(SAMPLE_NAME) %>%
  summarise(baseTotal = sum(binTotal)) %>%
  left_join(binnedHist) %>%
  mutate(binFraction = binTotal / baseTotal) %>%
  left_join(meta, by = 'SAMPLE_NAME')


# plot percentiles of coverage
percentilesCov <- ggplot(data = quantileInfo %>%
  filter(EXTRACTION_METHOD != 'Neg' & EXTRACTION_METHOD != 'chelex' ),
  aes(x = variable, y = value, colour = SWGA_COMBINATION)) +
  geom_point(alpha = 0.7) +
  geom_line(aes(group = SWGA_COMBINATION))+
  facet_grid(EXTRACTION_METHOD~PARASITE_DENSITY) +
  theme(
    axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 1)
  ) +
  scale_color_manual(values = brewer.pal(7, 'Set1'))



# Binary Coverage
pcNotCov <- ggplot(data = depthBool %>%
  filter(PARASITE_DENSITY != 'Neg') %>%
  filter(EXTRACTION_METHOD != 'chelex'),
  aes(x = PARASITE_DENSITY, y = (pcNotCovered), colour=SWGA_COMBINATION)) +
  geom_point() +
  geom_line(aes(group = SWGA_COMBINATION)) +
  facet_wrap(~EXTRACTION_METHOD, scales='free') +
  scale_color_manual(values =brewer.pal(7, 'Set1'))

ggsave('percentilesCovered.pdf', percentilesCov)
ggsave('pcNotCovered.pdf', pcNotCov)





# Binned histogram
unique(binnedHist$binnedDepth)
binnedHist$binnedDepth <- factor(binnedHist$binnedDepth, levels =
  c('zero', '1-TO-100', '100-TO-250', '250'))


unBinned <- covHist %>%
  left_join(meta) %>%
  filter(EXTRACTION_METHOD != 'Neg') %>%
  mutate(DEPTH = as.integer(DEPTH))

# Sorbet colored depth
ggplot(data = unBinned,
  aes(x = DEPTH , y = log10(COUNT), colour=SWGA_COMBINATION)) +
  # scale_x_continuous(breaks = seq(0, 250, 100)) +
  geom_line(aes(group = SWGA_COMBINATION), alpha=0.7) +
  facet_wrap(~EXTRACTION_METHOD) +
  scale_color_manual(values = brewer.pal(7, 'Set2'))


display.brewer.all()
