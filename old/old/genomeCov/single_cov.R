library(tidyverse)
library(wesanderson)
library(ggpubr)

setwd("~/bin/visuals/WGS/genomeCov")

# read meta
meta_s3 <- read_tsv("covData/wgs3/wgs3_mapStats.tab")
meta_m3 <- read_tsv("covData/wgs3/wgs3_m_mapStats.tab")
pCov <- read_tsv("covData/percentilesCoverage.tab") %>%
  filter(grepl('w3', sample)) %>%
  separate(sample, into = c('run', 'density', 'ext', 'swga', 'rep'), sep = '-', extra = 'drop') %>%
  select(-run)
maskCov <- read_tsv("covData/wgs3/maskedMetrics.tab") %>%
  separate(
    sample,
    into = c('density', 'ext', 'swga', 'rep'),
    sep = '-', extra = 'drop'
  )



# chromosome masking
mask <- read_tsv("covData/pf3d7_regiontypes.tab")

# sample depths
single_3 <- read_tsv("covData/wgs3/wgs3_sampleDepth.tab.gz") %>%
  left_join(meta_s3)
mcrbc_3 <- read_tsv("covData/wgs3/wgs3_mcrbc_sampleDepth.tab.gz") %>%
  left_join(meta_m3)

# separation values
vars <- c('density', 'ext', 'swga', 'rep', 'snum')
m_vars <- c('density', 'ext', 'digest', 'swga', 'rep', 'snum')


# read counts in map classes
meta3 <- meta_s3 %>%
  separate(sample, into = vars, sep = '-', extra = 'drop') %>%
  gather('mapGroup','count', -density, -ext, -swga, -rep, -snum) %>%
  filter(grepl('10', density))

# read counts in map classes for mcrbc
metam3 <- meta_m3 %>%
  separate(sample, into = m_vars, sep = '-', extra = 'drop') %>%
  gather('mapGroup','count', -density, -ext, -swga, -rep, -snum, -digest) %>%
  filter(density != 'Undetermined')

# plotting mapping counts
meta3$mapGroup <- factor(meta3$mapGroup, levels = c('paired','mapped', 'total'))
g <- ggplot(meta3, aes(x = mapGroup, y = count, fill = interaction(ext, swga), colour = rep)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  scale_fill_manual(values=wes_palette(n=4, name="Moonrise2")) +
  scale_color_manual(values=c("white", "white")) +
  facet_wrap(~density) +
  coord_flip() +
  theme_classic() +
  guides(color = F)

ggsave("plots/readCounts.png", g)


# sample coverage tibble
s3 <- single_3 %>%
  separate(sample, into = vars, sep = '-', extra = 'drop') %>%
  filter(grepl('10', density)) %>%
  left_join(mask) %>%
  mutate(pcd = depth / paired)

# mcrbc sample coverage tibble
m3 <- mcrbc_3 %>%
  separate(sample, into = m_vars, sep = '-', extra = 'drop') %>%
  filter(grepl('10', density)) %>%
  left_join(mask) %>%
  mutate(pcd = depth / paired)


# coverage across genome via chromPlot
ggplot(s3 %>%
    filter(chrom == '12') %>%
    filter(density == '1000') %>%
    filter(regionType == 'Core') %>%
    filter(pos >= start & pos <= end),
  aes(x = pos, y = log10(depth), colour = interaction(ext, rep))) +
  geom_line()+
  facet_wrap(ext~swga, ncol = 1) +
  scale_color_manual(values=wes_palette(n=4, name="Cavalcanti1"))

# coverage across mcrbc via chromPlot
ggplot(m3 %>%
    filter(chrom == '07') %>%
    filter(density == '1000'),htp
  aes(x = pos, y = log10(depth), colour = rep)) +
  geom_line() +
  facet_wrap(ext~swga, ncol = 1) +
  scale_color_manual(values=wes_palette(n=2, name="Moonrise2"))


#
unMasked <- pCov %>%
  gather(key = 'pct', value = 'genomeCov', -ext, -swga, -density, -rep ) %>%
  filter(grepl('PCT_[0-9]+X', pct)) %>%
  mutate(pct = gsub('PCT_','', pct)) %>%
  mutate(pct = gsub('X','', pct) %>% as.numeric())

Masked <- maskCov %>%
  gather(key = 'pct', value = 'genomeCov', -ext, -swga, -density, -rep ) %>%
  filter(grepl('PCT_[0-9]+X', pct)) %>%
  mutate(pct = gsub('PCT_','', pct)) %>%
  mutate(pct = gsub('X','', pct) %>% as.numeric())


# plotting percentils of masked and unmasked genome coverage
ggplot(unMasked %>% filter(grepl('10', density)), aes(x = pct, y = genomeCov, fill = interaction(ext, swga))) +
  geom_point(shape = 21, size = 4) +
  theme_classic() +
  geom_vline(xintercept = 10) +
  facet_wrap(~density)

g <- ggplot(Masked %>% filter(grepl('10', density)), aes(x = pct, y = genomeCov, fill = interaction(ext, swga))) +
  geom_point(shape = 21, size = 1, alpha = 0.5) +
  geom_line(size = .8, aes(group = interaction(rep, interaction(ext, swga)), colour = interaction(ext, swga))) +
  theme_classic() +
  geom_vline(xintercept = 10) +
  geom_hline(yintercept = 0.5, linetype = 'dashed') +
  facet_grid(~density) +
  scale_x_discrete(limits = seq(0,100,10)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 7)) +
  xlab("depth") +
  ylab("percent genome covered")

ggsave("plots/maskedPercentiles.png", g, width = 10, height = 7)


# calculate differences between replicates
# s3_diff <-
s3 %>%
  filter(regionType == 'Core') %>%
  filter(pos >= start & pos <= end) %>%
  select(-snum, -total, -mapped, -paired, -start, -end, -depth) %>%
  head()



s3_param <- s3 %>%
  filter(regionType == 'Core') %>%
  filter(pos >= start & pos <= end) %>%
  select(-snum, -total, -mapped, -paired, -start, -end, -depth) %>%
  group_by(
    chrom, pos, density, ext, swga
  ) %>%
  summarise(
    v = var(pcd),
    m = mean(pcd),
    s = sd(pcd)
  ) %>%
  gather(
    key = 'param',
    val = 'value',
    -chrom, -pos, -density, -ext, -swga
  )

ggplot(s3_param %>% filter(density == '1000'), aes(x = log10(value), fill = interaction(swga,ext))) +
  geom_density(position = 'identity', alpha = 0.5) +
  facet_wrap(~param, ncol = 1, scale = 'free')

# densities of log_difference between replicates
ggplot(s3_diff %>% filter(density == '1000'), aes(x = log10(diff), fill = interaction(swga, ext))) +
  geom_density(position = 'identity', alpha = 0.5) +
  facet_wrap(density~ext, ncol = 1) +
  scale_fill_manual(values=wes_palette(n=4, name="Moonrise2"))


# correlation plots between swga-extraction for variance between replicates
g <- ggplot(s3_diff %>% filter(density == '1000'), aes(x = log10(`1`), y=log10(`2`), fill = interaction(ext, swga)))+
  geom_point(size = 0.7, alpha = 0.2, shape = 21) +
  geom_abline() +
  facet_wrap(ext~swga) +
  stat_cor() +
  theme_classic()


ggsave("p1000_variance.png", g)
