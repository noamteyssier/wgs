library(tidyverse)

fs <- read_tsv("data/flagstat.tab")

fs <- fs %>%
  mutate(pc_mapped = mapped / total) %>%
  unite(sampleName, -total, -mapped, -paired, -pc_mapped) %>%
  left_join(fs) %>%
  filter(!grepl('Neg', sampleName))


# Comparison of Percentage Mapped by Condition
agg <- fs %>%
  group_by(density, extraction, digest, swga) %>%
  summarise(m = mean(pc_mapped), s = sd(pc_mapped))

# percent mapped bar plots
ggplot(agg, aes(
      x = interaction(swga, extraction),
      y = m,
      fill = digest)) +
    geom_bar(stat='identity', position='dodge') +
    facet_wrap(~density) +
    theme_classic() +
    scale_fill_brewer(palette='Set2') +
    ylab('percent mapped')

# percent mapped densities
ggplot(fs, aes(x = pc_mapped, fill=digest)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~density)


# total reads
ggplot(fs, aes(x = interaction(swga, extraction), y = total, fill=digest)) +
  geom_bar(stat='identity', position='dodge') +
  facet_wrap(~density) +
  scale_fill_brewer(palette='Set2')

# total read densities
ggplot(fs, aes(x = total, fill = digest)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~density) +
  scale_fill_brewer(palette='Set2')


fs %>%
  select(-mapped, -paired, -sampleName, -pc_mapped) %>%
  tidyr::spread(key = rep, value = total) %>%
  mutate(min_pair = pmin(`1`, `2`)) %>%
  gather('rep', 'total_reads', -density, -extraction, -digest, -swga, -min_pair) %>%
  mutate(ds_fraction = min_pair / total_reads)
