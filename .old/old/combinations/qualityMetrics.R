library(tidyverse)

setwd('~/bin/visuals/WGS/combinations/')
mapq <- read_tsv("insilico_mapQ.tab")
baseq <- read_tsv("insilico_baseq.tab")
bitflag <- read_tsv("bitflag.hist")
ds_bitflag <- read_tsv("ds_bitflags.tab")


ggplot(baseq %>%  filter(baseQuality > 20), aes(x = baseQuality, y = bq_count, colour = interaction(s1, s2))) +
  geom_bar(stat = 'identity', position = 'dodge') +
  facet_wrap(ext~density)

ggplot(mapq, aes(x = mapQuality, y = log10(mq_count), colour = density)) +
  geom_bar(stat='identity', position = 'dodge') +
  facet_wrap(ext~s1~s2, scales = 'free')


perfectMapping <- c(99, 147, 83, 163)
mate_unmapped <- c(73, 133, 89, 121, 165, 181, 101, 117, 153, 185, 69, 137)
both_unmapped <- c(77,141)
incorrectOrientation <- c(67,131,115,179)
uniquelyWrongSize <- c(81,161,97,145,65,129,113,177)

# group flags into status types
bitflag <- bitflag %>%
  mutate(flagStatus = case_when(
    flag %in% perfectMapping ~ 'perfectMapping',
    flag %in% mate_unmapped ~ 'mate_unmapped',
    flag %in% both_unmapped ~ 'both_unmapped',
    flag %in% incorrectOrientation ~ 'incorrectOrientation',
    flag %in% uniquelyWrongSize ~ 'uniquelyWrongSize',
    flag > 2000 ~ 'multimap'
    ))

ds_bitflag <- ds_bitflag %>%
  mutate(flagStatus = case_when(
    flag %in% perfectMapping ~ 'perfectMapping',
    flag %in% mate_unmapped ~ 'mate_unmapped',
    flag %in% both_unmapped ~ 'both_unmapped',
    flag %in% incorrectOrientation ~ 'incorrectOrientation',
    flag %in% uniquelyWrongSize ~ 'uniquelyWrongSize',
    flag > 2000 ~ 'multimap'
    ))

# summarize normalized frequencies of flag statuses
bitflag <- bitflag %>%
  group_by(ext, s1, s2, density) %>%
  summarise(total = sum(count)) %>%
  left_join(bitflag)

ds_bitflag <- ds_bitflag %>%
  group_by(ext, s1, s2, density) %>%
  summarise(total = sum(count)) %>%
  left_join(ds_bitflag)

bitflag <- bitflag %>%
  group_by(ext, s1, s2, density, flagStatus) %>%
  mutate(statusPercentage = sum(count) / total)

ds_bitflag <- ds_bitflag %>%
  group_by(ext, s1, s2, density, flagStatus) %>%
  mutate(statusPercentage = sum(count) / total)


## Comparison of flag status in different conditions
ggplot(bitflag,
    aes(x = flagStatus, y = statusPercentage, colour = interaction(s1, s2))) +
  geom_bar(stat='identity', position = 'dodge') +
  theme(axis.text.x = element_text(angle = 45)) +
  facet_wrap(ext~density)

ggplot(ds_bitflag,
    aes(x = flagStatus, y = statusPercentage, colour = interaction(s1, s2))) +
  geom_bar(stat='identity', position = 'dodge') +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_wrap(ext~density)



full <- bitflag %>%
  left_join(mapq) %>%
  left_join(baseq) %>%
  mutate(
    mqpc = mq_count / total,
    bqpc = bq_count / total
  )


wide_bitflag <- bitflag %>%
  select(-flag, -count) %>%
  unique() %>%
  spread(key = flagStatus, val = statusPercentage) %>%
  mutate(
    r_goodbad = perfectMapping / both_unmapped,
    r_multi = multimap / perfectMapping
  )

wide_ds_bitflag <- ds_bitflag %>%
  select(-flag, -count) %>%
  unique() %>%
  spread(key = flagStatus, val = statusPercentage) %>%
  mutate(
    r_goodbad = perfectMapping / both_unmapped,
    r_multi = multimap / perfectMapping
  )

write_tsv(wide_ds_bitflag, 'ds_swgacombometa.tab')

full_bf <- wide_bitflag %>%
  left_join(mapq) %>%
  left_join(baseq) %>%
  mutate(
    mqpc = mq_count / total,
    bqpc = bq_count / total
  )


ggplot(full_bf, aes(x = mqpc, y = bqpc, colour = interaction(s1, s2))) +
  geom_point(alpha = 0.5, shape = 21) +
  facet_wrap(density~ext)


ggplot(full_bf, aes(x = r_goodbad, y = r_multi, fill = interaction(s1, s2))) +
  geom_point(shape = 21) +
  facet_wrap(ext~density)

ggplot(full_bf, aes(x = r_multi, y = mqpc, colour = interaction(s1, s2))) +
  geom_point(shape = 21) +
  facet_wrap(ext~density)
