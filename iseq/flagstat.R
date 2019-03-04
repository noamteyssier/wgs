library(tidyverse)


flag <- read_tsv("filtered_flagstat.tab") %>%
  filter(!grepl('Undetermined', sample))

flag <- flag %>%
  mutate(density_flag = case_when(
    grepl('100-', sample) ~ '100',
    grepl('1k-', sample) ~ '1000',
    grepl('10-', sample) ~ '10',
    grepl('-100', sample) ~ '100',
    grepl('-1k', sample) ~ '1000',
    grepl('-10', sample) ~ '10'
    ))

ggplot(flag, aes(x=sample, y = total, fill=density_flag)) +
  geom_bar(stat='identity') +
  coord_flip()

ggplot(flag, aes(x=sample, y = pcPaired, fill=density_flag)) +
  geom_bar(stat='identity') +
  coord_flip()
  # coord_polar()
  
