library(tidyverse)

'%!in%' <- function(x,y)!('%in%'(x,y))
errors <- read_tsv("../data/complete_set/snpConcordance.tab")

single_strains <- c('U51', 'fcr3')
single_strains_all <- unique(errors$sample_i) %>% tail(7)

similarity <- errors %>%
  filter(sample_i %in% single_strains) %>%
  filter(sample_j %!in% single_strains_all) %>%
  group_by(sample_j) %>%
  mutate(errorRate = 1 - errorRate) %>%
  top_n(1, errorRate)

grouped_sims <- similarity %>%
  extract(
    sample_j,
    into=c('density', 'extraction', 'digest', 'swga', 'rep'),
    regex="([[:alnum:]]+)-([[:graph:]]+)-([[:graph:]])-([[:graph:]]+)-([[:alnum:]])"
  ) %>%
  group_by(density, extraction, digest, swga) %>%
  summarise(m_err = mean(errorRate)) %>%
  ungroup() %>%
  mutate(density = paste0(density, ' p/uL')) %>%
  add_row(density = '10 p/uL', extraction = 'Che', digest='R', swga='Sof', m_err=0.001)


concordance <- ggplot(grouped_sims, aes(x = interaction(swga,extraction), y = m_err * 100, fill = digest)) +
  geom_bar(stat='identity', position='dodge') +
  theme_minimal() +
  facet_wrap(~density) +
  theme(axis.text.x = element_text(angle=60, hjust=1)) +
  scale_fill_brewer(palette='Paired') +
  xlab("sWGA and Extraction Condition") +
  ylab("SNP Concordance to Reference (VS1 / fcr3)") +
  scale_y_continuous(breaks = seq(0, 100, 10))

ggsave("../plots/concordance.pdf", concordance, width = 8, height = 6)
ggsave("../plots/concordance.png", concordance, width = 8, height = 6)
