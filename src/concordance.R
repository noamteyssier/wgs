library(tidyverse)
library(scales)

'%!in%' <- function(x,y)!('%in%'(x,y))
errors <- read_tsv("../data/complete_set/concordance.tab")
pf3k_sites <- read_tsv("../data/complete_set/concordance.summary.tab")
single_strains <- c('U51', 'fcr3')
hex = c(
    '#e45c26',
    '#00B26e',
    '#efa89a',
    '#88D5A9',
    '#ffe6de',
    '#bfead0'
  ) %>% rev()

#########################
# pf3k site_concordance #
#########################

pf3k_sim <- pf3k_sites %>%
  group_by(CALL_SAMPLE) %>%
  top_n(1, GENOTYPE_CONCORDANCE) %>%
  extract(
    CALL_SAMPLE,
    into=c('density', 'extraction', 'digest', 'swga', 'rep'),
    regex="([[:alnum:]]+)-([[:graph:]]+)-([[:graph:]])-([[:graph:]]+)-([[:alnum:]])"
  )
grouped_pf3k_sim <- pf3k_sim %>%
  group_by(density, extraction, digest, swga) %>%
  summarise(
    m_concordance = mean(GENOTYPE_CONCORDANCE),
    s_concordance = sd(GENOTYPE_CONCORDANCE)
  )  %>%
  ungroup() %>%
  mutate(density = paste0(density, ' p/uL')) %>%
  add_row(density = '10 p/uL', extraction = 'Che', digest='R', swga='Sof', m_concordance=0.001)


grouped_pf3k_sim$digest <- factor(grouped_pf3k_sim$digest, levels=c('R', 'M'))
grouped_pf3k_sim$swga <- factor(grouped_pf3k_sim$swga, levels=c('Sof', 'Sang'))
full_size_plot <- ggplot(grouped_pf3k_sim,
  aes(
    x = interaction(digest, interaction(swga,extraction)),
    y = m_concordance,
    fill = extraction)) +
    geom_bar(stat='identity', position='dodge', size=.1, color='black', width=1) +
    theme_bw() +
    facet_grid(~density, scales='free', space='free') +
    theme(axis.text.x = element_text(angle=45, hjust=.95, vjust=.95)) +
    scale_fill_manual(values=hex %>% tail(2)) +
    xlab("sWGA and Extraction Condition") +
    ylab("SNP Concordance At Pf3k SNP Sites (N = 1189831)") +
    geom_hline(yintercept = 0.995, linetype='dashed') +
    scale_y_continuous(expand=c(0.001,0.01), limits=c(0,1), oob=rescale_none)
zoom_in_plot <- ggplot(grouped_pf3k_sim,
    aes(
      x = interaction(digest, interaction(swga,extraction)),
      y = m_concordance,
      fill = extraction)) +
  geom_bar(stat='identity', position='dodge', size=.1, color='black', width=1) +
  theme_bw() +
  facet_grid(~density, scales='free', space='free') +
  theme(axis.text.x = element_text(angle=45, hjust=.95, vjust=.95)) +
  scale_fill_manual(values=hex %>% tail(2)) +
  xlab("sWGA and Extraction Condition") +
  ylab("SNP Concordance At Pf3k SNP Sites (N = 1189831)") +
  scale_y_continuous(expand=c(0.0000001,0.000001), limits=c(0.995,1), oob=rescale_none)

ggsave("../plots/full_size.pdf", full_size_plot, width = 8, height = 6)
ggsave("../plots/full_size.png", full_size_plot, width = 8, height = 6)
ggsave("../plots/zoom_in.pdf", zoom_in_plot, width = 4, height = 3)
ggsave("../plots/zoom_in.png", zoom_in_plot, width = 4, height = 3)


#############################
# Non Reference Concordance #
#############################

similarity <- errors %>%
  filter(
    (s1 %in% single_strains) & (s2 %!in% single_strains)) %>%
  group_by(s2) %>%
  mutate(concordance = 1 - pc_error) %>%
  top_n(1, concordance)
grouped_sims <- similarity %>%
  extract(
    s2,
    into=c('density', 'extraction', 'digest', 'swga', 'rep'),
    regex="([[:alnum:]]+)-([[:graph:]]+)-([[:graph:]])-([[:graph:]]+)-([[:alnum:]])"
  ) %>%
  group_by(density, extraction, digest, swga) %>%
  summarise(m_concordance = mean(concordance)) %>%
  ungroup() %>%
  mutate(density = paste0(density, ' p/uL')) %>%
  add_row(density = '10 p/uL', extraction = 'Che', digest='R', swga='Sof', m_concordance=0.001)

grouped_sims$digest <- factor(grouped_sims$digest, levels=c('R', 'M'))
grouped_sims$swga <- factor(grouped_sims$swga, levels=c('Sof', 'Sang'))

concordance_plt <- ggplot(grouped_sims,
    aes(
      x = interaction(digest, interaction(swga,extraction)),
      y = m_concordance,
      fill = extraction)) +
  geom_bar(stat='identity', position='dodge', size=.1, color='black', width=1) +
  theme_bw() +
  facet_grid(~density, scales='free_x', space='free') +
  theme(axis.text.x = element_text(angle=45, hjust=.95, vjust=.95)) +
  scale_fill_manual(values=hex %>% tail(2)) +
  xlab("sWGA and Extraction Condition") +
  ylab("SNP Concordance to Reference (VS1 / fcr3)") +
  scale_y_continuous(breaks = seq(0, 1, .1), expand=c(0.001,0.01), limits=c(0,1))

ggsave("../plots/concordance.pdf", concordance_plt, width = 8, height = 6)
ggsave("../plots/concordance.png", concordance_plt, width = 8, height = 6)
