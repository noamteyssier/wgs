library(tidyverse)

## Processing
melt_pc_df <- function(percentiles_df){
  # function to melt percentiles df to longform
  percentiles_df %>%
    gather("depth", "percentage", -digest, -extraction, -swga, -density) %>%
    filter(grepl('PCT_[[:digit:]]+X', depth)) %>%
    mutate(
      depth = gsub('PCT_', '', depth),
      depth = as.numeric(gsub('X', '', depth))) %>%
    group_by(density, extraction, digest, swga, depth) %>%
    summarise(percentage = mean(percentage)) %>%
    ungroup() %>%
    return()
}
melt_flagstat <- function(flagstat){
  # function to melt total-mapped-paired statistic of flagstat
  flagstat %>%
    mutate(sampleName=paste(density, extraction, digest, swga, rep,  sep='-')) %>%
    filter(!grepl('Neg', sampleName)) %>%
    gather('statistic', 'value', -sampleName, -rep, -density, -swga, -extraction, -digest) %>%
    return()
}
percent_mapping <- function(melted_flagstat){
  # function to calculate flagstat pc_mapped from melted dataframe
  melted_flagstat %>%
    spread(statistic, value) %>%
    mutate(pc_mapped = mapping / total) %>%
    return()
}

## Plotting
plot_percentiles <- function(percentiles_df){
  # plot percentiles of coverage by group
  percentiles_df <- percentiles_df %>%
    mutate(
      swga = ifelse(swga == 'Sof', '6A + 10A', '10A'),
      digest = ifelse(digest == 'M', 'McrBC', 'No McrBC'),
      `sWGA primers` = as.factor(swga)
    )
  ggplot(percentiles_df, aes(x = depth, y = percentage, color=interaction(`sWGA primers`, extraction), linetype=digest)) +
    geom_line(size = 1) +
    facet_wrap(~density, ncol=1, scales='free') +
    theme_bw() +
    scale_y_continuous(limits=c(0,1))  +
    scale_color_brewer(palette="paired") %>%
    return()
}
plot_flagstat <- function(flagstat){
  # plot flagstat statistics by group
  flagstat$statistic <- factor(flagstat$statistic, levels=c('human', 'mapped', 'total'))

  ggplot(flagstat %>%
      filter(statistic != 'paired') %>%
      unite("sampleName", c(digest, swga, rep), remove=FALSE, sep='-'),
      aes(x = sampleName, y = value/2, fill=digest)) +
    geom_bar(stat='identity', position='identity', alpha = 0.5, aes(color=statistic), size=0.8) +
    facet_wrap(extraction~density, scales='free_x', ncol=3) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10)) +
    scale_fill_manual(values=c("firebrick4", "peru")) +
    scale_color_manual(values=c("#28377f", "#2c3f18", "black")) +
    ylab("Reads Mapped") +
    xlab("Sample") %>%
    return()
}
plot_mapping <- function(pc_flagstat){
  # plot percent mapping by group
  ggplot(pc_flagstat, aes(x = as.factor(density), y = pc_mapped, fill=digest, color=as.factor(rep))) +
    geom_bar(stat='identity', position='dodge') +
    facet_wrap(extraction~swga) +
    theme_classic() +
    xlab("Parasite Density") +
    ylab("Percent of Reads Mapping to Human") +
    scale_color_brewer(palette='Greys') +
    scale_fill_brewer(palette='Reds')
}

# load in data
percentiles_df <- read_tsv("../data/complete_set/metrics.tab")
flagstat <- read_tsv("../data/complete_set/flagstat.tab") %>%
  mutate(s_name = paste0(density, '-',extraction, '-',digest, '-',swga, '-',rep)) %>%
  mutate(
    falciparum = mapped_falc / total,
    human = mapped_human / total,
    unknown = 1 - (falciparum + human)
  ) %>%
  select(-mapped_falc, -mapped_human) %>%
  gather('organism', 'pc_mapped', -total, -density, -extraction, -digest, -swga, -rep, -s_name)
original_runs <- read_tsv("../data/complete_set/original_run_results.tab")

# Flagstat processing
summarised_flagstat <- flagstat %>%
  group_by(density, extraction, digest, organism) %>%
  summarise(
    mean_pc = mean(pc_mapped),
    std_pc = sd(pc_mapped)
  )
summarised_flagstat$organism <- factor(summarised_flagstat$organism, levels=c('unknown', 'human', 'falciparum'))
summarised_flagstat$digest <- factor(summarised_flagstat$digest, levels=c('R', 'M'))
pc_mapped_plt <- ggplot(summarised_flagstat, aes(x = interaction(extraction, digest), y = mean_pc, fill=organism)) +
  geom_bar(stat='identity', position='stack') +
  geom_errorbar(
    data = summarised_flagstat %>% filter(organism == 'falciparum'),
    aes(ymin=mean_pc - std_pc, ymax = mean_pc + std_pc),
    width=.2
    ) +
  facet_wrap(~density, scales='free_x') +
  theme_minimal() +
  scale_fill_brewer(palette = 'BuGn') +
  xlab('Extraction and Digest Condition') +
  ylab('Percentage Mapped') +
  ggtitle('Percentage Mapped to Organism by Parasite Density')
ggsave('../plots/pf_flagstat.pdf', pc_mapped_plt, width=8, height=8)
ggsave('../plots/pf_flagstat.png', pc_mapped_plt, width=8, height=8)

grouped_flags <- original_runs %>%
  mutate(pc_mapped = mapped / totalReads) %>%
  group_by(extraction, density, swga) %>%
  summarise(pc_mapped = mean(pc_mapped))

ggplot(grouped_flags, aes(y = as.factor(density), x = extraction, fill=pc_mapped)) +
  geom_tile() +
  facet_wrap(~swga) +
  theme_minimal() +
  ylab('Parasite Density') +
  xlab('Extraction Method') +
  theme(axis.text.x = element_text(angle = 45, vjust = .95, hjust=.95, size = 10)) +
  scale_fill_gradient(low = "gray80", high = "steelblue")


# process
percentiles_df <- melt_pc_df(percentiles_df)
# plot
percentile_cov_plot <- percentiles_df %>% plot_percentiles()

ggsave("../plots/percentile_coverage.pdf", percentile_cov_plot, width = 8, height = 6)
ggsave("../plots/percentile_coverage.png", percentile_cov_plot, width = 8, height = 6)
