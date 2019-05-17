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
    return()
}
melt_flagstat <- function(flagstat){
  # function to melt total-mapped-paired statistic of flagstat
  flagstat %>%
    mutate(sampleName=paste(density, extraction, digest, swga, rep, sep='-')) %>%
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
  ggplot(percentiles_df, aes(x = depth, y = percentage, color=as.factor(density), linetype=digest)) +
    geom_line(size = 1) +
    facet_wrap(extraction~swga, ncol=2) +
    theme_classic() +
    scale_color_brewer(palette = 'Dark2') +
    geom_hline(yintercept=c(0.5, 0.8), linetype=3, size= 0.2) +
    geom_vline(xintercept=40, linetype=3, size = 0.2) %>%
  return()
}
plot_flagstat <- function(flagstat){
  # plot flagstat statistics by group
  flagstat$statistic <- factor(flagstat$statistic, levels=c('paired', 'mapped', 'total'))
  ggplot(flagstat, aes(x = sampleName, y = value/2, fill=digest)) +
    geom_bar(stat='identity', position='identity', alpha = 0.5, aes(color=statistic)) +
    facet_wrap(extraction~swga, scales='free_x') +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10)) +
    scale_color_brewer(palette='RdBu') +
    scale_fill_brewer(palette='Blues') %>%
    return()
}
plot_mapping <- function(pc_flagstat){
  # plot percent mapping by group
  ggplot(pc_flagstat, aes(x = as.factor(density), y = pc_mapped, fill=digest, color=as.factor(rep))) +
    geom_bar(stat='identity', position='dodge') +
    facet_wrap(swga~extraction) +
    theme_classic() +
    xlab("Parasite Density") +
    ylab("Percent of Reads Mapping to Human") +
    scale_color_brewer(palette='Greys') +
    scale_fill_brewer(palette='Reds')
}

# load in data
percentiles_df <- read_tsv("../data/opt4/merged_metrics.tab")
human_flagstat <- read_tsv("../data/opt4/humanMappingFlagstat.tab")
flagstat <- read_tsv("../data/opt4/flagstat.tab")

# process
percentiles_df <- melt_pc_df(percentiles_df)
flagstat <- melt_flagstat(flagstat)
human_flagstat <- melt_flagstat(human_flagstat)
human_mapping <- percent_mapping(human_flagstat)

# plot
percentile_cov_plot <- percentiles_df %>% plot_percentiles()
flagstat_plot <- flagstat %>% plot_flagstat()
human_flagstat_plot <- human_flagstat %>% plot_flagstat()
human_mapping_plot <- human_mapping %>% plot_mapping()

ggsave("../plots/percentile_coverage.png", percentile_cov_plot, width = 8, height = 6)
ggsave("../plots/pf_flagstat.png", flagstat_plot, width = 8, height = 6)
ggsave("../plots/hum_flagstat.png", human_flagstat_plot, width = 8, height = 6)
ggsave("../plots/human_mapping.png", human_mapping_plot, width = 8, height = 6)
