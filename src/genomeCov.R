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
  ggplot(percentiles_df, aes(x = depth, y = percentage, color=as.factor(swga), linetype=digest)) +
    geom_line(size = 1) +
    facet_wrap(extraction~density, ncol=2) +
    theme_bw() +
    scale_color_manual(values = c("darkgreen", '#965e04')) %>%
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
    facet_wrap(extraction~density, scales='free_x') +
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
percentiles_df <- read_tsv("../data/opt4/merged_metrics.tab")
human_flagstat <- read_tsv("../data/opt4/humanMappingFlagstat.tab")
flagstat <- read_tsv("../data/opt4/flagstat.tab")

# process
percentiles_df <- melt_pc_df(percentiles_df)
flagstat <- melt_flagstat(flagstat)
human_flagstat <- melt_flagstat(human_flagstat)
pc_human_mapping <- percent_mapping(human_flagstat)

human_mapped <- human_flagstat %>%
  filter(statistic == 'mapping') %>%
  mutate(statistic = 'human')
flagstat <- rbind(flagstat, human_mapped)


# plot
percentile_cov_plot <- percentiles_df %>% plot_percentiles()
flagstat_plot <- flagstat %>% plot_flagstat()
human_mapping_plot <- pc_human_mapping %>% plot_mapping()

ggsave("../plots/percentile_coverage.png", percentile_cov_plot, width = 8, height = 6)
ggsave("../plots/pf_flagstat.png", flagstat_plot, width = 8, height = 6)
ggsave("../plots/human_mapping.png", human_mapping_plot, width = 8, height = 6)
