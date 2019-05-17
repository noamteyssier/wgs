library(tidyverse)
library(ggbeeswarm)
library(reshape2)
library(wesanderson)

setwd("~/projects/visuals/WGS/covStats/")
stats_premeta <- read_tsv("stats.tab.txt")
meta <- read_tsv("meta.tab.txt")

stats <- left_join(stats_premeta, meta, by="sampleName")

ext_stats <- stats %>%
  filter(extractionMethod != "na" & swgaMethod != 'na') %>%
  filter(owner == 'GH')


###############################################################################
#Mean and Median Coverage Analysis
###############################################################################
# extration vs MEAN_COVERAGE by swgaMethod
ggplot(data=ext_stats, aes(x = extractionMethod, y = MEAN_COVERAGE)) +
  geom_beeswarm(aes(colour=density)) +
  facet_grid(~swgaMethod) +
  scale_color_manual(values=wes_palette(n=4, name='IsleofDogs1')) +
  labs(title='Extraction vs MEAN_COVERAGE by sWGA Method')

# extraction + swga vs MEAN_COVERAGE
ggplot(data = ext_stats, aes(x = density, y = MEAN_COVERAGE)) +
  geom_beeswarm(aes(colour = swgaMethod, shape = extractionMethod),
    size = 2, dodge.width=0.2) +
  scale_color_manual(values=wes_palette(n=3, name='Darjeeling1')) +
  labs(title='MEAN_COVERAGE')

# extraction + swga vs MEDIAN_COVERAGE
ggplot(data = ext_stats, aes(x = density, y = MEDIAN_COVERAGE)) +
  geom_beeswarm(aes(colour = swgaMethod, shape = extractionMethod),
    size = 2, dodge.width=0.2) +
  scale_color_manual(values=wes_palette(n=3, name='Darjeeling1')) +
  labs(title='MEDIAN_COVERAGE')
###############################################################################


## Comparison of sWGA Method with Percentiles
pctStats <- stats %>%
  filter(owner == 'GH')

extPcts <- melt(pctStats, id.vars=c('extractionMethod', 'swgaMethod', 'density')) %>%
  filter(grepl('PCT', variable)) %>%
  filter(grepl('X$', variable)) %>%
  mutate(value = as.double(value))

###############################################################################
# Percentiles and Exclusion Analysis
###############################################################################
## Percentiles of Coverage by Density
ggplot(data = extPcts %>% filter(extractionMethod == 'qia'), aes(x = variable, y = value, colour = swgaMethod))+
  geom_point(alpha=0.5) +
  geom_hline(yintercept = 0.05) +
  facet_wrap(~density) +
  scale_color_manual(values=wes_palette(n=3, name='Darjeeling1')) +
  theme(axis.text.x = element_text(size = 11, angle = 270)) +
  labs(title='Percent at X Coverage by Density')

## Percent of Total Excluded Reads
ggplot(data = ext_stats %>% filter(grepl('qia', extractionMethod)) ,
  aes(x = density, y = PCT_EXC_TOTAL, colour = swgaMethod)) +
  geom_point() +
  scale_color_manual(values=wes_palette(n=3, name='Darjeeling1')) +
  facet_grid(~extractionMethod) +
  labs(title='Percent Excluded of Total Excluded Reads')

## Percent of excluded reads by extraction method at densities
ggplot(data = ext_stats %>% filter(grepl('qia', extractionMethod)) ,
  aes(x = density, y = PCT_EXC_UNPAIRED, colour = swgaMethod)) +
  geom_point() +
  scale_color_manual(values=wes_palette(n=3, name='Darjeeling1')) +
  facet_grid(~extractionMethod) +
  labs(title='Percent Excluded of Unpaired Reads')
###############################################################################

bhStats <- stats %>%
  filter(owner == 'BH')

###############################################################################
# Comparison of Repair vs No Repair in Qiagen and Chelex
###############################################################################
ggplot(data= bhStats %>% filter(extractionMethod != 'na'), aes(x = extractionMethod, y = MEAN_COVERAGE))+
  geom_point(aes(colour = repairMethod))
