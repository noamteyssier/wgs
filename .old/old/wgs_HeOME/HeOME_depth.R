library(tidyverse)
library(ggplot2)
library(wesanderson)
library(ggbeeswarm)

setwd("~/projects/visuals/wgs_HeOME")

raw <- read_tsv("xV_HeOME.tab")
chromLengths <- read_tsv("chromLengths.tab")
meta <- read_tsv("xV_meta.tab")
heome <- left_join(raw, meta, by = 'sample') %>%
  arrange(sample, depth) %>%
  filter(sWGA == 'yes') %>%
  mutate(window = paste(chrom, start, end, direction, sep = '.'))


## calculate median depth at each window
medianDepth <- heome %>%
  group_by(window, totalCoverage) %>%
  summarise(med = median(depth))

## create dataframe of mean/median/standard deviation of depth at each window
windowStats <- heome %>%
  group_by(window, totalCoverage) %>%
  summarise(mean = mean(depth), std = sd(depth), median = median(depth))

## calculate the number of samples each window is found at
windowCoverage <- heome %>%
  group_by(window, totalCoverage) %>%
  filter(depth > 5) %>% # with a depth greater than 5
  count()

## determine which windows perform poorly in all 6 samples
poor_windows <- heome %>%
  filter(depth < 10) %>%
  group_by(window, totalCoverage) %>%
  add_count() %>%
  filter(n == 6)

length((poor_windows %>% filter(totalCoverage == 0))$window) # 25 windows fall out on each sample

## calculate log10 variance of each window
windowVar <- heome %>%
  group_by(window, totalCoverage) %>%
  mutate(v = var(depth), m = mean(depth)) %>%
  mutate(dev = v - m)

## mean of all window variance
empSD <- windowVar$ls %>% mean(na.rm = T)

## Quantiles
quantile((medianDepth %>% filter(totalCoverage == 0))$med, probs = seq(0.1,.9, .1)) # quantiles for each window w/o total coverage
quantile((medianDepth %>% filter(totalCoverage == 1))$med, probs = seq(0.1,.9, .1)) # quantiles for each window w/ total coverage

################################################################################
###Plotting#####################################################################
################################################################################

## Histograms
medianDepth <- medianDepth %>% mutate(med = ifelse(med == 0, 0.01, med)) # transformation of just zeroes

# median depth of each window
hist(log10((medianDepth %>% filter(totalCoverage == 0))$med),
     xlab = 'Sequencing Depth (log10)', main ='Median Average Depth Across Windows')
# median depth of each window with totalCoverage
hist(log10((medianDepth %>% filter(totalCoverage == 1))$med),
     xlab = 'Sequencing Depth (log10)', main = 'Median Depth with Complete Coverage of Window')



hist(windowCoverage$n)        # sample coverage at each window over 10x
hist(log10(windowVar$v))      # variance of each window



## chromosomal visualization
ggplot(data = chromLengths %>% arrange(chrom), aes(x = chrom, y = length)) +
  geom_bar(stat = 'identity', width = 0.05) +
  geom_point(data = heome,
    aes(x = chrom, y = start,
      size = depth, color = depth), alpha = 0.5) +
  coord_flip()  +
  scale_color_gradient(low = 'blue', high = 'firebrick4') +
  theme_classic()

## sample coverage
ggplot(data = heome %>% filter(grepl('xV', sample)) %>% filter(depth < 4), aes(x = totalReads, y = log10(depth), group = interaction(chrom, start, end, direction))) +
  geom_jitter(size = 1, alpha = 0.5) +
  geom_line()

## depth at each window
ggplot(data = heome, aes(x = window, y = log10(depth), fill = window)) +
  geom_point(alpha = 0.5) +
  # geom_boxplot() +
  guides(fill = F)


## mean vs variance for each window across samples
ggplot(data = windowVar, aes(x = log10(m), y = log10(v))) +
  geom_jitter()

## variance - mean for each window across samples
ggplot(data = windowVar, aes(x = window, y = dev)) +
  geom_jitter()



# count data
# (assume : poisson model w/ underdispersion on meanxSD)
