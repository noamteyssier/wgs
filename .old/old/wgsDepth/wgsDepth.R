library(tidyverse)
library(ggplot2)
library(broom)
library(purrr)

depth <- read_tsv("~/projects/visuals/wgsDepth/sampleDepth.tab") 
totalReads <- read_tsv("~/projects/visuals/wgsDepth/totalReads.tab")


?quantile

## calculate quantiles of data
quantiles <- depth %>% 
  nest(-sample) %>% 
  mutate(Quantiles = map(data, ~ quantile(.$depth, probs = seq(0,1,0.1)))) %>% 
  unnest(map(Quantiles, tidy))

# quantiles as matrix
quantiles %>%  
 spread(sample,x)

## ddply version of applying quantile
percentiles <- ddply(depth, .(sample), function(x) quantile(x$depth, probs = seq(0,1,0.1)))


# join quantiles with total read data
df <- left_join(quantiles, totalReads, by = 'sample')

# reorder quantile names
#df$names = factor(df$names, levels = unique(df$names))

## Plotting
ggplot(data = df, aes(x = log10(totalReads), y = x, color = sWGA)) +
  geom_jitter() +
  scale_y_log10() +
  facet_wrap(~names)

df$x
