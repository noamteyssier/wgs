library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(broom)

setwd("~/bin/visuals/WGS/genomeCov")

## Sample Read In
# bin depth of original 20M samples
cov <- read_tsv("covData/wgs2/binnedDepth_qia-che_1K-100.tab.gz") %>%
  mutate(chrom = as.numeric(chrom)) %>%
  filter(!is.na(chrom))

# bin depth of double sWGA coverage
doubleCov <- read_tsv("covData/wgs2/doubleSWGA_coverage.tab.gz") %>%
  mutate(chrom = as.numeric(chrom)) %>%
  filter(!is.na(chrom))

filteredDouble <- read_tsv("covData/wgs2/swgaComboCoverage.tab.gz") %>%
  mutate(chrom = as.numeric(chrom)) %>%
  filter(!is.na(chrom))

ds_filt_dbl <- read_tsv("covData/wgs2/downsampledBinned.tab.gz") %>%
  mutate(chrom = as.numeric(chrom)) %>%
  filter(!is.na(chrom))

mask <- read_tsv("covData/pf3d7_regiontypes.tab")

comboMeta <- read_tsv("covData/swgaComboMeta.tab")
ds_comboMeta <- read_tsv("covData/ds_swgacombometa.tab")


# unbin depth of regional coverage (AMA1/GLURP/MSP/TRAP)
protCov <- read_tsv("covData/coverage_highDiverse.tab.gz")


## Summary stats on singleton coverage
sumCov <- cov %>%
  group_by(rep, chrom, swga, ext, den) %>%
  summarise(
    median = median(depth),
    mean = mean(depth),
    iqr = IQR(depth),
    std = sd(depth)
    ) %>%
    left_join(cov, by = c('rep','chrom','swga','ext','den'))

########################
# Chromosomal Coverage #
########################

# coverage across a chromosome and density
h <- ggplot(sumCov %>% filter(chrom == 7 & den == 1000), aes(x = pos, y = log10(depth), fill = interaction(rep,swga))) +
  geom_point(shape = 21, aes(alpha = 0.3)) +
  # geom_line(aes(colour = rep)) +
  geom_hline(aes(yintercept = log10(median), colour = rep)) +
  facet_grid(ext~swga, scale='free_x')

ggsave("chrom7_facet.tiff", h, width = 24, height = 12, device = 'tiff')
ggsave("chrom3_facet.tiff", h, width = 24, height = 12, device = 'tiff')
ggsave('chrom3_facet_100.tiff', h, width = 24, height = 12, device = 'tiff')

sumCov <- sumCov %>%
  group_by(rep, swga, ext, den) %>%
  mutate(index = row_number())

########################
# Full Genome Coverage #
########################

# full genome
massive <- ggplot(sumCov %>% filter(den == 1000), aes(x = index, y = log10(depth), fill = interaction(rep,swga))) +
  geom_point(shape = 21, aes(alpha = 0.3)) +
  # geom_line(aes(colour = rep)) +
  geom_hline(aes(yintercept = log10(median), colour = rep)) +
  facet_grid(ext~swga, scale='free_x')
ggsave("massiveGenome.pdf", massive, width = 30, height = 12)

# full genome violin
sumCov$depth[sumCov$depth == 0 ] <- 0.99
genomeViolin <- ggplot(sumCov %>% filter(den == 1000), aes(x = interaction(swga, ext), y = log10(depth), colour = interaction(swga,ext,rep))) +
  geom_jitter(alpha = 0.1, size = 0.02) +
  geom_violin(draw_quantiles = c(0.25,0.50,0.75), alpha = 0.3) +
  guides(colour = F)
ggsave('genViolin.pdf', genomeViolin)

#
plot(
  log10((sumCov %>% filter(rep == 'run1' & swga == 'Sanger' & ext == 'Chelex'))$depth),
  log10((sumCov %>% filter(rep == 'run2' & swga == 'Sanger' & ext == 'Chelex'))$depth)
)
plot(
  log10((sumCov %>% filter(rep == 'run1' & swga == 'Sof' & ext == 'Qiagen'))$depth),
  log10((sumCov %>% filter(rep == 'run2' & swga == 'Sof' & ext == 'Qiagen'))$depth),
  col = 'skyblue'
)
plot(
  log10((sumCov %>% filter(rep == 'run1' & swga == 'Sanger' & ext == 'Qiagen'))$depth),
  log10((sumCov %>% filter(rep == 'run2' & swga == 'Sanger' & ext == 'Qiagen'))$depth),
  col = 'purple'
)
plot(
  log10((sumCov %>% filter(rep == 'run2' & swga == 'Sof' & ext == 'Chelex'))$depth),
  log10((sumCov %>% filter(rep == 'run1' & swga == 'Sof' & ext == 'Chelex'))$depth),
  col = 'darkgreen'
)



# comparison of percentiles
sumCov %>%
  filter(chrom == 3 ) %>%
  mutate(depBool = ifelse(depth == 0, 0, 1)) %>%
  group_by(swga, ext, den, depBool) %>%
  tally() %>%
  spread(key = depBool, val = n) %>%
  mutate(pcBad = `0` / `1`)


# violin summaries
sumCov$depth[sumCov$depth == 0] <- 0.01
violin <- ggplot(data = sumCov, aes(x = ext, y = log10(depth), fill = interaction(ext,swga))) +
  geom_violin(draw_quantiles = c(0.25,0.50,0.75)) +
  facet_grid(chrom~swga) +
  theme_classic()
ggsave("violin.pdf", violin, width = 12, height = 24)





#####################
# Regional Coverage #
#####################

proteinCoverage <- ggplot(protCov %>%
    filter(den == 100),
    aes(x = pos, y = (depth), fill = interaction(swga, rep))) +
  # geom_point(shape = 21, alpha = 0.3, size = 0.5) +
  geom_line(aes(colour = interaction(swga,rep))) +
  facet_grid(ext~prot, scale = 'free_x')
ggsave("proteinCoverage.tiff", proteinCoverage, width = 24, height = 12, device = 'tiff')

quantile(protCov$depth, probs = seq(1, 0.1, -0.1))


protpc <- protCov %>%
  group_by(prot, rep, ext, den, swga) %>%
  do( tidy(t(quantile(.$depth, probs = seq(0.1, 1, 0.1))))) %>%
  gather('percentile', 'coverage', X10.:X100.) %>%
  mutate(percentile = gsub('X','',percentile)) %>%
  mutate(percentile = as.numeric(gsub('.$','',percentile)))

protPC_plot <- ggplot(protpc , aes(x = as.factor(percentile), y = coverage, colour = interaction(ext, swga))) +
  geom_point() +
  geom_line(alpha = 0.5, aes(group = interaction(ext, swga, rep))) +
  facet_grid(den~prot) +
  theme_classic()
ggsave('protpc.pdf', protPC_plot, width = 24, height = 12)


########################
# Double SWGA Coverage #
########################
dbl <- ggplot(data = doubleCov %>%
    filter(den == 1000 & chrom == 1),
    aes(x = pos, y = log10(dep), fill = interaction(swga1, swga2))) +
  geom_point(shape = 21, alpha = 0.3) +
  facet_wrap(~ext)
ggsave("doubleSWGA.tiff", dbl, width = 24, height = 12, device = 'tiff')

doubleCov %>% head()
doubleCov <- doubleCov %>%
  group_by(swga1, swga2, ext, den) %>%
  mutate(index = row_number())

## full genomes
dblmassive <- ggplot(doubleCov %>% filter(den == 1000), aes(x = index, y = log10(dep), fill = interaction(swga1,swga2))) +
  geom_point(shape = 21, aes(alpha = 0.3)) +
  # geom_line(aes(colour = rep)) +
  facet_grid(~ext, scale='free_x')

ggsave("dbl_massiveGenome.png", dblmassive, width = 30, height = 12, device = 'png')


doubleCov %>% head()
cov %>% head()

# violin plots double swga
doubleCov$dep[doubleCov$dep == 0] <- 0.01
dblViolin <- ggplot(data = doubleCov, aes(x = interaction(swga1, swga2), y = log10(dep))) +
  geom_violin(aes(colour = interaction(swga1,swga2)), draw_quantiles = c(0.25,0.50,0.75)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(colour = FALSE) +
  facet_wrap(~ext)
ggsave('doubleviolin.pdf', dblViolin)


##############################
# Filtered SWGA Combinations #
##############################

# add total mapping to combination meta data
comboMeta <- comboMeta %>%
  mutate(totalMapping = perfectMapping * total)

# normalize depth to total mapping reads
normFiltDouble <- filteredDouble %>%
  left_join(comboMeta) %>%
  mutate(normalizedDepth = depth / totalMapping)


# pick a chrom and compare double of same swga method
ggplot(normFiltDouble %>% filter(chrom == 3) %>% filter(nchar(s1) == nchar(s2)),
    aes(x = pos, y = log10(normalizedDepth), colour = interaction(s1, s2))) +
  geom_point(shape = 21, alpha = 0.5) +
  facet_wrap(ext~density)

nfd_dropouts <- normFiltDouble %>%
  mutate(dropout = ifelse(depth == 0, 1, 0)) %>%
  group_by(chrom, ext, s1, s2, density) %>%
  summarise(
    numDropout = sum(dropout),
    numPositions = n()
  ) %>%
  mutate(r_dropout = numDropout / numPositions)


# dropout comparison
ggplot(nfd_dropouts %>% filter(nchar(s1) == nchar(s2)) ,
    aes(x = ext, y = r_dropout, colour = interaction(s1,s2))) +
  geom_boxplot() +
  facet_wrap(~density)


# set dropouts to 10e-9 to keep dropout visible
normFiltDouble$normalizedDepth[normFiltDouble$normalizedDepth == 0] <- 0.000000001

# plot violins
ggplot(normFiltDouble, aes(x = interaction(s1, s2), y = log10(normalizedDepth))) +
  geom_violin(aes(colour = interaction(s1,s2)), draw_quantiles = c(0.25,0.5,0.75)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(colour = FALSE) +
  facet_wrap(density~ext)

##########################################
# Downsampled Filtered SWGA Combinations #
##########################################

dfd <- ds_filt_dbl %>%
  left_join(ds_comboMeta) %>%
  mutate(normDepth = depth / (total * perfectMapping))

dfd_dropouts <- dfd %>%
  mutate(dropout = ifelse(depth == 0, 1, 0)) %>%
  group_by(chrom, ext, s1, s2, density) %>%
  summarise(
    numDropout = sum(dropout),
    numPositions = n()
  ) %>%
  mutate(r_dropout = numDropout / numPositions)


masked_dfd <- mask %>%
  mutate(chrom = as.numeric(chrom)) %>%
  filter(regionType == 'Core') %>%
  right_join(dfd) %>%
  filter(pos <= end & pos >= start)

masked_dfd_dropouts <- masked_dfd %>%
  mutate(dropout = ifelse(depth == 0, 1, 0)) %>%
  group_by(chrom, ext, s1, s2, density) %>%
  summarise(
    numDropout = sum(dropout),
    numPositions = n()
  ) %>%
  mutate(r_dropout = numDropout / numPositions)

# downsampled dropouts
ggplot(dfd_dropouts %>% filter(nchar(s1) == nchar(s2)) ,
    aes(x = ext, y = r_dropout, colour = interaction(s1,s2))) +
  geom_boxplot() +
  facet_wrap(~density)

# masked downsampled dropouts
ggplot(masked_dfd_dropouts %>% filter(nchar(s1) == nchar(s2)) ,
       aes(x = ext, y = r_dropout, colour = interaction(s1,s2))) +
  geom_boxplot() +
  facet_wrap(~density)


# violin of normalized downsampled depth
dfd$normDepth[dfd$normDepth == 0] <- (1 * (10**-10))

ggplot(dfd, aes(x = interaction(s1,s2), y = log10(normDepth))) +
  geom_violin(aes(colour = interaction(s1,s2)), draw_quantiles = c(0.25,0.5,0.75)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(colour = FALSE) +
  facet_wrap(density~ext)


# violin of normalized downsampled depth post-masking
masked_dfd$normDepth[masked_dfd$normDepth == 0] <- (1 * (10**-10))

ggplot(masked_dfd, aes(x = interaction(s1,s2), y = log10(normDepth))) +
  geom_violin(aes(colour = interaction(s1,s2)), draw_quantiles = c(0.25,0.5,0.75)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(colour = FALSE) +
  facet_wrap(density~ext)

# histogram of conditions with masked downsampled combinations
ggplot(masked_dfd, aes(x = log10(normDepth))) +
  geom_histogram() +
  geom_vline(xintercept = log10(10**-7)) +
  facet_wrap(s1~s2~ext~density)

###########################################################
# Comparison of Downsampled swga merge and original swgas #
###########################################################

t <- cov %>%
  mutate(
    s1 = paste(swga, rep, sep = ''),
    s2 = paste(swga, rep, sep = '')
    ) %>%
  select(-swga, -rep) %>%
  bind_rows(dfd)


masked_t <- mask %>%
  mutate(chrom = as.numeric(chrom)) %>%
  filter(regionType == 'Core') %>%
  right_join(t) %>%
  filter(pos <= end & pos >= start)

kino <- masked_t %>%
  mutate(dropout = ifelse(depth == 0, 1, 0)) %>%
  group_by(chrom, ext, density, s1, s2) %>%
  summarise(r_dropout = sum(dropout) / n()) %>%
  mutate(combination = case_when(
    s1 == s2 ~ s1,
    TRUE ~ paste(s1, s2, sep = '.'))
  ) %>%
  mutate(colorVar = ifelse(grepl('sanger', combination), 'sanger', 'sof'))

ggplot(kino, aes(x = combination, y = r_dropout, colour = colorVar)) +
  geom_boxplot() +
  facet_wrap(ext~density, scales = 'free')
