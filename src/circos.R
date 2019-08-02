library(tidyverse)
library(circlize)


bed = read_tsv("../data/complete_set/collapse.tab")

total <- bed$count %>% sum()
bed$value <- bed$count / total
bed <- bed %>%
  filter(
    (chr != 'Pf3D7_API_v3') & (chr != 'Pf_M76611')
  ) %>%
  mutate(
    chr = gsub('_v3', '', chr),
    chr = gsub('Pf3D7_', 'chr_', chr)
  )

pdf("../plots/genomic_circos.pdf")
circos.genomicInitialize(bed)
circos.genomicTrack(bed,
    panel.fun = function(region, value, ...) {
        circos.genomicLines(region, value)
      }
)
dev.off()
