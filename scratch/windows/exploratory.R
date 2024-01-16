library(tidyverse)

#### Number of singletons per bin
# Subtype: A>G
# Population: EUR

df <- read_tsv("scratch/windows/data/singletons/bin_counts.txt", col_names = c("bin_id", "singletons"))

df %>%
  summary()

df %>%
  ggplot(aes(x = singletons)) +
  geom_density()
sum(df$singletons > 10) / length(df$singletons)
