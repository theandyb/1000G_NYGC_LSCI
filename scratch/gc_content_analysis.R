library(tidyverse)
library(vroom)
library(ggbeeswarm)

# Let's first look at the GC content in 10-mer windows
# We'll do our first comparison using the ALL population
# and compare the gc content of singletons to their matched controls

c_dir <- "/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/controls/ALL/gc_content/"
s_dir <- "/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/singletons/ALL/gc_content/"

df_s <- read_tsv(paste0(s_dir, "AT_GC5_cont.bed"),
                 skip = 1,
                 col_names = c("chr", "begin", "end", "pct_at", "pct_gc", "n_A", "n_C", "n_G", "n_T", "n_N", "n_O", "len")) %>%
  select(chr, begin, end, starts_with("n_")) %>%
  mutate(pct_gc = (n_C + n_G)/(n_A + n_C + n_G + n_T))

df_c <- vroom(paste0(c_dir, "AT_GC5_cont.bed"),
              delim="\t",
                 skip = 1,
                 col_names = c("chr", "begin", "end", "pct_at", "pct_gc", "n_A", "n_C", "n_G", "n_T", "n_N", "n_O", "len")) %>%
  select(chr, begin, end, starts_with("n_")) %>%
  mutate(pct_gc = (n_C + n_G)/(n_A + n_C + n_G + n_T))

df_c_meta <- vroom("/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/controls/ALL/AT_GC.csv",
                      col_names = c("chr", "s_pos", "ref1", "ref2", "window", "dist", "c_pos", "motif"),
                   delim = ",") %>%
  select(chr, s_pos, window, dist)
df_s_meta <- vroom("/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/singletons/ALL/AT_GC.txt",
                      col_names = c("chr", "pos", "ref"), delim = "\t") %>%
  select(-ref)

df_c_meta$c_gc <- df_c$pct_gc
df_s_meta$s_gc <- df_s$pct_gc

rm(df_s, df_c)
gc()

df_c_meta <- df_c_meta %>%
  left_join(df_s_meta, by = c("chr", "s_pos" = "pos"))

df_c_meta %>%
  filter(chr == "chr22") %>%
  ggplot(aes(x = s_gc, y = c_gc)) +
  stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    legend.position='none'
  ) +
  ggtitle("Singleton GC Content vs Matched Control", "A → G, 500-mer Window, Chr 22") +
  xlab("Singleton GC Content") +
  ylab("Control GC Content")

df_c_meta %>%
  filter(chr == "chr22") %>%
  select(c_gc, s_gc) %>%
  rename(Singletons = s_gc, Controls = c_gc) %>%
  pivot_longer(everything(), names_to = "status", values_to = "gc") %>%
  ggplot(aes(x = status, y = gc)) +
  geom_boxplot() +
  xlab("Status") +
  ylab("GC Content") +
  ggtitle("Singleton GC Content vs Matched Control", "A → G, 10-0mer Window, chr22")

c12 <- df_c_meta %>%
  filter(chr == "chr12") %>%
  pull(s_gc)
c22 <- df_c_meta %>%
  filter(chr == "chr22") %>%
  pull(s_gc)
t.test(c12, c22)

t.test(df_c_meta$s_gc, df_c_meta$c_gc)

for(i in 1:22){
  chrom_str = paste0("chr",i)
  mod_obj <- df_c_meta %>%
    filter(chr == chrom_str) %>%
    {t.test(.$s_gc, .$c_gc)}
  print(paste0(i, " - ", mod_obj$statistic, " - p = ", mod_obj$p.value))
  rm(mod_obj)
}

mod_obj <- df_c_meta %>%
  filter(chr == "chr22") %>%
  {t.test(.$s_gc, .$c_gc)}


mod_obj <- df_c_meta %>%
  filter(chr == "chr1") %>%
  lm(c_gc ~ s_gc + dist - 1, data = .)

for(i in 1:22){
  chrom_str = paste0("chr",i)
  mod_obj <- df_c_meta %>%
    filter(chr == chrom_str) %>%
    lm(c_gc ~ s_gc + dist, data = .) %>%
    summary()
  print(paste0(i, " - ", mod_obj$r.squared ))
  rm(mod_obj)
}

df_c_meta %>%
  select(c_gc, s_gc) %>%
  rename(Singletons = s_gc, Controls = c_gc) %>%
  mutate(diff_gc = Controls - Singletons) %>%
  ggplot(aes(x = diff_gc)) + geom_density()


### !!! From here on out we will only be looking at chr2 !!!
chrom_str = "chr2"
df_s_meta2 <- df_s_meta %>%
  filter(chr == chrom_str)

df_c_meta2 <- df_c_meta %>%
  filter(chr == chrom_str)

df_s_meta2$c_gc_min <- df_c_meta2 %>%
  select(s_pos, dist, c_gc, s_gc) %>%
  group_by(s_pos) %>%
  slice(which.min(dist)) %>%
  pull(c_gc)

df_s_meta2$c_gc_max <- df_c_meta2 %>%
  select(s_pos, dist, c_gc, s_gc) %>%
  group_by(s_pos) %>%
  slice(which.max(dist)) %>%
  pull(c_gc)

df_s_meta2$max_d <- df_c_meta2 %>%
  select(s_pos, dist, c_gc, s_gc) %>%
  group_by(s_pos) %>%
  slice(which.max(dist)) %>%
  pull(dist)

df_s_meta2$min_d <- df_c_meta2 %>%
  select(s_pos, dist, c_gc, s_gc) %>%
  group_by(s_pos) %>%
  slice(which.min(dist)) %>%
  pull(dist)

df_s_meta2 %>%
  pivot_longer(contains("gc"), names_to = "Category", values_to = "GC") %>%
  mutate(Category = fct_recode(factor(Category), Singleton = "s_gc", `Near Control` = "c_gc_min", `Far Control` = "c_gc_max") ) %>%
  ggplot(aes(x = Category, y = GC)) +
  geom_boxplot() +
  ggtitle(chrom_str)

df_s_meta2 %>%
  pivot_longer(contains("gc"), names_to = "Category", values_to = "GC") %>%
  mutate(Category = fct_recode(factor(Category), Singleton = "s_gc", `Near Control` = "c_gc_min", `Far Control` = "c_gc_max") ) %>%
  ggplot(aes(x = GC, fill = Category)) +
  geom_histogram(position="dodge") +
  ggtitle(chrom_str)

t.test(df_s_meta2$s_gc, df_s_meta2$c_gc_min)
t.test(df_s_meta2$s_gc, df_s_meta2$c_gc_max)


### AFR - EUR Comp
eur_dir <- "/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/singletons/EUR/gc_content/"
afr_dir <- "/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/singletons/AFR/gc_content/"

df_eur <- read_tsv(paste0(eur_dir, "AT_GC5_cont.bed"),
                   skip = 1,
                   col_names = c("chr", "begin", "end", "pct_at", "pct_gc", "n_A", "n_C", "n_G", "n_T", "n_N", "n_O", "len")) %>%
  select(chr, begin, end, starts_with("n_")) %>%
  mutate(pct_gc = (n_C + n_G)/(n_A + n_C + n_G + n_T))

df_afr <- read_tsv(paste0(afr_dir, "AT_GC5_cont.bed"),
                   skip = 1,
                   col_names = c("chr", "begin", "end", "pct_at", "pct_gc", "n_A", "n_C", "n_G", "n_T", "n_N", "n_O", "len")) %>%
  select(chr, begin, end, starts_with("n_")) %>%
  mutate(pct_gc = (n_C + n_G)/(n_A + n_C + n_G + n_T))

df_eur_meta <- vroom("/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/singletons/EUR/AT_GC.txt",
                     col_names = c("chr", "pos", "ref"), delim = "\t") %>%
  select(-ref)
df_afr_meta <- vroom("/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/singletons/AFR/AT_GC.txt",
                     col_names = c("chr", "pos", "ref"), delim = "\t") %>%
  select(-ref)

df_eur_meta$s_gc <- df_eur$pct_gc
df_afr_meta$s_gc <- df_afr$pct_gc

rm(df_eur, df_afr)
gc()

df_afr_meta$pop <- "AFR"
df_eur_meta$pop <- "EUR"

df <- bind_rows(df_afr_meta, df_eur_meta)

t.test(df_afr_meta$s_gc, df_eur_meta$s_gc)

df %>%
  filter(chr == "chr22") %>%
  ggplot(aes(x = pop, y = s_gc)) +
  geom_boxplot()
