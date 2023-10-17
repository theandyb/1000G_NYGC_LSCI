library(tidyverse)

nuc_comp <- function(x){
  x <- toupper(x)
  if(x == "A") return("T")
  else if(x == "C") return("G")
  else if(x == "G") return("C")
  else if(x == "T") return("A")
  else return(x)
}

recode_st <- function(x){
  if(str_sub(x,1,1) %in% c("A", "C")){
    return(x)
  } else{
    p1 <- str_sub(x,1,1)
    p2 <- str_sub(x,3,3)
    q1 <- nuc_comp(p1)
    q2 <- nuc_comp(p2)
    return(paste0(q1, ">", q2))
  }
}

subj_info <- read_tsv("data/1000G_samples.tsv")

singleton_df <- read_csv("data/chr22_sing.ann.csv", col_names = c("chr", "pos", "motif", "st", "alt", "Sample", "ref", "motif2"))

singleton_df <-
  left_join(singleton_df, {subj_info %>% select(Sample, SuperpopCode)}) %>%
  filter(alt %in% c("A", "C", "G", "T")) %>%
  rowwise() %>%
  filter(alt != ref)

singleton_df <- singleton_df %>%
  rowwise() %>%
  mutate(st2 = recode_st(st))

singleton_df$st2 %>%
  table()

count_df <- singleton_df %>%
  select(Sample, SuperpopCode, motif2, st2) %>%
  group_by(SuperpopCode, Sample, motif2, st2) %>%
  summarize(count = n()) %>%
  ungroup() %>%
  group_by(Sample) %>%
  mutate(prop_c = count / sum(count)) %>%
  ungroup()

count_df %>%
  filter(Sample == "HG01879") %>%
  arrange(motif2) %>%
  mutate(st = paste0(motif2, ",", st2)) %>%
  ggplot(aes(x = st, y = prop_c)) +
  geom_col()

enrich_df <- count_df %>%
  select(SuperpopCode, motif2, st2, prop_c) %>%
  group_by(SuperpopCode, motif2, st2) %>%
  summarize(mean_prop = mean(prop_c)) %>%
  rowwise() %>%
  mutate(subtype = paste0(motif2, "," ,st2)) %>%
  ungroup() %>%
  select(-motif2, -st2) %>%
  pivot_wider(names_from = SuperpopCode, values_from = mean_prop)

### AFR - EUR enrichments
enrich_df %>%
  mutate(en_eur_afr = EUR/AFR)
