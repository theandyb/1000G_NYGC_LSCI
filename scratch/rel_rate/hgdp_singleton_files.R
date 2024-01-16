library(tidyverse)

# id files
id_dir <- "/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/data/hgdp/"
pops <- c("Africa", "America", "Central_South_Asian", "East_Asia", "Europe", "Middle_East", "Oceania")

af_id <- read_csv(paste0(id_dir, pops[1], ".csv"))
am_id <- read_csv(paste0(id_dir, pops[2], ".csv"))
ca_id <- read_csv(paste0(id_dir, pops[3], ".csv"))
ea_id <- read_csv(paste0(id_dir, pops[4], ".csv"))
eu_id <- read_csv(paste0(id_dir, pops[5], ".csv"))
me_id <- read_csv(paste0(id_dir, pops[6], ".csv"))
oc_id <- read_csv(paste0(id_dir, pops[7], ".csv"))


# generate population-specific lists of singletons
sing_dir <- "/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/data/hgdp/annotated/"

df <- read_csv(paste0(sing_dir, "chr1.csv"), col_names = c("chr", "pos", "ref", "alt", "cat", "id", "motif")) %>%
  filter(cat == "AT_GC") %>% select(chr, pos, ref, id, motif)

af_df <- df %>% filter(id %in% af_id$id) %>% select(-id)
am_df <- df %>% filter(id %in% am_id$id) %>% select(-id)
ca_df <- df %>% filter(id %in% ca_id$id) %>% select(-id)
ea_df <- df %>% filter(id %in% ea_id$id) %>% select(-id)
eu_df <- df %>% filter(id %in% eu_id$id) %>% select(-id)
me_df <- df %>% filter(id %in% me_id$id) %>% select(-id)
oc_df <- df %>% filter(id %in% oc_id$id) %>% select(-id)

for(i in 2:22){
  df <- read_csv(paste0(sing_dir, "chr", i,".csv"), col_names = c("chr", "pos", "ref", "alt", "cat", "id", "motif")) %>%
    filter(cat == "AT_GC") %>% select(chr, pos, ref, id, motif)
  af_df <- af_df %>% bind_rows({df %>% filter(id %in% af_id$id) %>% select(-id)})
  am_df <- am_df %>% bind_rows({df %>% filter(id %in% am_id$id) %>% select(-id)})
  ca_df <- ca_df %>% bind_rows({df %>% filter(id %in% ca_id$id) %>% select(-id)})
  ea_df <- ea_df %>% bind_rows({df %>% filter(id %in% ea_id$id) %>% select(-id)})
  eu_df <- eu_df %>% bind_rows({df %>% filter(id %in% eu_id$id) %>% select(-id)})
  me_df <- me_df %>% bind_rows({df %>% filter(id %in% me_id$id) %>% select(-id)})
  oc_df <- oc_df %>% bind_rows({df %>% filter(id %in% oc_id$id) %>% select(-id)})
}

out_dir <- "/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/scratch/rel_rate/hgdp_singletons/"

write_csv(af_df, paste0(out_dir, pops[1], "_A_G.csv"))
write_csv(am_df, paste0(out_dir, pops[2], "_A_G.csv"))
write_csv(ca_df, paste0(out_dir, pops[3], "_A_G.csv"))
write_csv(ea_df, paste0(out_dir, pops[4], "_A_G.csv"))
write_csv(eu_df, paste0(out_dir, pops[5], "_A_G.csv"))
write_csv(me_df, paste0(out_dir, pops[6], "_A_G.csv"))
write_csv(oc_df, paste0(out_dir, pops[7], "_A_G.csv"))
