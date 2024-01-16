library(tidyverse)
library(ggbeeswarm)

# Get relative rates
## genome-wide

df_gw_3mer <- read_tsv("/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/scratch/rel_rate/gw_A.txt") %>%
  rename(gw_n = n) %>%
  mutate(gw_rr = gw_n / sum(gw_n))
df_gw_3mer_aug <- read_tsv("/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/scratch/rel_rate/gw_A_aug.txt") %>%
  rename(gw_n = n) %>%
  mutate(gw_rr = gw_n / sum(gw_n))

## control rates
# pop = ALL
df_c_3mer <- read_tsv("/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/scratch/rel_rate/1kgp_rates/EUR/A_G_control_3mer_count.txt") %>%
  rename(c_n = n) %>%
  filter(c_n > 100) %>%
  mutate(c_rr = c_n / sum(c_n))

df_c_3mer %>%
  mutate(p1 = str_sub(motif, 1, 1),
         p2 = str_sub(motif, 3, 3)) %>%
  ggplot(aes(x = p2, y = p1, fill = c_rr)) +
  geom_tile() +
  scale_fill_distiller(palette = "Blues", direction = 1)

df_c_3mer_aug <- read_tsv("/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/scratch/rel_rate/1kgp_rates/EUR/A_G_control_3mer_count_aug.txt") %>%
  rename(c_n = n) %>%
  filter(c_n > 100) %>%
  mutate(c_rr = c_n / sum(c_n))

df_c_3mer_aug %>%
  mutate(p1 = str_sub(motif, 1, 1),
         p2 = str_sub(motif, 2, 2)) %>%
  ggplot(aes(x = p2, y = p1, fill = c_rr)) +
  geom_tile() +
  scale_fill_distiller(palette = "Blues", direction = 1)

## singleton rates
df_s_3mer <- read_tsv("/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/scratch/rel_rate/1kgp_rates/EUR/A_G_3mer_count.txt") %>%
  rename(s_n = n) %>%
  filter(s_n > 10) %>%
  mutate(s_rr = s_n / sum(s_n))

df_s_3mer_aug <- read_tsv("/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/scratch/rel_rate/1kgp_rates/EUR/A_G_3mer_count_aug.txt") %>%
  rename(s_n = n) %>%
  filter(s_n > 10) %>%
  mutate(s_rr = s_n / sum(s_n))

df_s_3mer %>%
  mutate(type = "central") %>%
  bind_rows({df_s_3mer_aug %>% mutate(type = "nonCen")}) %>%
  ggplot(aes(x = type, y = s_rr)) + geom_boxplot()

## Generate relative rates dfs
df_rr <- full_join(df_gw_3mer, df_s_3mer, by = "motif") %>%
  full_join(df_c_3mer, by = "motif") %>%
  mutate(rr_g = s_rr / gw_rr,
         rr_c = s_rr / c_rr)

df_rr_aug <- full_join(df_gw_3mer_aug, df_s_3mer_aug, by = "motif") %>%
  full_join(df_c_3mer_aug, by = "motif") %>%
  mutate(rr_g = s_rr / gw_rr,
         rr_c = s_rr / c_rr)

df_rr %>%
  mutate(type = "centrol") %>%
  bind_rows({df_rr_aug %>% mutate(type="nonCen")}) %>%
  ggplot(aes(x = type, y = rr_c)) + geom_beeswarm()

# Load singletons and matched controls from HGDP
singleton_dir <- "/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/scratch/rel_rate/hgdp_singletons/"
df_sing_afr <- read_csv(paste0(singleton_dir, "Africa_A_G.csv"))
df_sing_amr <- read_csv(paste0(singleton_dir, "America_A_G.csv"))
df_sing_csa <- read_csv(paste0(singleton_dir, "Central_South_Asian_A_G.csv"))
df_sing_eas <- read_csv(paste0(singleton_dir, "East_Asia_A_G.csv"))
df_sing_eur <- read_csv(paste0(singleton_dir, "Europe_A_G.csv"))
df_sing_mes <- read_csv(paste0(singleton_dir, "Middle_East_A_G.csv"))
df_sing_oce <- read_csv(paste0(singleton_dir, "Oceania_A_G.csv"))

control_dir <- "/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/scratch/rel_rate/hgdp_controls/annotated/"
df_control_afr <- read_csv(paste0(control_dir, "Africa_A_G.csv"), col_names = c("chr", "pos", "ref", "motif"))
df_control_amr <- read_csv(paste0(control_dir, "America_A_G.csv"), col_names = c("chr", "pos", "ref", "motif"))
df_control_csa <- read_csv(paste0(control_dir, "Central_South_Asian_A_G.csv"), col_names = c("chr", "pos", "ref", "motif"))
df_control_eas <- read_csv(paste0(control_dir, "East_Asia_A_G.csv"), col_names = c("chr", "pos", "ref", "motif"))
df_control_eur <- read_csv(paste0(control_dir, "Europe_A_G.csv"), col_names = c("chr", "pos", "ref", "motif"))
df_control_mes <- read_csv(paste0(control_dir, "Middle_East_A_G.csv"), col_names = c("chr", "pos", "ref", "motif"))
df_control_oce <- read_csv(paste0(control_dir, "Oceania_A_G.csv"), col_names = c("chr", "pos", "ref", "motif"))

df_control_afr$status <- "control"
df_control_amr$status <- "control"
df_control_csa$status <- "control"
df_control_eur$status <- "control"
df_control_eas$status <- "control"
df_control_mes$status <- "control"
df_control_oce$status <- "control"

df_sing_afr$status <- "singleton"
df_sing_amr$status <- "singleton"
df_sing_csa$status <- "singleton"
df_sing_eas$status <- "singleton"
df_sing_eur$status <- "singleton"
df_sing_mes$status <- "singleton"
df_sing_oce$status <- "singleton"

df_test_afr <- bind_rows(df_sing_afr, df_control_afr)
df_test_amr <- bind_rows(df_sing_amr, df_control_amr)
df_test_csa <- bind_rows(df_sing_csa, df_control_csa)
df_test_eas <- bind_rows(df_sing_eas, df_control_eas)
df_test_eur <- bind_rows(df_sing_eur, df_control_eur)
df_test_mes <- bind_rows(df_sing_mes, df_control_mes)
df_test_oce <- bind_rows(df_sing_oce, df_control_oce)

# Construct the two types of 3mers in df_test
df_test_afr <- df_test_afr %>%
  mutate(cat1 = str_sub(motif, 2, 4),
         cat2 = str_sub(motif, 1, 3))
df_test_amr <- df_test_amr %>%
  mutate(cat1 = str_sub(motif, 2, 4),
         cat2 = str_sub(motif, 1, 3))
df_test_csa <- df_test_csa %>%
  mutate(cat1 = str_sub(motif, 2, 4),
         cat2 = str_sub(motif, 1, 3))
df_test_eas <- df_test_eas %>%
  mutate(cat1 = str_sub(motif, 2, 4),
         cat2 = str_sub(motif, 1, 3))
df_test_eur <- df_test_eur %>%
  mutate(cat1 = str_sub(motif, 2, 4),
         cat2 = str_sub(motif, 1, 3))
df_test_mes <- df_test_mes %>%
  mutate(cat1 = str_sub(motif, 2, 4),
         cat2 = str_sub(motif, 1, 3))
df_test_oce <- df_test_oce %>%
  mutate(cat1 = str_sub(motif, 2, 4),
         cat2 = str_sub(motif, 1, 3))

# attach relative rates to each singleton/control
df_test_afr <- df_test_afr %>%
  left_join({df_rr %>% select(motif, rr_g, rr_c)}, by = c("cat1" = "motif")) %>%
  rename(std_rr_c = rr_c,
         std_rr_g = rr_g) %>%
  left_join({df_rr_aug %>% select(motif, rr_g, rr_c)}, by = c("cat2" = "motif")) %>%
  rename(aug_rr_c = rr_c,
         aug_rr_g = rr_g) %>%
  remove_missing()

df_test_amr <- df_test_amr %>%
  left_join({df_rr %>% select(motif, rr_g, rr_c)}, by = c("cat1" = "motif")) %>%
  rename(std_rr_c = rr_c,
         std_rr_g = rr_g) %>%
  left_join({df_rr_aug %>% select(motif, rr_g, rr_c)}, by = c("cat2" = "motif")) %>%
  rename(aug_rr_c = rr_c,
         aug_rr_g = rr_g) %>%
  remove_missing()

df_test_csa <- df_test_csa %>%
  left_join({df_rr %>% select(motif, rr_g, rr_c)}, by = c("cat1" = "motif")) %>%
  rename(std_rr_c = rr_c,
         std_rr_g = rr_g) %>%
  left_join({df_rr_aug %>% select(motif, rr_g, rr_c)}, by = c("cat2" = "motif")) %>%
  rename(aug_rr_c = rr_c,
         aug_rr_g = rr_g) %>%
  remove_missing()

df_test_eas <- df_test_eas %>%
  left_join({df_rr %>% select(motif, rr_g, rr_c)}, by = c("cat1" = "motif")) %>%
  rename(std_rr_c = rr_c,
         std_rr_g = rr_g) %>%
  left_join({df_rr_aug %>% select(motif, rr_g, rr_c)}, by = c("cat2" = "motif")) %>%
  rename(aug_rr_c = rr_c,
         aug_rr_g = rr_g) %>%
  remove_missing()

df_test_eur <- df_test_eur %>%
  left_join({df_rr %>% select(motif, rr_g, rr_c)}, by = c("cat1" = "motif")) %>%
  rename(std_rr_c = rr_c,
         std_rr_g = rr_g) %>%
  left_join({df_rr_aug %>% select(motif, rr_g, rr_c)}, by = c("cat2" = "motif")) %>%
  rename(aug_rr_c = rr_c,
         aug_rr_g = rr_g) %>%
  remove_missing()

df_test_mes <- df_test_mes %>%
  left_join({df_rr %>% select(motif, rr_g, rr_c)}, by = c("cat1" = "motif")) %>%
  rename(std_rr_c = rr_c,
         std_rr_g = rr_g) %>%
  left_join({df_rr_aug %>% select(motif, rr_g, rr_c)}, by = c("cat2" = "motif")) %>%
  rename(aug_rr_c = rr_c,
         aug_rr_g = rr_g) %>%
  remove_missing()

df_test_oce <- df_test_oce %>%
  left_join({df_rr %>% select(motif, rr_g, rr_c)}, by = c("cat1" = "motif")) %>%
  rename(std_rr_c = rr_c,
         std_rr_g = rr_g) %>%
  left_join({df_rr_aug %>% select(motif, rr_g, rr_c)}, by = c("cat2" = "motif")) %>%
  rename(aug_rr_c = rr_c,
         aug_rr_g = rr_g) %>%
  remove_missing()

# some figures
df_test_afr %>%
  select(status, std_rr_c, aug_rr_c) %>%
  pivot_longer(!status, names_to = "type", values_to = "measure") %>%
  ggplot(aes(x = type, y = measure, colour = status)) +
  geom_boxplot() +
  labs(color = "Relative Rate Type", x = NULL, y = "Relative Rate") +
  scale_color_discrete(labels = c("Control", "Singleton")) +
  scale_x_discrete(labels = c("(-2, -1)", "Centered")) +
  ggtitle("Africa")

df_test_amr %>%
  select(status, std_rr_c, aug_rr_c) %>%
  pivot_longer(!status, names_to = "type", values_to = "measure") %>%
  ggplot(aes(x = status, y = measure, colour = type)) +
  geom_boxplot() +
  labs(color = "Relative Rate Type", x = NULL, y = "Relative Rate") +
  scale_color_discrete(labels = c("(-2, -1)", "Centered")) +
  ggtitle("America")

df_test_csa %>%
  select(status, std_rr_c, aug_rr_c) %>%
  pivot_longer(!status, names_to = "type", values_to = "measure") %>%
  ggplot(aes(x = status, y = measure, colour = type)) +
  geom_boxplot() +
  labs(color = "Relative Rate Type", x = NULL, y = "Relative Rate") +
  scale_color_discrete(labels = c("(-2, -1)", "Centered")) +
  ggtitle("Central South Asian")

df_test_eas %>%
  select(status, std_rr_c, aug_rr_c) %>%
  pivot_longer(!status, names_to = "type", values_to = "measure") %>%
  ggplot(aes(x = status, y = measure, colour = type)) +
  geom_boxplot() +
  labs(color = "Relative Rate Type", x = NULL, y = "Relative Rate") +
  scale_color_discrete(labels = c("(-2, -1)", "Centered")) +
  ggtitle("East Asia")

df_test_eur %>%
  select(status, std_rr_c, aug_rr_c) %>%
  pivot_longer(!status, names_to = "type", values_to = "measure") %>%
  ggplot(aes(x = type, y = measure, colour = status)) +
  geom_boxplot() +
  labs(color = "Relative Rate Type", x = NULL, y = "Relative Rate") +
  scale_color_discrete(labels = c("Control", "Singleton")) +
  scale_x_discrete(labels = c("(-2, -1)", "Centered")) +
  ggtitle("Europe")

df_test_mes %>%
  select(status, std_rr_c, aug_rr_c) %>%
  pivot_longer(!status, names_to = "type", values_to = "measure") %>%
  ggplot(aes(x = status, y = measure, colour = type)) +
  geom_boxplot() +
  labs(color = "Relative Rate Type", x = NULL, y = "Relative Rate") +
  scale_color_discrete(labels = c("(-2, -1)", "Centered")) +
  ggtitle("Middle East")

df_test_oce %>%
  select(status, std_rr_c, aug_rr_c) %>%
  pivot_longer(!status, names_to = "type", values_to = "measure") %>%
  ggplot(aes(x = status, y = measure, colour = type)) +
  geom_boxplot() +
  labs(color = "Relative Rate Type", x = NULL, y = "Relative Rate") +
  scale_color_discrete(labels = c("(-2, -1)", "Centered")) +
  ggtitle("Oceania")

df_test_afr %>%
  mutate(pop="Africa") %>%
  bind_rows({df_test_eur %>% mutate(pop = "Europe")}) +
  select(pop, status, std_rr_c, aug_rr_c)


df_test_eur %>%
  select(status, std_rr_c, aug_rr_c) %>%
  pivot_longer(!status, names_to = "type", values_to = "measure") %>%
  ggplot(aes(x = measure, fill = status)) +
  geom_histogram(position="dodge") +
  facet_wrap( ~ type)
  # ggplot(aes(x = status, y = measure, colour = type)) +
  # geom_boxplot() +
  # labs(color = "Relative Rate Type", x = NULL, y = "Relative Rate") +
  # scale_color_discrete(labels = c("(-2, -1)", "Centered")) +
  # ggtitle("Europe")

df_test_eur %>%
  select(status, std_rr_c) %>%
  group_by(status, std_rr_c) %>%
  summarize(freq_rr_c = n()) %>%
  group_by(status) %>%
  mutate(pct_rr_c = freq_rr_c / sum(freq_rr_c)) %>%
  ggplot(aes(x = std_rr_c, y = pct_rr_c, color = status)) +
  geom_point() + geom_line()

df_test_eur %>%
  select(status, aug_rr_c) %>%
  group_by(status, aug_rr_c) %>%
  summarize(freq_rr_c = n()) %>%
  group_by(status) %>%
  mutate(pct_rr_c = freq_rr_c / sum(freq_rr_c)) %>%
  ggplot(aes(x = aug_rr_c, y = pct_rr_c, color = status)) +
  geom_point() + geom_line()
# logistic regression models
df_test_sub <- df_test %>%
  slice_sample(prop = 0.01)

mod_aug <- df_test_sub %>%
  mutate(is_sing = status == "singleton") %>%
  glm(is_sing ~ aug_rr_c, data = ., family=binomial)

mod_std <- df_test_sub %>%
  mutate(is_sing = status == "singleton") %>%
  glm(is_sing ~ std_rr_c, data = ., family=binomial)

df_test$pred_aug <- predict(mod_aug, newdata = df_test, type = "response")
df_test$pred_std <- predict(mod_std, newdata = df_test, type = "response")

df_test %>%
  slice_sample(prop = 0.01) %>%
  select(status, pred_aug, pred_std) %>%
  pivot_longer(-status, names_to = "pred_type", values_to = "prob") %>%
  ggplot(aes(x = status, y = prob, color = pred_type)) +
  geom_boxplot() +
  labs(color = "Relative Rate Type", x = NULL, y = "Pr(Singleton)") +
  scale_color_discrete(labels = c("(-2, -1)", "Centered"))

df_test_afr %>%
  slice_sample(prop = 0.1) %>%
  ggplot(aes(x = std_rr_g, y = std_rr_c)) +
  geom_point()
