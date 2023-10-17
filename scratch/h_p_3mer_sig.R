library(tidyverse)

res_dir <- "/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/singletons/"

load_res <- function(pop, res_dir){
  f_name <- paste0(res_dir, pop, "/motifs/3mer/all.txt")
  df <- read_tsv(f_name) %>%
    arrange(subtype, motif) %>%
    filter(str_detect(motif, "[ACGT]{3}"))
  return(df)
}

# load all the counts
df_afr <- load_res("AFR", res_dir) %>%
  rename(AFR = count)
df_amr <- load_res("AMR", res_dir) %>%
  rename(AMR = count)
df_eas <- load_res("EAS", res_dir) %>%
  rename(EAS = count)
df_eur <- load_res("EUR", res_dir) %>%
  rename(EUR = count)
df_sas <- load_res("SAS", res_dir) %>%
  rename(SAS = count)

df <- full_join(df_afr, df_amr) %>%
  full_join(df_eas) %>%
  full_join(df_eur) %>%
  full_join(df_sas)

df <- df %>%
  mutate(nAFR = sum(AFR) - AFR,
         nAMR = sum(AMR) - AMR,
         nEAS = sum(EAS) - EAS,
         nEUR = sum(EUR) - EUR,
         nSAS = sum(SAS) - SAS)

res <- data.frame(motif = df$motif,
                  subtype = df$subtype,
                  EUR = rep(0, 96),
                  AFR = rep(0,96),
                  stat = rep(0, 96),
                  pval = rep(0, 96))

for(i in 1:96){
  d_mat <- matrix(c(df$EUR[i], df$nEUR[i],
                    df$AFR[i], df$nAFR[i]),
                  byrow = TRUE, nrow = 2)
  obj <- chisq.test(d_mat)
  res$stat[i] <- obj$statistic
  res$pval[i] <- obj$p.value
  res$EUR[i] <- df$EUR[i]
  res$AFR[i] <- df$AFR[i]
}

res <- res %>% arrange(desc(stat))
res$tAFR <- sum(res$AFR)
res$tEUR <- sum(res$EUR)
res$pAFR <- res$tAFR
res$pEUR <- res$tEUR

for(i in 2:95){
  res$pEUR[i] <- sum(res$EUR[(i+1):96])
  res$pAFR[i] <- sum(res$AFR[(i+1):96])
}

res$pEUR[96] <- res$EUR[95]
res$pAFR[96] <- res$AFR[95]

res$nstat <- res$stat
res$pnew <- res$pval

for(i in 2:96){
  d_mat <- matrix(c(res$EUR[i], res$pEUR[i],
                    res$AFR[i], res$pAFR[i]),
                  byrow = TRUE, nrow = 2)
  obj <- chisq.test(d_mat)
  res$nstat[i] <- obj$statistic
  res$pnew[i] <- obj$p.value
}

sum(res$pnew < 0.05 / 96)

sum(res$pnew < 1e-5)
sum(res$pval < 0.05)

res$propEUR <- res$EUR / sum(res$EUR)
res$propAFR <- res$AFR / sum(res$AFR)

res <- res %>%
  select(motif, subtype, EUR, AFR, nstat, pnew, propEUR, propAFR) %>%
  mutate(full_cat = paste0(motif, ",", subtype))

res$enrich <- res$propEUR / res$propAFR
max(res$enrich)
min(res$enrich)
res$type <- "EUR_AFR"

res[1:10,] %>%
  ggplot(aes(y =full_cat, x = type, fill = enrich)) +
  geom_tile() +
  scale_fill_gradient2(midpoint = 1, low = "green", mid = "white", high = "red") +
  scale_y_discrete(limits = res$full_cat[10:1]) +
  theme_classic()
