library(tidyverse)

# start with 7-mers
load_singletons <- function(pop, subtype){
  base_dir <- "/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/singletons/"
  f_name <- paste0(base_dir, pop, "/motifs/", subtype, ".csv")
  awk_tmpl <- "awk -F, '{count[substr($4,8,7)]++}END{for(key in count)print(key\"\\t\"count[key])}' "
  awk_cmd <- paste0(awk_tmpl, f_name)
  df <- vroom::vroom(pipe(awk_cmd), col_names = c("motif", "singletons"), show_col_types = FALSE) %>%
    separate(motif, into = c("p1", "p2", "p3", "p4", "p5", "p6", "p7"), sep=1:7)
  return(df)
}

load_controls <- function(pop, subtype){
  base_dir <- "/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/controls/"
  f_name <- paste0(base_dir, pop, "/motifs/", subtype, ".csv")
  awk_tmpl <- "awk -F, '{count[substr($4,8,7)]++}END{for(key in count)print(key\"\\t\"count[key])}' "
  awk_cmd <- paste0(awk_tmpl, f_name)
  df <- vroom::vroom(pipe(awk_cmd), col_names = c("motif", "controls"), show_col_types = FALSE) %>%
    separate(motif, into = c("p1", "p2", "p3", "p4", "p5", "p6", "p7"), sep=1:7)
  return(df)
}

df_s <- load_singletons("ALL", "AT_GC")
df_c <- load_controls("ALL", "AT_GC")

df <- full_join(df_s, df_c, by = c("p1", "p2", "p3", "p4", "p5", "p6", "p7")) %>%
  pivot_longer(!starts_with("p"), names_to = "status", values_to = "n") %>%
  filter(p1 %in% c("A","C","G","T"),
         p2 %in% c("A","C","G","T"),
         p3 %in% c("A","C","G","T"),
         p5 %in% c("A","C","G","T"),
         p6 %in% c("A","C","G","T"),
         p7 %in% c("A","C","G","T"))

df <- as.data.frame(unclass(df),stringsAsFactors=TRUE)

#mod_0 <- glm(n ~ (p1 + p2 + p3 + p5 + p6 + p7)^2 + p2*p3*p5*p6 + status, data = df, family = poisson())

ptm <- proc.time()
mod_0 <- glm(n ~ (p1 * p2 * p3 * p5 * p6 * p7) + status, data = df, family = poisson())
proc.time() - ptm


summary(mod_0)
logLik(mod_0)

ptm <- proc.time()
mod_1_1 <- update(mod_0, . ~ . + p5 * status)
proc.time() - ptm

summary(mod_1_1)
logLik(mod_1_1)

mod_1_2 <- update(mod_0, . ~ . + p2 * status)
summary(mod_1_2)
logLik(mod_1_2)

mod_1_3 <- update(mod_0, . ~ . + p3 * status)
summary(mod_1_3)
logLik(mod_1_3)

################################################
mod_2_1 <- update(mod_1_1, . ~ . + p2 * status)
summary(mod_2_1)
logLik(mod_2_1)

mod_2_2 <- update(mod_1_1, . ~ . + p3 * status)
summary(mod_2_2)
logLik(mod_2_2)

mod_2_3 <- update(mod_1_1, . ~ . + p5 * p6 * status - p6*status)
summary(mod_2_3)
logLik(mod_2_3)

mod_2_4 <- update(mod_1_1, . ~ . + p7 * status)
summary(mod_2_4)
logLik(mod_2_4)

#################################################
mod_3_1 <- update(mod_2_1, . ~ . + p3*status)
summary(mod_3_1)
logLik(mod_3_1)


##############################################################################
##
## glmnet idea (reduction of base model)
##
##############################################################################
library(glmnet)

X <- model.matrix(~ status + (p1 * p2 * p3 * p5 * p6 * p7) - 1, df)
y <- df$n

mod_glmnet <- glmnet(X, y, family = "poisson")
plot(mod_glmnet)
coef(mod_glmnet, s = 100) %>% head

# cross validation
cvfit <- cv.glmnet(X, y, family = "poisson")
plot(cvfit)

# lambda value with best fit in regards to minimum mean cross-validated error
l_min <- cvfit$lambda.min
coef_min_l <- coef(cvfit, s = "lambda.min") %>% as.matrix() %>% .[,1]
coef_non_zero <- coef_min_l[coef_min_l > 0]
##############################################################################
##
## Stepwise Approach
##
##############################################################################

library(MASS)
select <- dplyr::select

base_mod <- glm(n ~ (p1 + p2 + p3 + p5 + p6 + p7)^3 + status, data = df, family = poisson())
mod.stp <- stepAIC(base_mod,
                   scope = list(upper = ~(p1+p2+p3+p5+p6+p7+status)^7, lower = ~1),
                   trace = TRUE)

##############################################################################
##
## Code for sequentially fitting all models
##
##############################################################################

base_mod <-glm(n ~ p1 + p2 + p3 + p5 + p6 + p7 + status, data = df, family = poisson())
