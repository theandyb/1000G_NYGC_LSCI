library(tidyverse)
library(patchwork)

neutral <- c(0.59/2, 0.41/2, 0.41/2, 0.59/2)
N1 <- 14411404
N2 <- 435216

sim_dev <- function(N, p1, p2){
  df <- data.frame(n = c(rmultinom(1, N, p1) %>% as.vector(),
                   rmultinom(1, N*5, p2) %>% as.vector()),
             status = rep(c("singletons", "controls"), each = 4),
             nuc = rep(c("A", "C", "G", "T"), 2) )
  mod_obj <- glm(n ~ nuc + status, family = poisson, data = df)
  return(deviance(mod_obj))
}

sim_res <- data.frame(dev = replicate(1000, sim_dev(N1, neutral, neutral)))
p1 <- sim_res %>%
  ggplot(aes(x = dev)) + geom_density() +
  ggtitle("Deviance Statistics under No Difference in Distribution",
          paste0("Number of singletons: ", N))


sim_res <- data.frame(dev = replicate(1000, sim_dev(N2, neutral, neutral)))
p2 <- sim_res %>%
  ggplot(aes(x = dev)) + geom_density() +
  ggtitle("Deviance Statistics under No Difference in Distribution",
          paste0("Number of singletons: ", N))

p1 + p2

## Add some perturbation to singleton density
non_neutral <- neutral + c(0.002, -0.0005, -0.0005, -0.001)
l1_n <- sum(abs(non_neutral - neutral))
sim_res <- data.frame(dev = replicate(1000, sim_dev(N1, non_neutral, neutral)))
p1 <- sim_res %>%
  mutate(re = dev / (12*N1)) %>%
  ggplot(aes(x = dev)) + geom_density() +
  ggtitle("Deviance Statistics under Perturb.",
          paste0("Number of singletons: ", N1, "; dist: ", l1_n))

non_neutral <- neutral + c(0.002, -0.0005, -0.0005, -0.001)
l1_n <- sum(abs(non_neutral - neutral))
sim_res <- data.frame(dev = replicate(1000, sim_dev(N2, non_neutral, neutral)))
p2 <- sim_res %>%
  mutate(re = dev / (12*N2)) %>%
  ggplot(aes(x = dev)) + geom_density() +
  ggtitle("Deviance Statistics under Perturb.",
          paste0("Number of singletons: ", N2, "; dist: ", l1_n))

p1 + p2
