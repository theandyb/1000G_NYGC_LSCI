# Single-position result figures

library(tidyverse)
#library(sjPlot)


base_dir <- "/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/single_pos/"
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
plot_dir <- "/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/figures/single_pos/"
populations <- c("ALL", "AFR", "AMR", "EAS", "EUR", "SAS")
subtypes <- c("AT_CG", "AT_GC", "AT_TA",
              "GC_AT", "GC_TA", "GC_CG",
              "cpg_GC_AT", "cpg_GC_TA", "cpg_GC_CG")

theme_andy <- function(base_size = 10, base_family = "Helvetica"){
  (theme_bw(base_size = base_size, base_family = base_family) +
     theme(axis.line = element_line(color='black'),
           plot.background = element_blank(),
           panel.grid.minor = element_blank(),
           panel.grid.major = element_blank()))
}

# Position-level Results

#' Load position-level results for a population-subtype pair
#'
#' @param dir Location of the result files
#' @param st The mutation sub type
#' @param pop Which super-population
#' @returns A data frame with columns dev, singletons, controls, and offset
load_st_pop_res <- function(dir, st, pop, max_d = 1000){
  file_loc <- paste0(dir, pop, "/", st, ".csv")
  df <- read_csv(file_loc, show_col_types = FALSE) %>%
    filter(abs(offset) <= max_d)
  return(df)
}

subtype_text <- function(st){
  if(str_starts(st, "cpg")){
    x <- "(CpG) C → "
    x <- paste0(x, str_sub(st,9,9))
  } else if (str_starts(st, "GC")){
    x <- "C → "
    x <- paste0(x, str_sub(st, 5, 5))
  } else {
    x <- "A → "
    x <- paste0(x, str_sub(st, 4, 4))
  }
  return(x)
}

#' Plot the position-level results for a population-subtype pair
#'
#' @param dir Location of the result files
#' @param st The mutation sub type
#' @param pop Which super-population
#' @param re Whether to convert the deviance to relative entropy value
#' @returns ggplot object with the figure
plot_st_pop_res <- function(dir, st, pop, re = TRUE, max_d = 1000){
  df <- load_st_pop_res(dir, st, pop, max_d)
  st_text <- subtype_text(st)
  if(re){
    df <- df %>% mutate(stat_val = dev / (2*(singletons + controls)))
    stat_type <- "Relative Entropy"
  } else {
    df <- df %>% rename(stat_val = dev)
    stat_type <- "Deviance"
  }
  g_obj <- df %>%
    ggplot(aes(x = offset, y = stat_val)) +
    geom_point() +
    geom_line() +
    labs(x = "Relative Position",
         y = stat_type) +
    ggtitle("Influence by Position", paste0("Population: ", pop, "; Sub type: ", st_text))
  return(g_obj)
}

pop <- "ALL"
plot_st_pop_res(base_dir, "AT_CG", pop, max_d = 10)
plot_st_pop_res(base_dir, "AT_GC", pop, max_d = 10)
plot_st_pop_res(base_dir, "AT_TA", pop, max_d = 10)
plot_st_pop_res(base_dir, "GC_AT", pop, max_d = 10)
plot_st_pop_res(base_dir, "GC_TA", pop, max_d = 10)
plot_st_pop_res(base_dir, "GC_CG", pop, max_d = 10)
plot_st_pop_res(base_dir, "cpg_GC_AT", pop, max_d = 10)
plot_st_pop_res(base_dir, "cpg_GC_TA", pop, max_d = 10)
plot_st_pop_res(base_dir, "cpg_GC_CG", pop, max_d = 10)

#####
##### GENERATE ALL SINGLE SUBTYPE PLOTS: STATISTIC AT EACH POSITION
#####

for(pop in populations){
  for(st in subtypes){
    f_name <- paste0(plot_dir, "statistics/", pop, "/", st, ".png")
    ggplot_formatted <- plot_st_pop_res(base_dir, st, pop, max_d = 10) +
      theme_classic(base_family = "Helvetica", base_size = 10)
    #save_plot(f_name, width = 8.3, height = (8.3 * 5)/7, dpi = 1200)
    ggsave(filename=f_name, plot=ggplot_formatted, device="png", width=8.3, height=(8.3*5)/7, units="cm", dpi=500)
    #ggsave(f_name, p, dpi = 600, width = 8.3, height = (8.3*5)/7, units = "cm")
  }
}

#' Plot the position-level results for a population-subtype pair, with significance
#'
#' @param dir Location of the result files
#' @param st The mutation sub type
#' @param pop Which super-population
#' @param re Whether to convert the deviance to relative entropy value
#' @returns ggplot object with the figure
plot_st_pop_dev <- function(dir, st, pop, max_d = 1000, crit_val = qchisq(0.95, df = 3), re = TRUE){
  df <- load_st_pop_res(dir, st, pop, max_d)
  st_text <- subtype_text(st)

  df <- df %>%
    rowwise() %>%
    mutate(is_sig = dev > crit_val)

  if(re){
    df <- df %>% mutate(stat_val = dev / (2*(singletons + controls)))
    stat_type <- "Relative Entropy"
  } else {
    df <- df %>% rename(stat_val = dev)
    stat_type <- "Deviance"
  }
  g_obj <- df %>%
    ggplot(aes(x = offset, y = stat_val, color = factor(is_sig))) +
    geom_point() +
    labs(x = "Relative Position",
         y = stat_type) +
    ggtitle("Influence by Position", paste0("Population: ", pop, "; Sub type: ", st_text))
  return(g_obj)
}

###
pop <- "ALL"
plot_st_pop_dev(base_dir, "AT_CG", "ALL", re=FALSE) + scale_y_log10() + theme_andy() + labs(color = "Significant\nChi Sq") + ylab("Log10 Deviance")
plot_st_pop_dev(base_dir, "AT_GC", "ALL", re=FALSE) + scale_y_log10() + theme_andy() + labs(color = "Significant\nChi Sq") + ylab("Log10 Deviance")
plot_st_pop_dev(base_dir, "AT_TA", "ALL", re=FALSE) + scale_y_log10() + theme_andy() + labs(color = "Significant\nChi Sq") + ylab("Log10 Deviance")
plot_st_pop_dev(base_dir, "GC_AT", "ALL", re=FALSE) + scale_y_log10() + theme_andy() + labs(color = "Significant\nChi Sq") + ylab("Log10 Deviance")
plot_st_pop_dev(base_dir, "GC_TA", "ALL", re=FALSE) + scale_y_log10() + theme_andy() + labs(color = "Significant\nChi Sq") + ylab("Log10 Deviance")
plot_st_pop_dev(base_dir, "GC_CG", "ALL", re=FALSE) + scale_y_log10() + theme_andy() + labs(color = "Significant\nChi Sq") + ylab("Log10 Deviance")
plot_st_pop_dev(base_dir, "cpg_GC_AT", "ALL", re=FALSE) + scale_y_log10() + theme_andy() + labs(color = "Significant\nChi Sq") + ylab("Log10 Deviance")
plot_st_pop_dev(base_dir, "cpg_GC_TA", "ALL", re=FALSE) + scale_y_log10() + theme_andy() + labs(color = "Significant\nChi Sq") + ylab("Log10 Deviance")
plot_st_pop_dev(base_dir, "cpg_GC_CG", "ALL", re=FALSE) + scale_y_log10() + theme_andy() + labs(color = "Significant\nChi Sq") + ylab("Log10 Deviance")

#' Load position-level results for subtype across all populations
#'
#' @param dir Location of the result files
#' @param st The mutation sub type
#' @returns A data frame with columns dev, singletons, controls, and offset
load_st_res <- function(dir, st, max_d = 20){
  df_ALL <- load_st_pop_res(dir, st, "ALL", max_d = max_d)
  df_AFR <- load_st_pop_res(dir, st, "AFR", max_d = max_d)
  df_AMR <- load_st_pop_res(dir, st, "AMR", max_d = max_d)
  df_EAS <- load_st_pop_res(dir, st, "EAS", max_d = max_d)
  df_EUR <- load_st_pop_res(dir, st, "EUR", max_d = max_d)
  df_SAS <- load_st_pop_res(dir, st, "SAS", max_d = max_d)
  df_ALL$pop <- "ALL"
  df_AFR$pop <- "AFR"
  df_AMR$pop <- "AMR"
  df_EAS$pop <- "EAS"
  df_EUR$pop <- "EUR"
  df_SAS$pop <- "SAS"
  df <- bind_rows(df_ALL, df_AFR) %>%
    bind_rows(df_AMR) %>%
    bind_rows(df_EAS) %>%
    bind_rows(df_EUR) %>%
    bind_rows(df_SAS) %>%
    arrange(offset, pop)
  return(df)
}

#' Plot the position-level results for a subtype across all populations
#'
#' @param dir Location of the result files
#' @param st The mutation sub type
#' @param re Whether to convert the deviance to relative entropy value
#' @returns ggplot object with the figure
plot_st_res <- function(dir, st, re = TRUE, max_d = 20, exclude_all = FALSE){
  df <- load_st_res(dir, st, max_d)
  st_text <- subtype_text(st)
  if(re){
    df <- df %>%
      rowwise() %>%
      mutate(stat_val = dev / (2*(singletons + controls)))
    stat_type <- "Relative Entropy"
  } else {
    df <- df %>% rename(stat_val = dev)
    stat_type <- "Deviance"
  }
  if(exclude_all) {
    df <- df %>% filter(pop != "ALL")
  }
  g_obj <- df %>%
    ggplot(aes(x = offset, y = stat_val, colour = pop)) +
    geom_point() +
    geom_line() +
    labs(x = "Relative Position",
         y = stat_type,
         color = "Population") +
    ggtitle("Influence by Position", paste0("Sub type: ", st_text))
  return(g_obj)
}

#' Load position-level results for all subtypes for a population
#'
#' @param dir Location of the result files
#' @param pop Which super-population
#' @returns A data frame with columns dev, singletons, controls, and offset
load_all_st_pop <- function(dir, pop, max_d = 1000){
  df_A_C <- load_st_pop_res(dir, "AT_CG", pop, max_d = max_d)
  df_A_G <- load_st_pop_res(dir, "AT_GC", pop, max_d = max_d)
  df_A_T <- load_st_pop_res(dir, "AT_TA", pop, max_d = max_d)
  df_C_A <- load_st_pop_res(dir, "GC_TA", pop, max_d = max_d)
  df_C_T <- load_st_pop_res(dir, "GC_AT", pop, max_d = max_d)
  df_C_G <- load_st_pop_res(dir, "GC_CG", pop, max_d = max_d)
  df_CpG_A <- load_st_pop_res(dir, "cpg_GC_TA", pop, max_d = max_d)
  df_CpG_T <- load_st_pop_res(dir, "cpg_GC_AT", pop, max_d = max_d)
  df_CpG_G <- load_st_pop_res(dir, "cpg_GC_CG", pop, max_d = max_d)

  df_A_C$type <- "AT_CG"
  df_A_G$type <- "AT_GC"
  df_A_T$type <- "AT_TA"
  df_C_A$type <- "GC_TA"
  df_C_T$type <- "GC_AT"
  df_C_G$type <- "GC_CG"
  df_CpG_A$type <- "cpg_GC_TA"
  df_CpG_T$type <- "cpg_GC_AT"
  df_CpG_G$type <- "cpg_GC_CG"

  df <- bind_rows(df_A_C, df_A_G) %>%
    bind_rows(df_A_T) %>%
    bind_rows(df_C_A) %>%
    bind_rows(df_C_T) %>%
    bind_rows(df_C_G) %>%
    bind_rows(df_CpG_A) %>%
    bind_rows(df_CpG_T) %>%
    bind_rows(df_CpG_G)
  return(df)
}

#' Plot the position-level results for all subtypes for a population
#'
#' @param dir Location of the result files
#' @param st The mutation sub type
#' @param re Whether to convert the deviance to relative entropy value
#' @returns ggplot object with the figure
plot_all_st_res <- function(dir, pop, re = TRUE, max_d = 20, exclude_all = FALSE){
  df <- load_all_st_pop(dir, pop, max_d = max_d)
  if(re){
    df <- df %>%
      rowwise() %>%
      mutate(stat_val = dev / (2*(singletons + controls)))
    stat_type <- "Relative Entropy"
  } else {
    df <- df %>% rename(stat_val = dev)
    stat_type <- "Deviance"
  }
  if(exclude_all) {
    df <- df %>% filter(pop != "ALL")
  }
  g_obj <- df %>%
    ggplot(aes(x = offset, y = stat_val, colour = type)) +
    geom_point() +
    geom_line() +
    labs(x = "Relative Position",
         y = stat_type,
         color = "Sub Type") +
    ggtitle("Influence by Position", paste0("Population: ", pop))
  return(g_obj)
}

# Nucleotide level results
#' Load nucleotide-level results for a population-subtype pair at a rp
#'
#' @param dir Location of the result files
#' @param st The mutation sub type
#' @param pop Which super-population
#' @param rp Relative position
#' @returns A data frame with columns dev, singletons, controls, and offset
load_st_pop_resid <- function(dir, st, pop, rp){
  file_loc <- paste0(dir, "resid/" , pop, "/", st, "_rp_", rp, ".csv")
  df <- read_csv(file_loc, show_col_types = FALSE)
  return(df)
}

plot_st_pop_resid <- function(dir, st, pop, rp){
  df <- load_st_pop_resid(dir, st, pop, rp)
  p <- df %>%
    ggplot(aes(x = status, y = nuc, fill = res)) +
    geom_tile() +
    scale_fill_distiller(palette = "PiYG") +
    ggtitle("Position Level Residuals",
            paste0("Subtype: ", subtype_text(st), "; Position: ", rp)) +
    xlab("") +
    ylab("Nucleotide")
  return(p)
}

#### Paper plots

#### Single subtype plots with significance colored in
plot_st_pop_dev(base_dir, "AT_GC", "ALL", re=FALSE) + scale_y_log10() + theme_andy() + labs(color = "Significance\nChi Sq")

#### All sub-types within a given population
df <- load_all_st_pop(base_dir, "ALL", max_d = 10)
df <- df %>%
  rowwise() %>%
  mutate(from = ifelse(str_starts(type, "AT"), "A", ifelse(str_starts(type, "cpg"), "CpG", "C")),
         to = ifelse(str_starts(type, "AT"), str_sub(type, 4, 4), ifelse(str_starts(type, "cpg"), str_sub(type, 9, 9), str_sub(type, 5, 5))),
         re = dev / (2 * (singletons + controls))) %>%
  mutate(st2 = subtype_text(type))

# Original spaghetti plot
df %>%
  filter(from == "A") %>%
  ggplot(aes(x = offset, y = dev, color = st2)) +
  geom_point() +
  geom_line() +
  ggtitle("Marginal Influence of Individual Positions") +
  xlab("Relative Position") +
  ylab("Deviance") +
  scale_colour_manual(values=cbPalette) +
  labs(color = "Subtype") +
  scale_y_log10() +
  theme_andy()

df %>%
  filter(from != "A") %>%
  ggplot(aes(x = offset, y = dev, color = st2)) +
  geom_point() +
  geom_line() +
  ggtitle("Marginal Influence of Individual Positions") +
  xlab("Relative Position") +
  ylab("Deviance") +
  scale_colour_manual(values=cbPalette) +
  labs(color = "Subtype") +
  scale_y_log10() +
  theme_andy()

# Grid arrangement
df %>%
  ggplot(aes(x = offset, y = re, group = type)) + geom_point() + geom_line() +
  facet_grid(rows = vars(to), cols = vars(from)) +
  ggtitle("Marginal Influence of Individual Positions") +
  xlab("Relative Position") +
  ylab("Relative Entropy") +
  theme_bw()

df %>%
  filter(from == "A") %>%
  ggplot(aes(x = offset, y = re, group = type)) + geom_point() + geom_line() +
  facet_grid(rows = vars(to), cols = vars(from)) +
  ggtitle("Marginal Influence of Individual Positions") +
  xlab("Relative Position") +
  ylab("Relative Entropy") +
  theme_andy()

df %>%
  filter(from != "A") %>%
  ggplot(aes(x = offset, y = re, group = type)) + geom_point() + geom_line() +
  facet_grid(rows = vars(to), cols = vars(from)) +
  ggtitle("Marginal Influence of Individual Positions") +
  xlab("Relative Position") +
  ylab("Relative Entropy") +
  theme_andy()
