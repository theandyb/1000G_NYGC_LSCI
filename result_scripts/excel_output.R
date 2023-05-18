## Generate Single Position Excel files
## One document per population
## One sheet per subtype x population pair

library(tidyverse)
library(writexl)

subtypes <- c("AT_CG", "AT_GC", "AT_TA",
              "GC_AT", "GC_TA", "GC_CG",
              "cpg_GC_AT", "cpg_GC_TA", "cpg_GC_CG")

result_dir <- "/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/single_pos/"

# example
df <- read_csv("output/single_pos/ALL/AT_CG.csv")
df %>% arrange(offset)

# Get all subtypes for ALL
res_list <- vector("list", length(subtypes))
names(res_list) <- subtypes

for(st in subtypes){
  res_list[[st]] <- read_csv(paste0(result_dir, "ALL", "/", st, ".csv"), show_col_types = FALSE) %>%
    arrange(offset)
}

write_xlsx(res_list, paste0(result_dir, "ALL", ".xlsx"))

for(pop in c("AFR", "AMR", "EAS", "EUR", "SAS")){
  res_list <- vector("list", length(subtypes))
  names(res_list) <- subtypes

  for(st in subtypes){
    res_list[[st]] <- read_csv(paste0(result_dir, pop, "/", st, ".csv"), show_col_types = FALSE) %>%
      arrange(offset)
  }
  write_xlsx(res_list, paste0(result_dir, pop, ".xlsx"))
}



### Generate the same files for two position residuals

result_dir <- "/net/snowwhite/home/beckandy/research/1000G_NYGC_LSCI/output/two_pos/"

res_list <- vector("list", length(subtypes))
names(res_list) <- subtypes

for(st in subtypes){
  res_list[[st]] <- read_csv(paste0(result_dir, "ALL", "/", st, ".csv"), show_col_types = FALSE) %>%
    arrange(rp1, rp2)
}

write_xlsx(res_list, paste0(result_dir, "ALL", ".xlsx"))

for(pop in c("AFR", "AMR", "EAS", "EUR", "SAS")){
  res_list <- vector("list", length(subtypes))
  names(res_list) <- subtypes

  for(st in subtypes){
    res_list[[st]] <- read_csv(paste0(result_dir, pop, "/", st, ".csv"), show_col_types = FALSE) %>%
      arrange(rp1, rp2)
  }
  write_xlsx(res_list, paste0(result_dir, pop, ".xlsx"))
}
