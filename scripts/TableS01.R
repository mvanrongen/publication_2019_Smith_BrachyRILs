library(tidyverse)

# read data
peaks <- read_csv("./data/processed/ril_qtl/qtl_peaks.csv")

# pretty output
peaks_pretty <- peaks %>% 
  mutate(pheno = str_remove(pheno, "_mean_hn|_mean_ln|_mean_plas")) %>% 
  mutate(pheno = case_when(pheno == "z68_shoots" ~ "Shoots",
                           pheno == "senescence_height" ~ "Height", 
                           pheno == "flowering_senescence_interval" ~ "Flowering to\nSenescence",
                           pheno == "flowering_time" ~ "Flowering", 
                           pheno == "lifespan" ~ "Lifespan")) %>% 
  mutate(pheno = ifelse(covariate, paste(pheno, "flowering", sep = "~"), pheno)) %>% 
  mutate(ci_lo_bp = round(ci_lo_bp/1e6, 2), 
         ci_hi_bp = round(ci_hi_bp/1e6, 2),
         pos_bp = round(pos_bp/1e6, 2)) %>% 
  select(Trait = pheno, 
         Nitrate = nitrate,
         Chromosome = chrom, 
         `Position (cM)` = pos, 
         `Lower CI (cM)` = ci_lo,
         `Upper CI (cM)` = ci_hi,
         `Position (Mbp)` = pos_bp,
         `Lower CI (Mbp)` = ci_lo_bp,
         `Upper CI (Mbp)` = ci_hi_bp) %>% 
  arrange(Trait, Nitrate, Chromosome, `Position (cM)`)

# print as a figure
peaks_pretty %>% 
  gridExtra::grid.table()

