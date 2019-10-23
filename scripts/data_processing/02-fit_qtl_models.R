##################
# Fit QTL models #
##################
# This script can be run non-interactively from the working directory

# load libraries
library(qtl2)
library(qtl2helper) # https://github.com/tavareshugo/qtl2helper
library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(tibble)


#
# Read phenotype data -----------------------------------------------------
#

# read summarised phenotype data
pheno <- read_csv("./data/processed/ril_phenotypes/ril_phenotypic_data_averaged.csv",
                  col_types = cols(
                    .default = col_double(),
                    nitrate_level = col_factor(levels = c("LN", "HN")),
                    ril_id = col_character()))

# Get HN and LN into separate columns as well as plasticity
pheno_summary <- pheno %>%
  filter(grepl("RIL", ril_id)) %>%  # retain only RILs (not the parents)
  select(ril_id, nitrate_level, matches("_mean")) %>%
  drop_na() %>%
  group_by(ril_id) %>%
  summarise_at(vars(matches("_mean")),
               .funs = list(hn = ~ .[nitrate_level == "HN"],
                            ln = ~ .[nitrate_level == "LN"],
                            plas = ~ .[nitrate_level == "HN"] - .[nitrate_level == "LN"])) %>%
  # scale the data (this is probably not really needed...)
  mutate_at(vars(matches("_mean")),
            list(~ (. - mean(.))/sd(.)))


#
# Read genotype data ------------------------------------------------------
#

# # prepare directory structure
# if(!dir.exists("./data/external/")) dir.create("./data/external", recursive = TRUE)
# 
# # Download genotype data
# if(!file.exists("./data/external/Cui2012_qtl2cross.zip")){
#   download.file("https://github.com/tavareshugo/qtl2data/raw/master/brachy_Cui2012/Cui2012_qtl2cross.zip",
#                 destfile = "./data/external/Cui2012_qtl2cross.zip",
#                 method = "wget")
# }

# read data as a cross2 object
ril_cross <- read_cross2("./data/external/Cui2012_qtl2cross.zip")

# Add phenotypes to cross object
ril_cross <- add_pheno(ril_cross, pheno_summary)

# insert pseudo-markers for evenly spaced scan
geno_map <- insert_pseudomarkers(ril_cross$gmap, step = 0.5)

# calculate genotype probabilities - assume 1% error rate
geno_prob <- calc_genoprob(ril_cross, geno_map, error_prob = 0.01, cores = 3)


#
# Simple QTL scan ------------------------------------------------------------
#

# Run scan
pheno_scan <- scan1(geno_prob, ril_cross$pheno, cores = 3)

# Run permutations (takes a while)
pheno_perm <- scan1perm(geno_prob, ril_cross$pheno, n_perm = 1000, cores = 3)

# get peaks using permutation thresholds
pheno_peaks <- find_peaks(pheno_scan, geno_map, 
                          threshold = summary(pheno_perm)[1, ], 
                          drop = 1, peakdrop = 1,
                          expand2markers = TRUE)


#
# QTL scan with covariate -------------------------------------------------
#

# Run scan
covar_scan <- scan1(geno_prob, 
                    pheno = ril_cross$pheno[, c("z68_shoots_mean_hn"), drop = FALSE], 
                    addcovar = ril_cross$pheno[, c("flowering_time_mean_hn"), drop = FALSE], 
                    cores = 3)

covar_perm <- scan1perm(geno_prob, 
                        pheno = ril_cross$pheno[, c("z68_shoots_mean_hn"), drop = FALSE], 
                        addcovar = ril_cross$pheno[, c("flowering_time_mean_hn"), drop = FALSE], 
                        n_perm = 1000,
                        cores = 3)

covar_peaks <- find_peaks(covar_scan, geno_map, 
                          threshold = summary(covar_perm)[1, ], 
                          drop = 1, peakdrop = 1,
                          expand2markers = TRUE)


#
# Tidy scans -------------------------------------------------------------
#
# coerce scan1 object to data.frame
pheno_scan_df <- pheno_scan %>% 
  # tidy (from qtl2helper package)
  tidy(map = geno_map) %>% 
  # add nitrate variable based on pheno column
  mutate(nitrate = toupper(str_remove(pheno, ".*_")),
         covariate = FALSE)

# add covariate scan
pheno_scan_df <- covar_scan %>% 
  tidy(geno_map) %>% 
  # add nitrate variable based on pheno column
  mutate(nitrate = toupper(str_remove(pheno, ".*_")),
         covariate = TRUE) %>% 
  bind_rows(pheno_scan_df)


# coerce scan1perm object to data.frame
pheno_perm_df <- pheno_perm %>% 
  # tidy (from qtl2helper package)
  tidy(alpha = c(0.1, 0.05, 0.01)) %>% 
  # create variable with nitrate treatment
  mutate(nitrate = toupper(str_remove(pheno, ".*_")),
         covariate = FALSE)

# add covariate permutation result
pheno_perm_df <- covar_perm %>% 
  # tidy (from qtl2helper package)
  tidy(alpha = c(0.1, 0.05, 0.01)) %>% 
  # create variable with nitrate treatment
  mutate(nitrate = toupper(str_remove(pheno, ".*_")),
         covariate = TRUE) %>% 
  bind_rows(pheno_perm_df)


#
# Tidy peaks --------------------------------------------------------------
#
# get markers into a named vector (useful for later use)
markers <- get_markers(ril_cross) %>% 
  select(marker, physical_pos) %>% 
  deframe()

# tidy peaks table
pheno_peaks_df <- pheno_peaks %>% 
  rename(pheno = lodcolumn, chrom = chr) %>% 
  mutate(nitrate = toupper(str_remove(pheno, ".*_")),
         covariate = FALSE)

# add peaks from covariate scan
pheno_peaks_df <- covar_peaks %>% 
  rename(pheno = lodcolumn, chrom = chr) %>% 
  mutate(nitrate = toupper(str_remove(pheno, ".*_")),
         covariate = TRUE) %>% 
  bind_rows(pheno_peaks_df)

# add marker IDs and bp positions
pheno_peaks_df <- pheno_peaks_df %>% 
  mutate(peak_marker = find_marker(ril_cross$gmap, chrom, pos),
         ci_lo_marker = find_marker(ril_cross$gmap, chrom, ci_lo),
         ci_hi_marker = find_marker(ril_cross$gmap, chrom, ci_hi)) %>% 
  mutate(pos_bp = markers[peak_marker],
         ci_lo_bp = markers[ci_lo_marker],
         ci_hi_bp = markers[ci_hi_marker]) %>% 
  select(-ci_lo_marker, -ci_hi_marker)


#
# Save output -------------------------------------------------------------
#
pheno_scan_df %>% 
  write_csv("./data/processed/ril_qtl/qtl_scans.csv")

pheno_perm_df %>% 
  write_csv("./data/processed/ril_qtl/qtl_permutation_thresholds.csv")

pheno_peaks_df %>% 
  write_csv("./data/processed/ril_qtl/qtl_peaks.csv")
