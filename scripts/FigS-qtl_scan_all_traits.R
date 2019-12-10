#
# Plot all QTL scans 
#

library(tidyverse)
theme_set(theme_bw() + theme(text = element_text(size = 12), panel.grid = element_blank()))

# read QTL scan results
qtl_scan <- read_csv("./data/processed/ril_qtl/qtl_scans.csv") %>% 
  # tidy nitrate variable
  mutate(nitrate = ifelse(nitrate == "PLAS", "Plasticity", nitrate)) %>% 
  mutate(nitrate = factor(nitrate, levels = c("LN", "HN", "Plasticity"))) %>% 
  # tidy pheno column
  mutate(pheno = str_remove(pheno, "_mean_hn|_mean_ln|_mean_plas")) %>% 
  mutate(pheno = case_when(pheno == "z68_shoots" ~ "Shoots",
                           pheno == "z68_greenness" ~ "Greenness", 
                           pheno == "senescence_height" ~ "Height",
                           pheno == "seed_weight_g" ~ "Seed weight", 
                           pheno == "flowering_senescence_interval" ~ "Flowering to Senescence",
                           pheno == "flowering_time" ~ "Flowering", 
                           pheno == "lifespan" ~ "Lifespan"))

# read permutation thresholds
qtl_thr <- read_csv("./data/processed/ril_qtl/qtl_permutation_thresholds.csv") %>% 
  # tidy nitrate variable
  mutate(nitrate = ifelse(nitrate == "PLAS", "Plasticity", nitrate)) %>% 
  mutate(nitrate = factor(nitrate, levels = c("LN", "HN", "Plasticity"))) %>% 
  # tidy pheno column
  mutate(pheno = str_remove(pheno, "_mean_hn|_mean_ln|_mean_plas")) %>% 
  mutate(pheno = case_when(pheno == "z68_shoots" ~ "Shoots",
                           pheno == "z68_greenness" ~ "Greenness", 
                           pheno == "senescence_height" ~ "Height",
                           pheno == "seed_weight_g" ~ "Seed weight", 
                           pheno == "flowering_senescence_interval" ~ "Flowering to Senescence",
                           pheno == "flowering_time" ~ "Flowering", 
                           pheno == "lifespan" ~ "Lifespan"))
  

# read peak positions
qtl_peaks <- read_csv("./data/processed/ril_qtl/qtl_peaks.csv") %>% 
  # tidy nitrate variable
  mutate(nitrate = ifelse(nitrate == "PLAS", "Plasticity", nitrate)) %>% 
  mutate(nitrate = factor(nitrate, levels = c("LN", "HN", "Plasticity"))) %>% 
  # tidy pheno column
  mutate(pheno = str_remove(pheno, "_mean_hn|_mean_ln|_mean_plas")) %>% 
  mutate(pheno = case_when(pheno == "z68_shoots" ~ "Shoots",
                           pheno == "z68_greenness" ~ "Greenness", 
                           pheno == "senescence_height" ~ "Height",
                           pheno == "seed_weight_g" ~ "Seed weight", 
                           pheno == "flowering_senescence_interval" ~ "Flowering to Senescence",
                           pheno == "flowering_time" ~ "Flowering", 
                           pheno == "lifespan" ~ "Lifespan"))
  


#
# make plot ---------------------------------------------------------------
#

plot_qtl <- function(trait){
  
  if(!(trait %in% qtl_scan$pheno)) stop("trait is not valid.")
  
  # get trait of interest
  scan <- qtl_scan[which(qtl_scan$pheno == trait), ]
  threshold <- qtl_thr[which(qtl_thr$pheno == trait & qtl_thr$alpha == 0.05), ]
  peaks <- qtl_peaks[which(qtl_peaks$pheno == trait), ]
  
  # make plot
  p <- ggplot(scan, aes(pos, LOD)) +
    geom_line(aes(linetype = covariate), show.legend = FALSE) +
    geom_hline(data = threshold, aes(yintercept = threshold), linetype = 2) +
    facet_grid(nitrate ~ chrom, scales = "free_x", space = "free_x") +
    scale_x_continuous(breaks = seq(0, 500, 100)) +
    scale_y_continuous(limits = c(0, 6.5), breaks = seq(0, 6, 2)) +
    labs(x = "Position (cM)", y = "LOD", title = trait)
  
  # Add peak intervals (if there's any)
  if(nrow(peaks) != 0){
    p <- p +
      geom_rect(data = peaks,
                   aes(xmin = ci_lo, xmax = ci_hi, ymin = -Inf, ymax = Inf), 
                   fill = "brown", alpha = 0.3, inherit.aes = FALSE)
  }
  
  print(p)
}

pdf("./figures/FigS-qtl_scans_all_traits.pdf", width = 7, height = 3.5)
qtl_scan %>% 
  pull(pheno) %>% 
  unique() %>% 
  sort() %>% 
  walk(plot_qtl)
dev.off()
