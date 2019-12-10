#
# Plot QTL scans for main paper figure
#

library(tidyverse)
library(patchwork)

# change ggplot2 theme
theme_set(theme_bw() + theme(text = element_text(size = 12), panel.grid = element_blank()))


#
# read data ---------------------------------------------------------------
#
# read QTL scan results
qtl_scan <- read_csv("./data/processed/ril_qtl/qtl_scans.csv") %>% 
  filter(pheno %in% c("z68_shoots_mean_hn", "flowering_time_mean_hn", "flowering_time_mean_ln"))

# read permutation thresholds
qtl_thr <- read_csv("./data/processed/ril_qtl/qtl_permutation_thresholds.csv") %>% 
  filter(alpha == 0.05) %>% 
  filter(pheno %in% c("z68_shoots_mean_hn", "flowering_time_mean_hn", "flowering_time_mean_ln"))
  

# read peak positions
qtl_peaks <- read_csv("./data/processed/ril_qtl/qtl_peaks.csv")


#
# make plot ---------------------------------------------------------------
#
qtl_peaks %>% 
  ggplot() +
  geom_blank(data = qtl_scan, aes(x = pos)) +
  geom_segment(aes(x = ci_lo, xend = ci_hi, y = pheno, yend = pheno),
               size = 3) +
  facet_grid(. ~ paste0("Chr", chrom), scales = "free", space = "free") +
  scale_x_continuous(breaks = seq(0, 400, 100)) +
  theme(panel.grid.major = element_line(colour = "lightgrey", size = 0.2))




