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
  filter(str_detect(pheno, "PC1"))

# read permutation thresholds
qtl_thr <- read_csv("./data/processed/ril_qtl/qtl_permutation_thresholds.csv") %>% 
  filter(alpha == 0.05) %>% 
  filter(str_detect(pheno, "PC1"))
  
# read peak positions
qtl_peaks <- read_csv("./data/processed/ril_qtl/qtl_peaks.csv")


#
# make plot ---------------------------------------------------------------
#
p1 <- qtl_scan %>% 
  ggplot(aes(pos, LOD)) +
  geom_rect(data = qtl_peaks %>% filter(str_detect(pheno, "PC1")),
            aes(xmin = ci_lo, xmax = ci_hi, ymin = -Inf, ymax = Inf), 
            fill = "brown", alpha = 0.3, inherit.aes = FALSE) +
  geom_line(aes(linetype = covariate), show.legend = FALSE) +
  geom_hline(data = qtl_thr, 
             aes(yintercept = threshold), linetype = 2) +
  facet_grid(nitrate ~ chrom, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0, 500, 100)) +
  scale_y_continuous(limits = c(0, 6.5), breaks = seq(0, 6, 2)) +
  labs(x = "Position (cM)", y = "LOD", tag = "A")

p2 <- qtl_peaks %>% 
  # mutate(pheno = str_remove(pheno, "_mean_hn|_mean_ln|_mean_plas")) %>% 
  mutate(pheno = str_remove(pheno, "_mean")) %>% 
  ggplot() +
  geom_blank(data = qtl_scan, aes(x = pos)) +
  geom_segment(aes(x = ci_lo, xend = ci_hi, y = pheno, yend = pheno, colour = nitrate),
               size = 3) +
  facet_grid(. ~ paste0("Chr", chrom), scales = "free", space = "free") +
  scale_x_continuous(breaks = seq(0, 400, 150)) +
  theme(panel.grid.major = element_line(colour = "lightgrey", size = 0.2)) +
  scale_colour_manual(values = c("HN" = "brown", "LN" = "steelblue", "PLAS" = "black")) +
  labs(x = "Position", y = "", tag = "B")

gridExtra::grid.arrange(p1, p2, ncol = 1)



