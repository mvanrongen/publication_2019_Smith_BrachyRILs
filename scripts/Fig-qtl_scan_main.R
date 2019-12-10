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
qtl_peaks <- read_csv("./data/processed/ril_qtl/qtl_peaks.csv") %>% 
  filter(pheno %in% c("z68_shoots_mean_hn", "flowering_time_mean_hn", "flowering_time_mean_ln"))
  


#
# make plot ---------------------------------------------------------------
#
# panel A
p1 <- qtl_scan %>% 
  filter(str_detect(pheno, "shoots")) %>% 
  ggplot(aes(pos, LOD)) +
  geom_rect(data = qtl_peaks %>% filter(str_detect(pheno, "shoots")),
            aes(xmin = ci_lo, xmax = ci_hi, ymin = -Inf, ymax = Inf), 
            fill = "brown", alpha = 0.3, inherit.aes = FALSE) +
  geom_line(aes(linetype = covariate), show.legend = FALSE) +
  geom_hline(data = qtl_thr %>% filter(str_detect(pheno, "shoots")), 
             aes(yintercept = threshold), linetype = 2) +
  facet_grid(nitrate ~ chrom, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0, 500, 100)) +
  scale_y_continuous(limits = c(0, 6.5), breaks = seq(0, 6, 2)) +
  labs(x = "Position (cM)", y = "LOD", tag = "A")

# panel B
p2 <- qtl_scan %>% 
  filter(str_detect(pheno, "flowering")) %>% 
  ggplot(aes(pos, LOD)) +
  geom_rect(data = qtl_peaks %>% filter(str_detect(pheno, "flowering")),
            aes(xmin = ci_lo, xmax = ci_hi, ymin = -Inf, ymax = Inf), 
            fill = "brown", alpha = 0.3, inherit.aes = FALSE) +
  geom_line(aes(linetype = covariate), show.legend = FALSE) +
  geom_hline(data = qtl_thr %>% filter(str_detect(pheno, "flowering")), 
             aes(yintercept = threshold), linetype = 2) +
  facet_grid(nitrate ~ chrom, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0, 500, 100)) +
  scale_y_continuous(limits = c(0, 6.5), breaks = seq(0, 6, 2)) +
  labs(x = "Position (cM)", y = "LOD", tag = "B")


pdf("./figures/Fig-qtl_scans_main.pdf", width = 7, height = 4)
p1 + p2 + plot_layout(ncol = 1, heights = c(1, 2))
dev.off()


