library(tidyverse)
library(broom)
library(patchwork)
theme_set(theme_classic() + theme(text = element_text(size = 16)))


# Read data ---------------------------------------------------------------

# read individual data
ril_phenotypes <- read_csv("data/processed/ril_phenotypes/ril_phenotypic_data_individual.csv",
                           col_types = cols(
                             .default = col_double(),
                             sheet = col_character(),
                             identifier = col_character(),
                             block = col_character(),
                             shelf = col_character(),
                             nitrate_level = col_factor(levels = c("LN", "HN")),
                             id = col_character(),
                             batch = col_factor(),
                             position_grid_location = col_character(),
                             ril_id = col_character(),
                             sow_date = col_date(format = "%Y-%m-%d"),
                             z68_date = col_date(format = "%Y-%m-%d"),
                             senescence_date = col_date(format = "%Y-%m-%d")
                           ))

# read QTL scan results
qtl_scan <- read_csv("./data/processed/ril_qtl/qtl_scans.csv") %>% 
  filter(str_detect(pheno, "PC1"))

# read permutation thresholds
qtl_thr <- read_csv("./data/processed/ril_qtl/qtl_permutation_thresholds.csv") %>% 
  filter(alpha == 0.05) %>% 
  filter(str_detect(pheno, "PC1"))

# read peak positions
qtl_peaks <- read_csv("./data/processed/ril_qtl/qtl_peaks.csv")


# Run PCA -----------------------------------------------------------------

pca_data <- ril_phenotypes %>% 
  select(batch, nitrate_level, ril_id, 
         z68_shoots, flowering_time, z68_greenness, senescence_height, 
         seed_weight_g, lifespan, flowering_senescence_interval) %>% 
  drop_na()

pca_result <- pca_data %>% 
  select(-batch, -nitrate_level, -ril_id) %>% 
  prcomp(scale = TRUE)

library(ggfortify)
autoplot(pca_result, data = pca_data,
         colour = "nitrate_level", loadings = TRUE, loadings.label = TRUE) +
  scale_colour_manual(values = c("steelblue", "firebrick")) +
  theme_classic() +
  coord_fixed(ratio = 0.3)

# eigenvalues
eigenvalues <- pca_result %>% 
  tidy(matrix = "eigenvalues") %>% 
  mutate(percent = percent*100, cumulative = cumulative*100)

# loadings
loadings <- pca_result %>% 
  tidy(matrix = "loadings")

# screeplot
eigenvalues %>% 
  ggplot(aes(factor(PC), percent)) +
  geom_col(fill = "grey") +
  geom_text(aes(label = round(percent, 1)), vjust = -0.2) +
  theme_classic()


# Figures -----------------------------------------------------------------

# PC plot
p1 <- pca_result %>% 
  augment(data = pca_data) %>% 
  ggplot(aes(.fittedPC1, .fittedPC2)) +
  geom_point(aes(colour = nitrate_level)) +
  scale_colour_manual(values = c("HN" = "brown", 
                                 "LN" = "steelblue", 
                                 "PLAS" = "black")) +
  labs(x = paste0("PC1 (", round(eigenvalues$percent[1], 1), "%)"), 
       y = paste0("PC2 (", round(eigenvalues$percent[2], 1), "%)"), 
       colour = "Nitrate") +
  theme_classic() +
  coord_fixed(ratio = 0.3)

# PC1 loadings
p2 <- loadings %>% 
  filter(PC %in% 1) %>% 
  mutate(column = fct_reorder(column, value)) %>% 
  ggplot(aes(value, column)) +
  geom_point() +
  geom_vline(xintercept = 0) +
  theme_classic() +
  theme(panel.grid = element_blank()) +
  labs(x = "PC1 Loading", y = "")

# QTL scan for PC1
p3 <- qtl_scan %>% 
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
  labs(x = "Position (cM)", y = "LOD")

p4 <- qtl_peaks %>% 
  # mutate(pheno = str_remove(pheno, "_mean_hn|_mean_ln|_mean_plas")) %>% 
  mutate(pheno = str_remove(pheno, "_mean")) %>% 
  ggplot() +
  geom_blank(data = qtl_scan, aes(x = pos)) +
  geom_segment(aes(x = ci_lo, xend = ci_hi, y = pheno, yend = pheno, 
                   colour = nitrate),
               size = 3) +
  facet_grid(. ~ paste0("Chr", chrom), scales = "free", space = "free") +
  scale_x_continuous(breaks = seq(0, 400, 150)) +
  theme(panel.grid.major = element_line(colour = "lightgrey", size = 0.2)) +
  scale_colour_manual(values = c("HN" = "brown", 
                                 "LN" = "steelblue", 
                                 "PLAS" = "black")) +
  labs(x = "Position", y = "")


pdf("figures/FigQTL4-qtl_scans_pca.pdf", width = 14, height = 10)
(
  ((p1 | p2) + plot_layout(widths = c(2, 1))) / 
    (p3 / p4)
) +
  plot_layout(heights = c(1, 3)) +
  plot_annotation(tag_levels = "A")
dev.off()
