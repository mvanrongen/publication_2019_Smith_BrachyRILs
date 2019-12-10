#
# Plot QTL scans for main paper figure
#

library(qtl2)
library(tidyverse)
library(ggbeeswarm)
library(patchwork)
library(rsample)


# change ggplot2 theme
theme_set(theme_bw() + theme(text = element_text(size = 12), panel.grid = element_blank()))


#
# read data ---------------------------------------------------------------
#
# read phenotype summary
pheno <- read_csv("./data/processed/ril_phenotypes/ril_phenotypic_data_averaged.csv",
                  col_types = cols(
                    .default = col_double(),
                    nitrate_level = col_factor(levels = c("LN", "HN")),
                    ril_id = col_character()))

# read QTL peak positions
qtl_peaks <- read_csv("./data/processed/ril_qtl/qtl_peaks.csv") %>% 
  filter(pheno %in% c("flowering_time_mean_ln", "flowering_time_mean_hn"))
  
# read cross object
ril_cross <- read_cross2("./data/external/Cui2012_qtl2cross.zip")

# calculate genotype probabilities at each marker - assume 1% error rate
geno_prob <- calc_genoprob(ril_cross, ril_cross$gmap, error_prob = 0.01, cores = 3)


#
# Get marker genotypes --------------------------------------------------------
#
# Function to get marker genotypes in a tabular format
get_genoprobs <- function(target_marker, geno_prob){
  pull_genoprobpos(geno_prob, target_marker) %>% 
    as_tibble(rownames = "ril_id") %>% 
    mutate(genotype = case_when(aa <= 0.25 ~ "Bd3-1",
                                aa >= 0.75 ~ "Bd21", 
                                TRUE ~ as.character(NA))) %>% 
    select(ril_id, aa, genotype) %>% 
    mutate(marker = target_marker)
}

genos <- qtl_peaks %>% 
  pull(peak_marker) %>% 
  map_dfr(get_genoprobs, geno_prob = geno_prob)


#
# merge genotypes with phenotypes ----------------------------
#
parentals <- pheno %>% 
  filter(ril_id %in% c("Bd21", "Bd3-1")) %>% 
  select(ril_id, nitrate_level, z68_shoots_mean, flowering_time_mean) %>% 
  mutate(nitrate_level = tolower(nitrate_level)) %>% 
  pivot_wider(names_from = nitrate_level, 
              values_from = c(z68_shoots_mean, flowering_time_mean)) %>% 
  rename(genotype = ril_id)

rils <- pheno %>% 
  # retain only RILs with genotypes
  filter(ril_id %in% ind_ids(ril_cross)) %>% 
  select(ril_id, nitrate_level, z68_shoots_mean, flowering_time_mean) %>% 
  mutate(nitrate_level = tolower(nitrate_level)) %>% 
  pivot_wider(names_from = nitrate_level, 
              values_from = c(z68_shoots_mean, flowering_time_mean)) %>% 
  # add genotypes
  inner_join(genos, by = "ril_id")


#
# Make plots ------------------
#
p1 <- rils %>% 
  filter(marker == "BD3660_1") %>% 
  drop_na(genotype) %>% 
  ggplot(aes(genotype, flowering_time_mean_hn, colour = genotype)) +
  geom_quasirandom(alpha = 0.5) +
  geom_linerange(stat = "summary", fun.data = "mean_cl_boot", colour = "black", size = 1) +
  #geom_label(data = parentals, aes(label = genotype), nudge_x = c(0.5, -0.5)) +
  geom_point(data = parentals, shape = 8, position = position_nudge(x = c(0.5, -0.5)), size = 3) +
  scale_colour_manual(values = c("Bd21" = "#0072B2", "Bd3-1" = "#CC79A7"), 
                      na.value = "grey") +
  labs(x = "Genotype", y = "Flowering HN", tag = "A") +
  theme(legend.position = "none")

p1.1 <- rils %>% 
  filter(marker == "BD3660_1") %>% 
  drop_na(genotype) %>% 
  bootstraps(times = 1000, strata = genotype) %>% 
  mutate(dif = map_dbl(splits, function(x){
    x <- as.data.frame(x)
    bd21 <- x$flowering_time_mean_hn[x$genotype == "Bd21"]
    bd3.1 <- x$flowering_time_mean_hn[x$genotype == "Bd3-1"]
    mean(bd3.1) - mean(bd21)
  })) %>% 
  ggplot(aes(dif)) + 
  geom_density(fill = "grey", colour = "grey") +
  geom_vline(xintercept = 39 - 29.2, size = 1) +
  geom_vline(xintercept = 0, linetype = 2) +
  coord_flip() +
  labs(x = "Bd3-1 minus Bd21", y = "") +
  theme(axis.text.x = element_blank())

p2 <- rils %>% 
  filter(marker == "BD0836_1") %>% 
  drop_na(genotype) %>% 
  ggplot(aes(genotype, flowering_time_mean_hn, colour = genotype)) +
  geom_quasirandom(alpha = 0.5) +
  geom_linerange(stat = "summary", fun.data = "mean_cl_boot", colour = "black", size = 1) +
  #geom_label(data = parentals, aes(label = genotype), nudge_x = c(0.5, -0.5)) +
  geom_point(data = parentals, shape = 8, position = position_nudge(x = c(0.5, -0.5)), size = 3) +
  scale_colour_manual(values = c("Bd21" = "#0072B2", "Bd3-1" = "#CC79A7"), 
                      na.value = "grey") +
  labs(x = "Genotype", y = "Flowering HN", tag = "B") +
  theme(legend.position = "none")

p2.1 <- rils %>% 
  filter(marker == "BD0836_1") %>% 
  drop_na(genotype) %>% 
  bootstraps(times = 1000, strata = genotype) %>% 
  mutate(dif = map_dbl(splits, function(x){
    x <- as.data.frame(x)
    bd21 <- x$flowering_time_mean_hn[x$genotype == "Bd21"]
    bd3.1 <- x$flowering_time_mean_hn[x$genotype == "Bd3-1"]
    mean(bd3.1) - mean(bd21)
  })) %>% 
  ggplot(aes(dif)) + 
  geom_density(fill = "grey", colour = "grey") +
  geom_vline(xintercept = 39 - 29.2, size = 1) +
  geom_vline(xintercept = 0, linetype = 2) +
  coord_flip() +
  labs(x = "Bd3-1 minus Bd21", y = "") +
  theme(axis.text.x = element_blank())

p3 <- rils %>% 
  filter(marker == "BD0635_1") %>% 
  drop_na(genotype) %>% 
  ggplot(aes(genotype, flowering_time_mean_ln, colour = genotype)) +
  geom_quasirandom(alpha = 0.5) +
  geom_linerange(stat = "summary", fun.data = "mean_cl_boot", colour = "black", size = 1) +
  #geom_label(data = parentals, aes(label = genotype), nudge_x = c(0.5, -0.5)) +
  geom_point(data = parentals, shape = 8, position = position_nudge(x = c(0.5, -0.5)), size = 3) +
  scale_colour_manual(values = c("Bd21" = "#0072B2", "Bd3-1" = "#CC79A7"), 
                      na.value = "grey") +
  labs(x = "Genotype", y = "Flowering LN", tag = "C") +
  theme(legend.position = "none")


p3.1 <- rils %>% 
  filter(marker == "BD0635_1") %>% 
  drop_na(genotype) %>% 
  bootstraps(times = 1000, strata = genotype) %>% 
  mutate(dif = map_dbl(splits, function(x){
    x <- as.data.frame(x)
    bd21 <- x$flowering_time_mean_hn[x$genotype == "Bd21"]
    bd3.1 <- x$flowering_time_mean_hn[x$genotype == "Bd3-1"]
    mean(bd3.1) - mean(bd21)
  })) %>% 
  ggplot(aes(dif)) + 
  geom_density(fill = "grey", colour = "grey") +
  geom_vline(xintercept = 31.00 - 25.75, size = 1) +
  geom_vline(xintercept = 0, linetype = 2) +
  coord_flip() +
  labs(x = "Bd3-1 minus Bd21", y = "") +
  theme(axis.text.x = element_blank())

pdf("./figures/FigS-qtl_flowering_effects.pdf", width = 10, height = 2.5)
{p1 + p1.1 + plot_layout(ncol = 2, widths = c(2, 1))} | 
  {p2 + p2.1 + plot_layout(ncol = 2, widths = c(2, 1))} |
  {p3 + p3.1 + plot_layout(ncol = 2, widths = c(2, 1))}
dev.off()


# Get R2 for paper
rils %>% 
  filter(marker == "BD3660_1") %>% 
  lm(flowering_time_mean_hn ~ aa, data = .) %>% 
  summary()
