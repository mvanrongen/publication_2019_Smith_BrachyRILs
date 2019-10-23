#
# Find genes overlapping QTL intervals
#

library(biomaRt)
library(readr)
library(dplyr)
library(valr)



#
# Fetch genes from BioMart ------------------------------------------------
#
# List of BioMart plant databases
listMarts(host="plants.ensembl.org")

# Specify the plants_mart
m <- useMart("plants_mart", host="plants.ensembl.org")

# What datasets exist (we want brachypodium distachyon)
listDatasets(m) 

# Define dataset
m <- useMart("plants_mart", host="plants.ensembl.org", dataset="bdistachyon_eg_gene")

# # define GO terms of interest
# nitrate_go_terms <- c("GO:0031667", "GO:0043602", "GO:0090352", "GO:0090408", "GO:0090548", "GO:0042126", "GO:0042128", "GO:0042594", "GO:1902025", "GO:0015706", "GO:0071249", "GO:0010167")
# shoot_go_terms <- c("GO:0090506", "GO:2000032", "GO:0010223", "GO:0048831", "GO:0010016", "GO:0048367", "GO:1900618")
# 
# # download from BioMart
# go <- getBM(attributes=c("ensembl_gene_id", "tair_symbol", "external_gene_name",
#                          "chromosome_name", "start_position", "end_position",
#                          "go_id", "go_linkage_type", "name_1006"), 
#             mart=m,
#             filters = "go",
#             values = c(nitrate_go_terms, shoot_go_terms)) 

go <- getBM(attributes=c("ensembl_gene_id", "external_gene_name",
                         "chromosome_name", "start_position", "end_position",
                         "go_id", "go_linkage_type", "name_1006",
                         "description"), 
            mart=m)

# collapse GO term information so that we have one gene per row
go <- go %>% 
  group_by_at(vars(ensembl_gene_id:end_position, description)) %>% 
  summarise(go_id = paste(unique(go_id), collapse = "; "),
            go_description = paste(unique(name_1006), collapse = "; ")) %>% 
  ungroup() %>% 
  mutate_all(~ na_if(., "")) %>% 
  rename(chrom = chromosome_name, start = start_position, end = end_position)


#
# intersect with QTL intervals --------------------------------------------
#
# read qtl intervals
peaks <- read_csv("./data/processed/ril_qtl/qtl_peaks.csv") %>% 
  select(pheno, nitrate, covariate, chrom, pos_bp, start = ci_lo_bp, end = ci_hi_bp)

# intersecting
genes_overlap <- peaks %>% 
  mutate(chrom = as.character(chrom)) %>% 
  group_by(pheno, nitrate, covariate, chrom, pos_bp) %>% 
  bed_intersect(go, suffix = c("", "_gene")) %>% 
  rename(ensembl_gene_id = ensembl_gene_id_gene, 
         external_gene_name = external_gene_name_gene,
         go_id = go_id_gene,
         go_description = go_description_gene) %>% 
  mutate(dist_peak = pmin(start - (end_gene - start_gene)/2, end - (end_gene - start_gene)/2))

# save table  
genes_overlap %>% 
  write_csv("./data/processed/ril_qtl/qtl_peaks_gene_overlap.csv")
