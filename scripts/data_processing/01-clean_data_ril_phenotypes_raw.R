###
# Script to process and tidy raw Excel file
###

# load packages
library(tidyverse)
library(lubridate)
library(readxl)


# read data ---------------------------------------------------------------

# Path to original data file
path <- file.path("./data/raw/2019_08_18_all_RIL_data_reworked.xlsx")

# data are spread over multiple worksheets
ril_phenotypes <- path %>%
  excel_sheets() %>%
  set_names() %>% 
  map_df(~ read_excel(path = path, sheet = .x), .id = "sheet")

# remove some empty columns that were exported from excel
ril_phenotypes <- ril_phenotypes %>% 
  select(-("...25":"...33"))

# tidy up columns
ril_phenotypes <- ril_phenotypes %>% 
  # make consistent column naming
  rename(id = Position_ID,
         ril_id = Line,
         z68_greenness = Greeness_Z68,
         z68_senescence_interval = `Z68-Senescence_interval`) %>% 
  # make all column names lowercase
  rename_all(tolower) %>% 
  # rename nitrate levels
  mutate(nitrate_level = replace(nitrate_level, nitrate_level == "L", "LN")) %>% 
  mutate(nitrate_level = replace(nitrate_level, nitrate_level == "H", "HN"))



# fix data quality issues-----------------------------------------------------

# fix one of the plants that was missing identifier
ril_phenotypes <- ril_phenotypes %>%
  mutate(batch = ifelse(identifier == "B8-HN-049", 8, batch))

# remove columns that relate to flowering, as the original had errors
# we re-calculated below
ril_phenotypes <- ril_phenotypes %>% 
  select(-z68_das, -senescence_das, -z68_senescence_interval)

# make date columns and calculate related traits
ril_phenotypes <- ril_phenotypes %>% 
  # make date columns
  mutate(sow_date = make_date(sow_year, sow_month, sow_day),
         senescence_date = make_date(senescence_year, senescence_month, senescence_day),
         z68_date = make_date(z68_year, z68_month, z68_day)) %>% 
  # calculate date-related traits
  mutate(lifespan = as.numeric(senescence_date - sow_date),
         flowering_time = as.numeric(z68_date - sow_date),
         flowering_senescence_interval = as.numeric(senescence_date - z68_date))

# remove extreme outlier for height 
# (over 60cm, which is 10 SDs away from the population's mean!)
ril_phenotypes <- ril_phenotypes %>% 
  mutate(senescence_height = ifelse(scale(senescence_height) > 10, NA, senescence_height))

# remove one line that only had data for one nitrate treatment (and one replicate)
ril_phenotypes <- ril_phenotypes %>% 
  group_by(ril_id) %>% 
  filter(n_distinct(nitrate_level) == 2) %>% 
  ungroup()

# remove line that was used to fill up space in the growth chamber
# (not part of the RIL population)
ril_phenotypes <- ril_phenotypes %>% 
  filter(ril_id != "Bd21-3")

# remove lines that didn't flower nor have shoots measured
ril_phenotypes <- ril_phenotypes %>%
  drop_na(sow_date, z68_shoots)


# tests of data quality ----------------------------------------------------

# check for duplicates
{
  temp <- ril_phenotypes %>% 
    count(nitrate_level, position_grid_location, ril_id, batch) %>% 
    filter(n != 1)
  
  if(nrow(temp) != 0) stop("there are duplicate entries in the data")
  rm(temp)
}

# check identifier columns are correct
if(nrow(ril_phenotypes) != length(unique(ril_phenotypes$identifier))){
  stop("identifier column contains duplicates")
}



# summarise data ----------------------------------------------------------

# Simple summary of data
ril_summary <- ril_phenotypes %>% 
  group_by(nitrate_level, ril_id) %>% 
  summarise_at(vars(z68_shoots, z68_greenness, senescence_height, seed_weight_g, 
                    lifespan, flowering_senescence_interval, flowering_time),
               .funs = list(mean = ~ mean(., na.rm = TRUE),
                            median = ~ median(., na.rm = TRUE),
                            n = ~ sum(!is.na(.)),
                            sd = ~ sd(., na.rm = TRUE))) %>% 
  # remove grouping
  ungroup()


# export data -------------------------------------------------------------

# all individual data
write_csv(ril_phenotypes, 
          path = "./data/processed/ril_phenotypes/ril_phenotypic_data_individual.csv", 
          na = "")

# summarised data
write_csv(ril_summary,
          path = "./data/processed/ril_phenotypes/ril_phenotypic_data_averaged.csv")





