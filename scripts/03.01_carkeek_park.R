##########################
### 03.01 Carkeek Park ###
##########################

# Park-level analysis of Neighborhood Bird Project data at Carkeek Park, Seattle

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(readxl)

# Custom bcs plotting theme 
source("functions/theme_bcs.R")

# load tidied NBP data (see script "01_tidy_raw_nbp.R" for details)
dat <- read_excel("data/processed/nbp_tidy_jan_24.xlsx")

# load non-overlapping
covs <- read_csv("data/processed/circ_no_overlap_covariates.csv")

# filter and prepare data for Carkeek Park analsyis
focal.park <- "Carkeek Park"

d.full <- dat %>% filter(park == focal.park)  # Observations from all stations, including spuhs
d.short <- dat %>% filter(park == focal.park, # Observations only from non-overlapping count circles and excluding spuhs
                          station.code %in% covs$station.code,
                          !str_detect(species, pattern = "sp\\.|Domestic|Spotted Owl")) %>%
  mutate(short = "y")

short.list <- d.short %>% dplyr::select(survey_id, short) %>% unique()

d.full <- left_join(d.full, short.list, join_by("survey_id" == "survey_id")) %>%
  mutate(short = ifelse(is.na(short), "n", short))

view(d.full)

# Survey stats

## How many surveys in the full dataset
nsurv.full <- 

