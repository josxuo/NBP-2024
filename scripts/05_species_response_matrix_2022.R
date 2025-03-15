########################################
### 05. 2022 species response matrix ###
########################################

#### The purpose of this script is to create species response matricies for
#### each NBP count circle in 2022 so they may be ordinated against traffic flow
#### data

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(readxl)
library(openxlsx)
library(data.table)

# load tidied nbp data
dat <- read_xlsx("data/b_intermediate_data/nbp_tidy_jan_24.xlsx")
traits <- read_xlsx("data/zz_miscellaneous_data/NBP_species_list+functional_traits.xlsx", sheet = "FunctionalTraits")
circs <- read.csv("data/c_analysis_data/nbp_circ_codes.csv") %>% mutate(sort.col = row_number())

# Filter and structure data to create 2022 species response matrices

## only want count circles with complete year of surveys for 2022
n12in22 <- dat %>%
  filter(year == 2022, exclusions == "none") %>% 
  group_by(station.code) %>%
  summarise(n = n_distinct(survey_id)) %>%
  ungroup() %>%
  filter(n == 12) %>%
  pull(station.code)

colnames(dat)

# generate max species response matrix for 2022
dsrm <- dat %>%
  ## filter out hybirds and sp. reports and filter for circles with 12 surveys in 2022
  filter(year == 2022, !str_detect(species, pattern = " sp\\.| x "),
         station.code %in% n12in22) %>%
  
  ## new colum summing observation by all modalities
  mutate(obs = seen + heard + fly) %>%
  
  ## group by count circle and species
  group_by(station.code, bird.code) %>%
  
  ## keep just the maximum abundance for each site doesn't make sense to add
  reframe(max.obs = max(obs)) %>%

  ## pivot wider with each species as a column and observed numbers as the value
  pivot_wider(names_from = bird.code, values_from = max.obs) %>%
  
  ## replace nas with zero (keeping in mind that zero doesn't necessarily mean zero...)
  replace(is.na(.), 0) %>%
  
  left_join(., select(circs, code, sort.col), join_by("station.code" == "code")) %>%
  arrange(sort.col) %>%
  select(order(colnames(.))) %>%
  select(station.code, everything(), -sort.col)

## inspect data frame
view(dsrm)  ## looks good



## ensure only exporting and analyzing count circles with covariates
covs <- read_csv("data/c_analysis_data/circ_no_overlap_covariates.csv") %>%
  na.exclude()

view(covs)

## add traffic flow data
traf <- read_xlsx("data/c_analysis_data/circ_traffic_flow.xlsx")

t <- select(traf, Station, m_arterial, ampk_arterial, traff_ind)

covars <- inner_join(covs, t)


srm.p <- dsrm %>%
  filter(station.code %in% covars$Station)

colind <- unname(apply(srm.p[, -1], 2, sum))

colind <- c(1, colind)
zs <- which(colind == 0)

srmf <- srm.p %>% select(-all_of(zs))

view(srmf)


covariates <- covars %>% filter(Station %in% srmf$station.code) %>%
  select(Station, everything(), -UniqueStationID, -Canopy.m2.2021, -Canopy.m2.2016, -CanopyProp2016) %>%
  left_join(., select(circs, code, sort.col), join_by("Station" == "code")) %>%
  arrange(sort.col) %>%
  select(-sort.col)

## functional trait data
func <- traits %>%
  na.exclude() %>% 
  filter(alpha.code %in% colnames(srmf)) %>%
  select(-com.name, -sci.name) %>%
  arrange(alpha.code)

func
colnames(srmf)
cbind(c(func$alpha.code, "hold"), colnames(srmf[, -1]))

## ready for export
write.xlsx(func, "data/c_analysis_data/npb_2022_traits.xlsx")
write.xlsx(srmf, "data/c_analysis_data/nbp_srm_2022.xlsx")
write.xlsx(covariates, "data/c_analysis_data/circle_covs_truncated.xlsx")


## What if I average across all years and take max?
# generate max species response matrix for 2022
dsrm <- dat %>%
  ## filter out hybirds and exclusions
  filter(year > 2004, !str_detect(species, pattern = " sp\\.| x "), exclusions == "none") %>%
  
  ## new colum summing observation by all modalities
  mutate(obs = seen + heard + fly) %>%
  
  ## group by count circle and species
  group_by(station.code, bird.code, year) %>%
  
  ## keep just the maximum abundance for each site doesn't make sense to add
  reframe(max.obs = max(obs)) %>%
  
  group_by(station.code, bird.code) %>%
  
  reframe(med.max.obs = median(max.obs)) %>%
  
  ## pivot wider with each species as a column and observed numbers as the value
  pivot_wider(names_from = bird.code, values_from = med.max.obs) %>%
  
  ## replace nas with zero (keeping in mind that zero doesn't necessarily mean zero...)
  replace(is.na(.), 0) %>%
  
  left_join(., select(circs, code, sort.col), join_by("station.code" == "code")) %>%
  arrange(sort.col) %>%
  select(order(colnames(.))) %>%
  select(station.code, everything(), -sort.col)

## inspect data frame
view(dsrm)  ## looks good



## ensure only exporting and analyzing count circles with covariates
covs <- read_csv("data/c_analysis_data/circ_no_overlap_covariates.csv") %>%
  na.exclude()

view(covs)

## add traffic flow data
traf <- read_xlsx("data/c_analysis_data/circ_traffic_flow.xlsx")

t <- select(traf, Station, m_arterial, ampk_arterial, traff_ind)

covars <- inner_join(covs, t)


srm.p <- dsrm %>%
  filter(station.code %in% covars$Station)

colind <- unname(apply(srm.p[, -1], 2, sum))

colind <- c(1, colind)
zs <- which(colind == 0)

srmf <- srm.p %>% select(-all_of(zs))

view(srmf)


covariates <- covars %>% filter(Station %in% srmf$station.code) %>%
  select(Station, everything(), -UniqueStationID, -Canopy.m2.2021, -Canopy.m2.2016, -CanopyProp2016) %>%
  left_join(., select(circs, code, sort.col), join_by("Station" == "code")) %>%
  arrange(sort.col) %>%
  select(-sort.col)

## functional trait data
func <- traits %>%
  na.exclude() %>% 
  filter(alpha.code %in% colnames(srmf)) %>%
  select(-com.name, -sci.name) %>%
  arrange(alpha.code)

func
colnames(srmf)
cbind(c(func$alpha.code, "hold"), colnames(srmf[, -1]))

## ready for export
write.xlsx(func, "data/c_analysis_data/npb_traits_median_max_obs.xlsx")
write.xlsx(srmf, "data/c_analysis_data/nbp_srm_median_max_obs.xlsx")
write.xlsx(covariates, "data/c_analysis_data/circle_covs_truncated_median_max_obs.xlsx")
