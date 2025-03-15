#############################
####### 00. occupancy #######
#############################

# clear environment

rm(list = ls())

# libraries
library(tidyverse)
library(unmarked)

# Use a consistent, branded theme for plots----
bcs_colors <- c(
  "dg" = "#0A3C23",
  "c" = "#FAF5F0",
  "yg" = "#E6FF55",
  "p" = "#FFB98C",
  "bg" = "#36BA3A"
)

## Create a custom theme function
theme_bcs <- function() {
  theme(
    # Backgrounds
    panel.background = element_rect(fill = bcs_colors["c"], color = NA),
    plot.background = element_rect(fill = bcs_colors["c"], color = NA),
    panel.grid.major = element_line(color = bcs_colors["dg"], linewidth = 0.5, linetype = "dotted"),
    panel.grid.minor = ggplot2::element_blank(),
    
    # Text
    text = element_text(color = bcs_colors["dg"]),
    axis.text = element_text(color = bcs_colors["dg"]),
    axis.title = element_text(color = bcs_colors["dg"]),
    plot.title = element_text(color = bcs_colors["dg"], face = "bold", size = 28),
    plot.subtitle = element_text(color = bcs_colors["dg"], face = "italic", size = 22),
    plot.caption = element_text(color = bcs_colors["dg"], size = 8),
    
    # Lines and borders
    axis.line = element_line(color = bcs_colors["dg"]),
    axis.ticks = element_line(color = bcs_colors["dg"]),
    panel.border = element_rect(color = bcs_colors["dg"], fill = NA),
    
    # Legends
    legend.background = element_rect(fill = bcs_colors["c"]),
    legend.key = element_rect(fill = bcs_colors["c"]),
    legend.text = element_text(color = bcs_colors["dg"]),
    legend.title = element_text(color = bcs_colors["dg"]),
    
    # Facets
    strip.background = element_rect(fill = bcs_colors["yg"]),
    strip.text = element_text(color = bcs_colors["dg"], face = "bold")
  )
}

# libraries
library(tidyverse)
library(unmarked)
library(readxl)

# load data----
nbp <- read_xlsx("data/b_intermediate_data/nbp_tidy_jan_24.xlsx")
circ <- read_csv("data/c_analysis_data/nbp_circ_codes.csv")
covs <- read_csv("data/c_analysis_data/circ_no_overlap_covariates.csv")

## associate each count circle with its code
circ$station_id <- paste(circ$park, circ$loop, circ$station, sep = "-")
nbp$station_id <- paste(nbp$park, nbp$loop, nbp$station, sep = "-")

nbp <- left_join(nbp, select(circ, station_id, code))

## make covs column names more consistent with other datasets
names(covs) <- str_to_lower(names(covs))
names(covs)[names(covs) == "station"] <- "code"
colnames(covs)

## standardize covs
covs$scalecanopyprop2016 <- scale(covs$canopyprop2016)
covs$scalecanopyprop2021 <- scale(covs$canopyprop2021)
covs$scalemdowntown <- scale(covs$mdowntown)
covs$scalemmajorwater <- scale(covs$mmajorwater)
view(covs)
colnames(covs)
# subset data such that we limit subsequent analysis to count circles that were
# surveyed in Jan, Feb, March of 2017 and 2022 and have site-level and observation-
# level covariates

t1 <- 2017
t2 <- 2022
months <- c(1, 2, 3)

non_na_covs <- covs %>%
  na.omit() %>%
  pull(code)


d <- nbp %>%
  filter(year %in% c(t1, t2),
         month %in% months,
         !is.na(weather),
         !is.na(dogs),
         code %in% non_na_covs)

whi <- d %>%
  group_by(year, code) %>%
  summarise(n = n_distinct(survey_id)) %>%
  pivot_wider(names_from = year, values_from = n)

#view(whi)

sus <- whi[whi[, 1] == 3 & whi[, 2] == 3, ]$code

whi[whi[, 2] == 3 & whi[, 3] ==3, ]$code

sus ## this should be a list of sample units that have all necessary data


#### SINGLE SEASON OCCUPANCY ######


sp <- "Red-breasted Nuthatch"  ## set your species of interest
t1 <- 2017 # set your baseline year
t2 <- 2022 # set a second year to compare against baseline
months <- c(1, 2, 3) # set months (visits)--will usually be Jan, Feb, March

d <- nbp %>%
  filter(code %in% sus,
         month %in% months,
         year %in% c(t1, t2))

srm_t1 <- d %>%   ## species response matrix for t1
  filter(year == t1) %>%
  group_by(code, month, species) %>%
  summarise(count = sum(seen, heard, fly)) %>%
  ungroup() %>%
  pivot_wider(names_from = species, values_from = count) %>%
  replace(is.na(.), 0) %>%
  mutate(month = month.name[month])

sp_t1 <- srm_t1 %>%  ## species detection matrix for t1
  select(code, month, all_of(sp)) %>%
  pivot_wider(names_from = month, values_from = sp) %>%
  mutate_if(is.numeric, ~1 * (. > 0))

srm_t2 <- d %>%   ## species response matrix for t2
  filter(year == t2) %>%
  group_by(code, month, species) %>%
  summarise(count = sum(seen, heard, fly)) %>%
  ungroup() %>%
  pivot_wider(names_from = species, values_from = count) %>%
  replace(is.na(.), 0) %>%
  mutate(month = month.name[month])

sp_t2 <- srm_t2 %>%  ## species detection matrix for t2
  select(code, month, all_of(sp)) %>%
  pivot_wider(names_from = month, values_from = sp) %>%
  mutate_if(is.numeric, ~1 * (. > 0))

# view(sp_t1)
# view(sp_t2)

y_t1 <- as.matrix(sp_t1[, -1])
y_t2 <- as.matrix(sp_t2[, -1])

unmarked_simple_t1 <- unmarkedFrameOccu(y = y_t1)
summary(unmarked_simple_t1) # 10/76 for DOWO

unmarked_simple_t2 <- unmarkedFrameOccu(y = y_t2)
summary(unmarked_simple_t2) # naive occupancy 15/76 for DOWO

occu.m.t1 <- occu(~1 ~1, data = unmarked_simple_t1)
occu.m.t2 <- occu(~1 ~1, data = unmarked_simple_t2)

predict(occu.m.t1, newdata = data.frame(site = 1), type = "state")
predict(occu.m.t2, newdata = data.frame(site = 1), type = "state")

predict(occu.m.t1, newdata = data.frame(site = 1), type = "det")
predict(occu.m.t2, newdata = data.frame(site = 1), type = "det")


## create observation covariates

obs_covs_t1 <- list()
cov_names <- c('weather', 'precipitation', 'wind', 'dogs', 'dogs_off_leash', 'walkers')

length(obs_covs_t1)

for(i in 1:length(cov_names)){
  covar <- cov_names[i]
  obs_covs_t1[[i]] <- nbp %>% select(code, all_of(covar), year, month) %>%
    filter(year == t1, month %in% months, code %in% sus, !duplicated(.)) %>%
    mutate(month = month.name[month]) %>%
    arrange(code) %>%
    pivot_wider(names_from = month, values_from = covar) %>%
    select(-year, -code)
}

names(obs_covs_t1) <- cov_names
obs_covs_t1


obs_covs_t2 <- list()

for(i in 1:length(cov_names)){
  covar <- cov_names[i]
  obs_covs_t2[[i]] <- nbp %>% select(code, all_of(covar), year, month) %>%
    filter(year == t2, month %in% months, code %in% sus, !duplicated(.)) %>%
    mutate(month = month.name[month]) %>%
    arrange(code) %>%
    pivot_wider(names_from = month, values_from = covar) %>%
    select(-year, -code)
}

names(obs_covs_t2) <- cov_names
obs_covs_t2

## Create circle covs
circ_covs_t1 <- covs %>%
  filter(code %in% sus) %>%
  select(park, scalecanopyprop2016, scalemdowntown, scalemmajorwater)

circ_covs_t2 <- covs %>%
  filter(code %in% sus) %>%
  select(park, scalecanopyprop2021, scalemdowntown, scalemmajorwater)

# build new frame
unmarked_cov_t1 <- unmarkedFrameOccu(y = y_t1, obsCovs = obs_covs_t1, siteCovs = circ_covs_t1)
summary(unmarked_cov_t1)

unmarked_cov_t2 <- unmarkedFrameOccu(y = y_t2, obsCovs = obs_covs_t2, siteCovs = circ_covs_t2)
summary(unmarked_cov_t2)


# build new occupancy models

# weather detection
occu.m2.t1 <- occu(formula = ~wind + weather + precipitation
                   ~1, 
                   data = unmarked_cov_t1)

# human-related detection 
occu.m3.t1 <- occu(formula = ~walkers + dogs + dogs_off_leash
                   ~1,
                   data = unmarked_cov_t1)

# environmental occupancy
occu.m4.t1 <- occu(formula = ~1 
                   ~canopyprop2016 + mdowntown + mmajorwater,
                   data = unmarked_cov_t1)

occu.m5.t1 <- occu(formula = ~1
                   ~scalecanopyprop2016,
                   data = unmarked_cov_t1)

summary(occu.m2.t1)
summary(occu.m3.t1)
summary(occu.m4.t1)
summary(occu.m5.t1)
