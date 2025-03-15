#################################################
### 08. community response to land use change ###
#################################################


# Question: Do sites that have experienced measurable land use change differ
# in species composition pre/post change?

# objective: export species response matrix for sites with complete data for years 2019 and 2023. 
# data will be used for 

# libraries
library(tidyverse)
library(readxl)
library(openxlsx)

# load data
d <- read_xlsx("data/b_intermediate_data/nbp_tidy_jan_24.xlsx") 
circs <- read_csv("data/c_analysis_data/nbp_circ_codes.csv")
hrcd <- read_xls("data/zz_miscellaneous_data/NBP_countCircles_Points_withNearHRCD.xls")
traits <- read_xlsx("data/zz_miscellaneous_data/NBP_species_list+functional_traits.xlsx")

# add station codes to nbp data
d <- left_join(d, circs) %>% select(obs_id, year, month, code, everything())

# add alpha code to nbp data
d <- left_join(d, select(traits, com.name, alpha.code), join_by("species" == "com.name"))

# create a sorting order for codes
sort.order <- data.frame(code = circs$code,
                       sort.ord = seq(1:length(circs$code))) %>% na.exclude()

# find circles that contain land use change
change <- hrcd %>% filter(NEAR_DIST < 50 & NEAR_DIST >= 0) %>% pull(Station)


dat <- d %>%
  filter(year %in% c(2019, 2023), 
         code %in% change,
         !str_detect(species, pattern = " sp.| x "))

completeness <- dat %>%
  group_by(code, year) %>%
  reframe(nsurv = n_distinct(survey_id)) %>%
  pivot_wider(names_from = year, values_from = nsurv) %>%
  filter(`2019` == 12 & `2023` == 12)

view(completeness)

<<<<<<< HEAD
## not really enough. Maybe look at just a few months?
=======
## not sufficient. Maybe look at just spring when birds are most vocal?
>>>>>>> 7347b67d927aa10299344d068e014c2ddcbb1568

dat2 <- d %>%
  filter(year %in% c(2019, 2023),
         month %in% c(4, 5, 6),
         code %in% change,
         !str_detect(species, pattern = " sp.| x "))

completeness2 <- dat2 %>%
  group_by(code, year) %>%
  reframe(nsurv = n_distinct(survey_id)) %>%
  pivot_wider(names_from = year, values_from = nsurv) %>%
  filter(`2019` == 3 & `2023` == 3)

view(completeness2)

view(filter(hrcd, Station %in% completeness2$code))

change.codes <- completeness2 %>% filter(code != "LFM1") %>% pull(code)

# now need a sample of units with no chnage

dat3 <- d %>%
  filter(year %in% c(2019, 2023),
         month %in% c(4, 5, 6),
         code %in% hrcd$Station,
         !str_detect(species, pattern = " sp.| x "))

completeness3 <- dat3 %>%
  group_by(code, year) %>%
  reframe(nsurv = n_distinct(survey_id)) %>%
  pivot_wider(names_from = year, values_from = nsurv) %>%
  filter(`2019` == 3 & `2023` == 3,
         !code %in% completeness2$code)

no.change.codes <- completeness3[sample(dim(completeness3)[1], 12), ]$code

sus <- c(change.codes, no.change.codes)

# this is a better list

df <- d %>%
  filter(year %in% c(2019, 2023),
         month %in% c(4, 5, 6),
         !str_detect(species, pattern = " sp.| x "),
<<<<<<< HEAD
         code %in% sus,
         seen > 0 | heard > 0)
=======
         code %in% completeness2$code,
         code != "LFM1") ## lake forest park is an oddball, exclude

>>>>>>> 7347b67d927aa10299344d068e014c2ddcbb1568

spp <- sort(unique(df$alpha.code))

srm <- df %>%
  group_by(code, year, alpha.code) %>%
  summarise(count = sum(seen, heard)) %>%
  pivot_wider(names_from = alpha.code, values_from = count) %>%
  replace(is.na(.), 0) %>%
  ungroup() %>%
  left_join(., sort.order) %>%
  arrange(year, sort.ord) %>%
  mutate(time = case_when(year == 2019 ~ 0,
                          TRUE ~ 1),
         change = case_when(code %in% change.codes ~ 1,
                            TRUE ~ 0),
         site = paste(code, change, time, sep = ".")) %>%
  select(site, all_of(spp))

view(srm)
summary(srm)

# ok I think this looks good.

# now filter hrcd data for export

hrcd.d <- hrcd %>%
  filter(Station %in% df$code) %>%
  select(Station, Park, DeltaCanopy, NEAR_DIST, NEAR_CHANGE)

hrcd.df <- hrcd.d %>%
  rbind(., hrcd.d) %>%
  mutate(time = c(rep(0, dim(hrcd.d)[1]), rep(1, dim(hrcd.d)[1])),
         change = case_when(Station %in% change.codes ~ 1,
                            TRUE ~ 0),
         site = paste(Station, change, time, sep = ".")) %>%
  select(site, time, change, NEAR_DIST, DeltaCanopy)

sort(hrcd.d$Station)
sort(unique(df$code))
view(hrcd.df)

## good--export files
write.xlsx(srm, "data/c_analysis_data/srm_2019&2023_hrcd.xlsx")
write.xlsx(hrcd.df, "data/c_analysis_data/covs_2019&2023_hrcd.xlsx")


## so the results look complcated by the inclusion of wetland/open water habitats...


