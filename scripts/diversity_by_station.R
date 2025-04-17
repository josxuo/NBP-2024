#########################
#### DIVERSITY INDEX ####
#########################


library(vegan)
nbp <- read_excel("data/processed/nbp_tidy_jan_24.xlsx")

nbp[nbp$species == "Glaucous-winged x Western Gull", ]$bird.code <- "GWGU"
nbp[nbp$species == "American x Eurasian Wigeon", ]$bird.code <- "AMWI"

covs <- read.csv("data/processed/circ_no_overlap_covariates.csv")

focal.parks <- c("Carkeek Park", "Cheasty Greenspace", "Discovery Park", "Genesee Park", 
                 "Golden Gardens Park", "Lincoln Park", "Magnuson Park", "Seward Park", "Washington Park Arboretum")


library(vegan)

diversity_by_station <- nbp %>% 
  filter(park %in% focal.parks,
         !str_detect(species, pattern = "sp\\.|Domestic|Spotted Owl")) %>%
  group_by(station.code, bird.code) %>%
  summarise(count = mean(seen + fly + heard), .groups = "drop") %>%
  pivot_wider(names_from = bird.code, values_from = count, values_fill = 0) %>%
  rowwise() %>%
  mutate(shannon = diversity(c_across(-station.code), index = "shannon")) %>%
  ungroup()


write.csv(diversity_by_station, "diversity_by_station.csv", row.names = FALSE)
