focal.parks <- c("Carkeek Park", "Discovery Park", "Genesee Park", "Golden Gardens Park", "Lake Forest Park", 
                 "Magnuson Park", "Seward Park", "Washington Park Arboretum")

nsurv <- dat %>% 
  filter(year %in% c(2005:2019, 2021, 2022, 2023),
         park %in% focal.parks) %>%
  group_by(year) %>%
  summarise(nsurv = n_distinct(survey_id), .groups = "drop")

srm <- dat %>% 
  filter(!str_detect(species, pattern = " x | sp\\.|Spotted Owl"),
         year %in% c(2005:2019, 2021, 2022, 2023),
         park %in% focal.parks) %>%
  mutate(observed = seen + heard + fly) %>% 
  group_by(year, bird.code) %>%
  summarise(count = sum(observed), .groups = "drop") %>%
  pivot_wider(names_from = bird.code, values_from = count, values_fill = 0) %>%
  pivot_longer(-c(1, 2), names_to = "bird.code", values_to = "count") %>%
  left_join(., nsurv, join_by("year" == "year")) %>%
  mutate(maps = count / nsurv)


# Create a time period column
srm$time_period <- ifelse(srm$year <= 2014, "early", "late")

# Calculate mean abundance for each species in each time period
mean_abundance <- srm %>%
  group_by(bird.code, time_period) %>%
  summarise(mean_maps = mean(maps, na.rm = TRUE))

# Spread the data to compare early vs late period
abundance_change <- mean_abundance %>%
  spread(key = time_period, value = mean_maps) %>%
  mutate(change = (late - early) / early)

# Check species with negative change (less frequent in the later period)
declining_species <- abundance_change %>%
  filter(change < 0) %>%
  arrange(change)

view(declining_species)
