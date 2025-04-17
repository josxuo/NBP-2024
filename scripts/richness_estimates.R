library(vegan)

covs <- read.csv("data/processed/circ_no_overlap_covariates.csv")

focal.parks <- c("Carkeek Park", "Cheasty Greenspace", "Discovery Park", "Genesee Park", 
                 "Golden Gardens Park", "Lincoln Park", "Magnuson Park", "Seward Park", "Washington Park Arboretum")


srm <- nbp %>% 
  filter(!str_detect(species, pattern = "sp\\.|Spotted Owl|Domestic")) %>%
         #year %in% c(2005:2019, 2021, 2022, 2023),
         #station.code %in% covs$station.code) %>%
  mutate(observed = seen + heard + fly) %>% 
  group_by(park, year, station.code, bird.code) %>%
  summarise(count = sum(observed), .groups = "drop") %>%
  pivot_wider(names_from = bird.code, values_from = count, values_fill = 0)

srm.focal <- srm %>% filter(park %in% focal.parks)

dat.list <- list()

for(i in 1:length(focal.parks)){
    
    subset_data <- srm.focal %>% filter(park == focal.parks[i]) %>% select(-park, -year, -station.code)
    
    dat.list[[paste0(focal.parks[i])]] <- subset_data
  }

richness.estimates <- list()

for(i in 1:length(dat.list)){
  est <- specpool(dat.list[[i]])
  richness.estimates[[names(dat.list[i])]] <- est
}

combined_estimates <- bind_rows(richness.estimates)
head(combined_estimates)

combined_estimates$park <- names(dat.list)

combined_estimates.df <- combined_estimates %>%
  separate(park.year, into = c("park", "year"), sep = "-") %>%
  mutate(year = as.numeric(year))

plot.df <- combined_estimates.df %>%
  select(-grep("\\.se|^n$|Lower|Upper", colnames(.))) %>%
  pivot_longer(cols = c(boot, chao, jack1, jack2, Species),
               names_to = "estimator",
               values_to = "richness_estimate")

plot.df[plot.df$estimator == "Species", ]$estimator <- "Observed richness"
plot.df[plot.df$estimator == "boot", ]$estimator <- "Bootstrapped estimate"
plot.df[plot.df$estimator == "chao", ]$estimator <- "Chao1 estimate"
plot.df[plot.df$estimator == "jack1", ]$estimator <- "First-order jackknife estimate"
plot.df[plot.df$estimator == "jack2", ]$estimator <- "Second-order jackknife estimate"

plot.df$estimator <- factor(plot.df$estimator, levels = c("Observed richness", "Chao1 estimate", 
                                                          "First-order jackknife estimate", "Second-order jackknife estimate",
                                                          "Bootstrapped estimate"))

# panel plot of species richness measures over time at focal parks

p.park_ests_over_time <- ggplot() +
  geom_col(data = plot.df %>% filter(estimator == "Observed richness"), 
           aes(x = year, y = richness_estimate, group = interaction(park, estimator), 
               fill = estimator)) +
  geom_line(data = plot.df %>% filter(estimator != "Observed richness"), 
            aes(x = year, y = richness_estimate, group = interaction(park, estimator), 
                color = estimator), size = 0.5, alpha = 0.8) +
  scale_color_manual(name = "Measure",
                     values = c("Chao1 estimate" = "#36BA3A", 
                                "First-order jackknife estimate" = "#FFB98C", 
                                "Second-order jackknife estimate" = "#0A3C23",
                                "Bootstrapped estimate"= "#E6FF55")) +
  scale_fill_manual(name = element_blank(),
                    values = c("Observed richness" = "#0A3C23")) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0)),
                     limits = c(0, 180)) +
  facet_wrap(~ park, ncol = 4, nrow = 2) +
  labs(title = "Estimates of species richness over time", y = "Number of species", x = "Year") +
  theme_bcs() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.spacing = unit(-1, "lines"),
        legend.key.size = unit(0.8, "lines"))
