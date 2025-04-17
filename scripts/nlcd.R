library(tidyverse)
nlcd <- read.csv("data/raw/NLCD_2005-2023.csv")
colnames(nlcd)


long <- nlcd %>% pivot_longer(-c(1, 2), names_to = "VALUES", values_to = "COUNTS")

head(long)
xwalk <- data.frame(VALUES = c("VALUE_11", "VALUE_21", "VALUE_22", "VALUE_23", "VALUE_24",
                     "VALUE_31", "VALUE_41", "VALUE_42", "VALUE_43", "VALUE_52",
                     "VALUE_71", "VALUE_81", "VALUE_90", "VALUE_95"),
           CLASS = c("Open Water", "Developed-Open Space", "Developed-Low Intensity",
                     "Developed-Medium Intensity", "Developed-High Intensity",
                     "Barren Land", "Deciduous Forest", "Evergreen Forest",
                     "Mixed Forest", "Shurb/Scrub", "Grassland/Herbaceous", 
                     "Pasture/Hay", "Woody Wetlands", "Emergent Herbaceous Wetlands"))

long <- left_join(long, xwalk) %>% mutate(PARK.VALUES = paste(PARK, VALUES, sep = "-"))

ggplot(long %>% filter(VALUES == "VALUE_95"), aes(x = YEAR, y = COUNTS)) +
  geom_line() +
  facet_wrap(~ PARK, scales = "free")



view(long)



study_areas <- data.frame(PARK = c("Carkeek Park", "Cheasty Greenspace", "Discovery Park",
                     "Genesee Park", "Golden Gardens", "Lake Forest Park", 
                     "Lincoln Park", "Magnuson Park", "Seward Park", "Washington Park Arboretum"),
                     AREA = c(20327246.333613, 11087744.649455, 49606576.253013, 8422232.510264,
                     11970831.954355, 4713976.629057, 14745041.309755, 25244983.226958,
                     20260635.330094, 12906043.486227))


t <- long %>% group_by(PARK, YEAR) %>% summarise(total = round(sum(COUNTS),0), .groups = "drop") %>% select(PARK, total) %>% unique()


long.t <- left_join(long, t) %>% mutate(perc = COUNTS / total)

pdat <- long.t %>% filter(PARK == "Discovery Park" & perc > .001 & CLASS != "Open Water") %>% 
  mutate(CLASS2 = fct_reorder(CLASS, desc(perc)))

p <- ggplot(pdat, aes(x = YEAR, y = perc, color = CLASS2)) +
         geom_line(linewidth = 0.5, alpha = 0.5) +  geom_point(size = 1, alpha = 0.5) + 
  labs(title = "Land Cover Composition over Time",
       subtitle = "Discovery Park", x = "Year" , y = "% land cover within 300 meters of count stations",
       color = "Land cover class") +
  scale_color_manual(values = unname(chart_colors)) +
  theme_bcs()+
  theme(legend.key = element_rect(color = NA, fill = NA),
                      legend.key.size = unit(0.2, "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "in"),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "in"))

png(filename = "results/figures/nlcd_discovery.png", height = 3, width = 3, units = "in", res = 300)                           
p
dev.off()


area.richness <- res %>% left_join(., stations, join_by("park" == "park"))


ggplot(area.richness, aes(x = log(n.station), y = jack1, color = park)) +
  geom_point() + 
  theme_bcs() + 
  labs(title = "Estimated species richess vs. area surveyed", x = "No. of count stations", y = "Estimate")+
  scale_color_manual(name = "Park",
                     values = unname(chart_colors)) +
  scale_y_continuous(limits = c(0, 210))

log(0)
