#############################
### 02. nbp summary stats ###
#############################

# clear environment
rm(list = ls())

# load packages
library(tidyverse)
library(readxl)

# Custom bcs plotting theme 
source("functions/theme_bcs.R")

# load tidied data (see script "01_tidy_raw_nbp.R" for details)
dat <- read_excel("data/processed/nbp_tidy_jan_24.xlsx")
dat[dat$species == "Glaucous-winged x Western Gull", ]$bird.code <- "GWGU"
dat[dat$species == "American x Eurasian Wigeon", ]$bird.code <- "AMWI"

# inspect data if desired
## str(dat)
## head(dat)


#############################
##### data completeness #####
#############################

# How complete are the data? do we have surveys for each loop for each month that we expect we would?
## What time range does the dataset cover?
range(dat$survey_date) # 1996-04-13 to 2024-01-20, so the dataset should include 8 months for 1996;
                       # 12 months for 1997-2019; in 2020 we shut down NBP in march and 1 month for 2024 (noting that we shut down NBP in 2020)


## Create dataframe with year and number of months NBP was active in that year
dmths <- data.frame(year = sort(unique(dat$year)),
                    Nmonth = c(8, rep(12, 27), 1))

## Join months data frame with a summary table of how many loops each park had per year, and how many surveys
## were done for each loop at each park each year, then create a "completeness" column
dcmp <- dat %>%
  group_by(year, park) %>%
  summarise(Nloop = n_distinct(loop), nloop = n_distinct(paste0(survey_date, loop))) %>%
  ungroup() %>%
  left_join(., dmths) %>%
  mutate(completeness = nloop / Nloop / Nmonth)

## some greater than 1. 
dcmp$completeness[dcmp$completeness > 1] <- 1

## visualize with a heat map
pcmp <- ggplot(dcmp, aes(x = year, y = park, fill= completeness)) + 
  geom_tile() + 
  scale_fill_gradient(low = bcs_colors["yellow green"], high = bcs_colors["dark green"]) + 
  labs(x = "", y = "") +
  ggtitle("Dataset completeness") + 
  theme_bcs() 

pcmp

## print plot if desired
# pdf("figures/02_data_completeness.pdf", width = 11)
# print(pcmp)
# dev.off()


 #############################
 ####### survey summary ######
 #############################
 
## HOW MANY SURVEYS IN DATASET?

length(unique(dat$survey_id))  # 39174 surveys BUT this includes the odd Magnuson surveys with non-existent station IDs
dat %>% filter(!is.na(station.code)) %>% summarise(nsurv = n_distinct(survey_id))  ## 38,771 surveys

## SURVEYS BY YEAR
### create dataframe with count of distinct surveys per year, then plot results
dyr <- dat %>%
  filter(!is.na(station.code)) %>%
  group_by(year) %>%
  summarise(nsurv = n_distinct(survey_id))

pyr <- ggplot(dyr, aes(x = year, y = nsurv)) +
  geom_col(fill = bcs_colors["dark green"]) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  xlab("Year") + ylab("Number of NBP point counts in dataset") + ggtitle("NBP point counts by year") +
  theme_bcs() +
  theme(panel.grid.major.x = element_blank())

pyr

### print plot to pdf if desired
# pdf("figures/02_nsurv_yr.pdf")
# print(pyr)
# dev.off()

## SURVEYS BY LOCATION
### create dataframe with count of distinct surveys by park, then plot results
dprk <- dat %>%
  filter(!is.na(station.code)) %>%
  group_by(park) %>%
  summarise(nsurv = n_distinct(survey_id)) %>%
  ungroup() %>%
  mutate(park = fct_reorder(park, nsurv))

pprk <- ggplot(dprk, aes(x = park, y = nsurv)) +
  geom_col(fill = bcs_colors["dark green"]) +
  coord_flip() +
  xlab("") + ylab("Number of NBP point counts in dataset") + ggtitle("NBP point counts by park") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_bcs() +
  theme(panel.grid.major.y = element_blank())

pprk

#### print chart to pdf if desired
# pdf("figures/02_nsurv_park.pdf", width = 11, height = 7)
# print(pprk)
# dev.off()

#############################
###### species summary ######
#############################

# SPECIES TOTAL
## remove spurious observations, observations not resolved to species level and remove hybrids
dsp <- dat %>% filter(!str_detect(species, pattern = " sp\\.|Spotted Owl| x "))
sort(unique(dsp$species)) # 211 species

## FREQUENTLY REPORTED SPECIES
### step 1, create species response matrix for each unique survey
dfreq <- dat %>% 
  mutate(observed = 1) %>%
  select(survey_id, species, observed) %>%
  pivot_wider(names_from = species, values_from = observed, values_fill = 0)

### step 2: create mean abundance per survey 
d <- colnames(dfreq[, -1]) ## vector of the species names
times <- colSums(dfreq[, -1], na.rm = TRUE) ## vector of the number of times each species has been counted
times <- enframe(times, name = "species", value = "nobs")

#### add names and number of reports to a data frame, then create a new column dividing number of surveys reporting by total surveys
times$prop <- times$nobs / dim(dfreq)[1]

#### reorder species factor levels to make plot more appealing
times$species <- fct_reorder(times$species, times$prop)

#### plot top 20 most frequenlty reported species
pprop <- ggplot(slice_max(times, prop, n = 20), aes(x = species, y = prop)) +
  geom_col(fill = bcs_colors["dark green"]) +
  coord_flip() + 
  xlab("") + ylab("Propotion of surveys reporting") + ggtitle("Most frequently reported species on NBP counts") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_bcs() +
  theme(panel.grid.major.y = element_blank())

pprop

#### print plot if desired
# pdf("figures/02_prop_top_20.pdf", width = 11, height = 7)
# print(pprop)
# dev.off()

# NOT ID'D TO SPECIES LEVEL
## create dataframe with just unresolved species, then count how many times each sp. occurs in the dataset
dunr <- dat %>% filter(str_detect(species, pattern = " sp\\.")) %>%
  group_by(species) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(species = fct_reorder(species, desc(n)))


punr <- ggplot(slice_max(dunr, order_by = n, n = 15), aes(x = species, y = n)) + 
  geom_col(fill = bcs_colors["dark green"]) +
  labs(title = "NBP observations not identified to species level", y = "Number of times reported", x = "") + 
  scale_y_continuous(breaks = seq(0, 1400, by = 50), expand = c(0, .2), limits = c(0, 1425)) +
  theme_bcs() +
  theme(
    #panel.grid.major.y = element_blank(),
    axis.text.x = element_text(angle = 60, 
                               hjust = 1, 
                               margin = margin(t = 2, r = 0, b = 0, l = 0, unit = "pt")))

punr


png("results/figures/02_spuh_reports.png", width = 4, height = 4.5, unit = "in", res = 300)
punr
dev.off()

### print plot
# pdf("figures/02_sp._reports.pdf", width = 11, height = 7)
# print(punr)
# dev.off()

# ABUNDANCE

## Which species have the highest mean abundance per survey (maps)
dabund <- dat %>% filter(!str_detect(species, pattern = "sp\\.| x ")) %>%
  group_by(species) %>%
  summarise(nobs = sum(seen, heard, fly), nsurv = n_distinct(survey_id), maps_n = nobs / nsurv,
            maps_N = nobs / length(unique(dat$survey_id))) %>%
  ungroup() %>%
  mutate(species_n = fct_reorder(species, maps_n), species_N = fct_reorder(species, maps_N))

pmaps_n <- ggplot(slice_max(dabund, maps_n, n = 20), aes(x = species_n, y = maps_n)) +
  geom_col(fill = bcs_colors["dark green"]) +
  coord_flip() + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "", title = "Mean observed flock size", y = "Number of individuals") +
  theme_bcs() +
  theme(panel.grid.major.y = element_blank())

pmaps_n  ## this stat is driven by bird group size, but species may be very infrequently observed (e.g., western bluebird)

pmaps_N <- ggplot(slice_max(dabund, maps_N, n = 20), aes(x = species_N, y = maps_N)) +
  geom_col(fill = bcs_colors["dark green"]) +
  coord_flip() + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "", y = "Mean abundance per survey", title = "Mean abundance per survey") +
  theme_bcs() +
  theme(panel.grid.major.y = element_blank())

pmaps_N  ## this stat is largely driven by frequency of reports

# Visualizing mean abundance per survey by species across years
## use mean abundance per survey N, stat mostly driven by observation frequency rather than group size
### new dataframe calculating number of surveys conducted each year
dNsurvyr <- dat %>%
  group_by(year) %>%
  summarise(Nsurv = n_distinct(survey_id))

### new data frame with mean abundance per survey per year
dmapsN <- dat %>%
  group_by(year, species) %>%
  summarise(nobs = sum(seen, fly, heard), .groups = "drop") %>%
  left_join(., dNsurvyr) %>%
  mutate(maps_N = nobs / Nsurv) %>%
  select(year, species, maps_N) %>%
  pivot_wider(names_from = species, values_from = maps_N, values_fill = 0) %>%
  pivot_longer(-1, names_to = "species", values_to = "maps_N")

### Visualize on heat map with sets of random species or by key word

spp <- d[sample(1:length(d), 20)]  # this will return 20 random
# spp <- d[str_detect(d, "Warbler")]  # choose a species or group of birds that share a common work (e.g, warbler, thrush, swallow)
# spp <- c("CHOOSE YOUR OWN")

ggplot(filter(dmapsN, species %in% spp), aes(x = year, y = species, fill = maps_N)) +
  geom_tile() +
  scale_fill_gradient(low = bcs_colors["yellow green"], high = bcs_colors["dark green"]) + 
  labs(x = "Year", y = "", fill = str_wrap("Mean abundance per survey", width = 10)) +
  ggtitle("Mean abundance per survey") + 
  theme_bcs()

# SPECIES TOTALS BY LOCATION
dS <- dat %>%
  filter(!str_detect(species, pattern = " sp\\.| x ",),  ## filter out sp and hybrids
         species != "Spotted Owl") %>%  ## truely don't think we saw the spotted owl at Magnuson
  group_by(park) %>% 
  summarise(S = n_distinct(species)) %>%
  ungroup() %>%
#  left_join(., active) %>%    ## could create a join to label parks by their current NBP status (active vs inactive)
#  replace(is.na(.), "Not Active") %>% 
  mutate(park = fct_reorder(park, S))

meanS <- mean(dS$S)  # average species reported by park

pmeanS <- ggplot(dS, aes(x = park, y = S)) +
  geom_col(fill = bcs_colors["dark green"]) + 
  ylab("Total species reported") + xlab("") + ggtitle("NBP total species reported") +
  geom_hline(yintercept = mean(dS$S), color = bcs_colors["bright green"], lty = 2) +
  geom_text(aes(label = S), color = bcs_colors["cream"], hjust = 1.4, family = "Archivo", size = 5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)),
                     breaks = c(75, 150)) + 
  coord_flip() + 
  annotate("text", x = 1, y = (meanS + 10), label = paste("Avg =", round(meanS, 0), "species"), 
           color = bcs_colors["dark green"], family = "Archivo", size = 5) +
  theme_bcs() +
  theme(panel.grid.major.y = element_blank())

pmeanS

## print plot if desired
# pdf("figures/02_mean_S_by_park.pdf", width = 11)
# print(pmeanS)
# dev.off()


library(vegan)
nbp <- dat
dat <- nbp


######### AGGREGATE SPECIES RICHNESS TRENDS ################
srm <- dat %>% 
  filter(!str_detect(species, pattern = "sp\\.|Spotted Owl|Domestic"),
         year %in% c(2005:2019, 2021, 2022, 2023), 
         station.code %in% covs$station.code) %>%
  mutate(observed = seen + heard + fly) %>% 
  group_by(park, year, station.code, bird.code) %>%
  summarise(count = sum(observed), .groups = "drop") %>%
  pivot_wider(names_from = bird.code, values_from = count, values_fill = 0) %>%
  mutate(station.year = paste(station.code, year, sep = "-"))


nsurv <- dat %>% 
  filter(!str_detect(species, pattern = "sp\\.|Spotted Owl|Domestic"),
         year %in% c(2005:2019, 2021, 2022, 2023)) %>%
  group_by(year, station.code) %>%
  summarise(nsurv = n_distinct(survey_id), .groups = "drop") %>% 
  mutate(station.year = paste(station.code, year, sep = "-"))


years <- unique(srm$year)

res.df <- data.frame(year = numeric(),
                     Species = numeric(),
                     jack1 = numeric())

for(i in 1:length(years)){
  year.dat <- srm %>% filter(year == years[i]) %>% select(-park, -station.code, -year,-station.year)
  res <- specpool(year.dat)
  res.df[i, 1] <- years[i]
  res.df[i, 2] <- res$Species
  res.df[i, 3] <- res$jack1
}

res.df <- left_join(res.df, nsurv.disco)

head(res.df)

plot(res.df$year, res.df$jack1)

combined <- ggplot(res.df, aes(x = year, y = jack1)) +
  geom_point(color = bcs_colors["dark green"]) + 
  geom_smooth(method = "lm", color = bcs_colors["peach"], linetype = "dashed") +
  labs(title = "Estimated combined species richness", subtitle = "All NBP parks", x = "Year", y = "Estimated species richness") +
  theme_bcs()


png(filename = "results/figures/combined_richness_over_time.png",
    height = 3, width = 3, units = "in", res = 300)
combined
dev.off()

covs <- read.csv("data/processed/circ_no_overlap_covariates.csv")

focal.parks <- c("Carkeek Park", "Discovery Park", "Genesee Park", "Golden Gardens Park", 
                 "Magnuson Park", "Seward Park", "Cheasty Greenspace", "Lake Forest Park", "Washington Park Arboretum")

srm <- dat %>% 
  filter(!str_detect(species, pattern = "sp\\.|Spotted Owl|Domestic"),
         year %in% c(2005:2019, 2021, 2022, 2023)) %>%
         #station.code %in% covs$station.code,
         #park == "Discovery Park", 
         #!loop %in% c("Capehart", "Cemetery", "Nike/500")) %>%
  mutate(observed = seen + heard + fly) %>% 
  group_by(park, year, station.code, bird.code) %>%
  summarise(count = sum(observed), .groups = "drop") %>%
  pivot_wider(names_from = bird.code, values_from = count, values_fill = 0)


nsurv.disco <- dat %>% 
  filter(!str_detect(species, pattern = "sp\\.|Spotted Owl|Domestic"),
         year %in% c(2005:2019, 2021, 2022, 2023),
         #station.code %in% covs$station.code,
         park == "Discovery Park", 
         !loop %in% c("Capehart", "Cemetery", "Nike/500")) %>%
  group_by(year) %>%
  summarise(nsurv = n_distinct(survey_id), .groups = "drop")


unique(nbp[nbp$park == "Discovery Park",]$loop)
disco.est <- specpool(select(srm, -park, -year, -station.code))

colnames(srm)

srm.focal <- srm %>% filter(park %in% focal.parks)

years <- unique(srm$year)

res.df <- data.frame(year = numeric(),
                     Species = numeric(),
                     jack1 = numeric())

for(i in 1:length(years)){
  year.dat <- srm %>% filter(year == years[i]) %>% select(-park, -station.code)
  res <- specpool(year.dat)
  res.df[i, 1] <- years[i]
  res.df[i, 2] <- res$Species
  res.df[i, 3] <- res$jack1
}

res.df <- left_join(res.df, nsurv.disco)

nbp %>% filter(park == "Discovery Park") %>% group_by(loop) %>%
  summarise(min = min(survey_date), max = max(survey_date))


disco.s.yr <- ggplot(data = res.df, aes(x = year, y = Species)) +
  #geom_line(aes(x = year, y = jack1, size = 0.5, alpha = 0.8)) +
  geom_point(color = "#0A3C23") +
  geom_smooth(method = "lm", color = "#FFB98C", linetype = "dashed") +
  #scale_color_manual(name = "Measure",
  #                   values = c("First-order jackknife estimate" = "#FFB98C")) +
  #scale_y_continuous(expand = expansion(mult = c(0, 0)),
  #                   limits = c(0, 180)) +
  labs(title = "Total species reported over time", subtitle = "Discovery Park",
       y = "Number of species", x = "Year", color = "Key") +
  theme_bcs() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.spacing = unit(-1, "lines"),
        legend.key.size = unit(0.8, "lines"))

png(filename = "results/figures/disco_S_over_time.png",
    height = 3, width = 3, units = "in", res = 300)
disco.s.yr
dev.off()


res.df$yr <- res.df$year - min(res.df$year)

mod <- lm(Species ~ yr, data = res.df)


summary(mod)


specpool(select(srm, -park, -station.code))

res <- data.frame(park = character(),
                  jack1 = numeric())

for(i in 1:length(focal.parks)) {
  rich.dat <- srm %>% filter(park == focal.parks[i]) %>% select(-park, -station.code)
  res[i, 1] <- focal.parks[i]
  res[i, 2] <- specpool(rich.dat)$jack1
}

view(res)

stations <- covs%>%group_by(park) %>%
  summarise(n.station = n_distinct(station.code))



year.list <- list()
for(i in 1:length(focal.parks)){
  year.list[[i]] <- sort(unique(srm.focal[srm.focal$park == focal.parks[i], ]$year))
}

names(year.list) <- focal.parks

dat.list <- list()

for(i in 1:length(focal.parks)){
  
  for(y in year.list[[focal.parks[i]]]) {
    
    subset_data <- srm.focal %>% filter(park == focal.parks[i] & year == y) %>% select(-park, -year, -station.code)
    
    dat.list[[paste0(focal.parks[i], "-", y)]] <- subset_data
  }
}

richness.estimates <- list()

for(i in 1:length(dat.list)){
  est <- specpool(dat.list[[i]])
  richness.estimates[[names(dat.list[i])]] <- est
}

combined_estimates <- bind_rows(richness.estimates)
head(combined_estimates)

combined_estimates$park.year <- names(dat.list)

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

#png(filename = "results/figures/02_estimates_spp_richness_over_time_by_park.png",
#    width = 7, height = 4, units = "in", res = 300)
#p.park_ests_over_time
#dev.off()


dplot <- plot.df %>% filter(park == "Discovery Park", 
                            (estimator == "Observed richness" | estimator == "First-order jackknife estimate"))

d_plot <- ggplot() +
  geom_col(data = dplot %>% filter(estimator == "Observed richness"), 
           aes(x = year, y = richness_estimate, group = interaction(park, estimator), 
               fill = estimator)) +
  
  geom_smooth(data = dplot %>% filter(estimator == "First-order jackknife estimate"), 
              aes(x = year, y = richness_estimate, linetype = "Best fit trend"), 
              color = "#FFB98C", linewidth = 0.5, alpha = 0.7, method = "lm", se = FALSE) +
  
  geom_line(data = dplot %>% filter(estimator == "First-order jackknife estimate"), 
            aes(x = year, y = richness_estimate, group = interaction(park, estimator), 
                color = estimator), size = 0.7) +
  
  scale_fill_manual(name = "Measure",
                    values = c("Observed richness" = "#0A3C23"),
                    limits = c("Observed richness")) + 
  scale_color_manual(name = element_blank(),
                     values = c(#"Chao1 estimate" = "#FFB98C", 
                                "First-order jackknife estimate" = "#36BA3A" 
                                #"Second-order jackknife estimate" = "#0A3C23",
                                #"Bootstrapped estimate"= "#E6FF55",
                                #"Best fit trend" = "black"),
                     ),
                     limits = c(#"Chao1 estimate", 
                                "First-order jackknife estimate"
                                #"Second-order jackknife estimate",
                               #"Bootstrapped estimate",
                                #"Best fit trend"))
                               )) +
  
  scale_linetype_manual(name = element_blank(),
                        values = c("Best fit trend" = "solid"),
                        limits = c("Best fit trend")) + # This keeps it last
  guides(#fill = guide_legend(order = 1),
           color = guide_legend(order = 1), 
         linetype = guide_legend(order = 2)) + # Ensures linetype is last
  
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)))+ 
  labs(title = "Estimates of Species Richness over Time", 
       subtitle = "Discovery Park", y = "Number of Species", x = "Year") +
  theme_bcs() +
  theme(legend.spacing = unit(-1, "lines"),
        legend.key.size = unit(0.2, "cm"), 
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "in"))



png(filename = "results/figures/discovery_species_over_time.png", height = 3.5, 
    width = 3.5, units = "in", res = 300)
d_plot
dev.off()

moddat <- filter(dplot, estimator == "First-order jackknife estimate") %>% 
  mutate(yr = as.numeric(scale(year)))

mod.null <- glm(round(richness_estimate, 0) ~ 1, data = moddat, family = poisson)
mod.1 <- glm(round(richness_estimate, 0) ~ yr, data = moddat, family = poisson) 

summary(mod.null)
summary(mod.1)

library(DHARMa)

plot(simulateResiduals(mod.1))

mod_qp <- glm(round(jack1, 0) ~ year + park, family = quasipoisson, data = combined_estimates.df)

summary(mod.null)
summary(mod.1)
summary(mod_qp)



diversity_by_year <- dat %>% 
  filter(park == "Discovery Park",
         !str_detect(species, pattern = "sp\\.|Domestic|Spotted Owl")) %>%
  group_by(station.code, year, bird.code) %>%
  summarise(count = mean(seen + fly + heard), .groups = "drop") %>%
  pivot_wider(names_from = bird.code, values_from = count, values_fill = 0) %>%
  rowwise() %>%
  mutate(shannon = diversity(c_across(-c(year, station.code)), index = "shannon")) %>%
  ungroup()

view(diversity_by_year[diversity_by_year$shannon == 0, ])

pdat <- diversity_by_year %>% filter(!year %in% c(2020, 2024)) %>% mutate(yr = as.numeric(scale(year)),
                                                                          shannon2 = ifelse(shannon == 0, 0.01, shannon)) %>%
  dplyr::select(year, yr, station.code, shannon, shannon2)

ggplot(pdat, aes(x = year, y = shannon)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)
  theme_bcs()

summary(pdat)
class(pdat$month)
hist(pdat$shannon)
shapiro.test(log(pdat$shannon + 1))
## data not normally distributed

pdat %>% group_by(year) %>% mutate(variance = var(shannon)) %>%
  ggplot(., aes(y = variance, group = year))+ geom_boxplot()

## variance not constant (but better aggregagted to year

mod <- glmmTMB(shannon2 ~ yr + station.code, family = tweedie(link = "log"), data = pdat)


summary(mod)


resids <- simulateResiduals(mod)
plot(resids)



mod <- glmmTMB((() ~ poly(yrmodmod <- glmmTMB((() ~ poly(yr,2) + as.factor(month) + (1 | station.code), zi = ~ as.factor(month), family = tweedie, data = pdat)

testZeroInflation(mod)
residuals <- simulateResiduals(mod)
plot(residuals)


pdat[pdglmmTMB()pdat[pdat$shannon == 0, ]

ggplot(pdat, aes(x = year, y = shannon)) + geom_point()

shapiro.test(pdat$shannon)

hist(log(pdat$shannon + 1))

hist(pdat$shannon)

summary(pdat$year)

library(MASS)
library(glmmTMB)

mod.gamma.null <- glmmTMB(shannon ~ 1, family = Gamma(link = "log"), data = pdat)
mod.gamma.1 <- glmmTMB(shannon ~ yr, family = Gamma(link = "log"), data = pdat)
mod.gamma.2 <- glmmTMB(shannon ~ yr + as.factor(month), family = Gamma(link = "log"), data = pdat)
mod.tweedie.null <- glmmTMB(shannon ~ 1, family = tweedie, data = pdat)
mod.tweedie.1 <- glmmTMB(shannon ~ yr, family = tweedie, data = pdat)
mod.tweedie.2 <- glmmTMB(shannon ~ yr + as.factor(month), family = tweedie, data = pdat)

mod.gaus.2 <- glmmTMB((1 / shannon) ~ yr + (1 | month), family = gaussian, data = pdat)

AIC(mod.gamma.null, mod.gamma.1, mod.gamma.2, mod.tweedie.null, mod.tweedie.1, mod.tweedie.2, mod.gaus, mod.gaus.2)

library(statmod)



mod.tweedie.2.1 <- glm( ~ yr + as.factor(month), family = tweedie(var.power = 6, link.power = 0), data = pdat)


par(mfrow = c(2, 2))
plot(mod.gaus)

predicted <- predict(mod.gaus.2, type = "response")

# Plot observed vs. predicted
plot((1/ pdat$shannon), predicted, main = "Observed vs. Predicted", 
     xlab = "Observed", ylab = "Predicted")
abline(a = 0, b = 1, col = "red")

residuals_gaus <- residuals(mod.gaus.2)
fitted_gaus <- fitted(mod.gaus.2)

qqnorm(residuals_gaus)
qqline(residuals_gaus, col = "red")

# Plot residuals vs. fitted values
plot(fitted_gaus, residuals_gaus, 
     xlab = "Fitted values", 
     ylab = "Residuals", 
     main = "Residuals vs Fitted Values")
abline(h = 0, col = "red")  # Add a horizontal line at 0 for reference

resids <- simulateResiduals(mod.gaus.2)
plot(resids)

mod.gaus <- glm(shannon ~ yr + as.factor(month), family = gaussian, data = pdat)


boxcox_res <- boxcox(mod.gaus)

plot(boxcox_res)

lambda <- boxcox_res$x[which.max(boxcox_res$y)]
pdat$shannon_bc <- (pdat$shannon ^ lambda - 1) / lambda  # Apply Box-Cox


# Fit the model again
mod.bc <- glm(shannon_bc ~ yr + as.factor(month), family = gaussian(), data = pdat)

# Check diagnostics
summary(mod.bc)


predicted <- predict(mod, type = "response")

# Plot observed vs. predicted
plot((sqrt(pdat$shannon)), predicted, main = "Observed vs. Predicted", 
     xlab = "Observed", ylab = "Predicted")
abline(a = 0, b = 1, col = "red")

residuals_gaus <- residuals(mod)
fitted_gaus <- fitted(mod)

qqnorm(residuals_gaus)
qqline(residuals_gaus, col = "red")

# Plot residuals vs. fitted values
plot(fitted_gaus, residuals_gaus, 
     xlab = "Fitted values", 
     ylab = "Residuals", 
     main = "Residuals vs Fitted Values")
abline(h = 0, col = "red")  # Add a horizontal line at 0 for reference

resids <- simulateResiduals(mod.bc)
plot(resids)

mod.bc2 <- glm(shannon_bc ~ poly(yr, 2) + as.factor(month), family = gaussian(), data = pdat)
mod.whatever <- glm((1/shannon) ~ poly(yr, 2) + as.factor(month), family = gaussian(), data = pdat)

resids <- simulateResiduals(mod.whatever)
plot(resids)



summary(mod.bc2)


AIC(mod.gaus, mod.gaus.2, mod.bc, mod.bc2, mod.whatever)



shapiro.test(residuals_gaus)
testOutliers(residuals)
testZeroInflation(resids)

residuals <- residuals(mod)
threshold <- 3 * sd(residuals)


# Find observations that exceed the threshold
outlier_indices <- which(abs(residuals) > threshold)

# Extract the specific rows of data that are outliers
outlier_data <- pdat[outlier_indices, ]

## try other estimators

plot(y = combined_estimates.df$jack1, x = combined_estimates.df$year)

combined_estimates.df$yr <- as.numeric(scale(combined_estimates.df$year))


mod.null <- glm(round(jack1, 0) ~ 1, data = combined_estimates.df, family = poisson)
mod.1 <- glm(round(jack1, 0) ~ yr + park, data = combined_estimates.df, family = poisson) 
mod_qp <- glm(round(jack1, 0) ~ year + park, family = quasipoisson, data = combined_estimates.df)

summary(mod.null)
summary(mod.1)
summary(mod_qp)


dispersion <- sum(residuals(mod.1, type = "pearson")^2) / mod.1$df.residual

#slight UNDERdispersion


## NOTE: You can run the models on each of the different estimates; always a small but significant
## negative trend


# Create a new data frame for prediction
predict_data <- expand.grid(
  yr = seq(min(combined_estimates.df$yr), max(combined_estimates.df$yr), length.out = 100),
  park = levels(as.factor(combined_estimates.df$park))
)

# Generate predictions on the link scale (logit)
predicted <- predict(mod.1, newdata = predict_data, type = "link", se.fit = TRUE)

# Convert predictions to a data frame and transform back to response scale
predicted_df <- cbind(predict_data, mu = predicted$fit, se = predicted$se.fit) %>%
  mutate(
    lower = exp(mu - 1.96 * se),
    upper = exp(mu + 1.96 * se),
    mu = exp(mu)
  )


ggplot() +
  geom_col(data = combined_estimates.df, aes(x = yr, y = Species)) +
  geom_line(data = predicted_df, aes(x = yr, y = mu)) +
  geom_ribbon(data = predicted_df, aes(x = yr, ymin = lower, ymax = upper), 
              alpha = 0.2, fill = bcs_colors["dark green"], inherit.aes = FALSE) +
  facet_wrap(~ park)

mean_predicted <- predicted_df %>%
  group_by(yr) %>%
  summarize(mean_mu = mean(mu),
            mean_lower = mean(lower),
            mean_upper = mean(upper))

mean_predicted$year <- (mean_predicted$yr * year_sd) + year_mean

year_mean <- mean(combined_estimates.df$year)
year_sd <- sd(combined_estimates.df$year)


combined_mean <- combined_estimates.df %>%
  group_by(year) %>%
  summarise(mean_jack = mean(jack1), .groups = "drop")

p.S_trend <- ggplot(combined_mean, aes(x = year, y = mean_jack)) +
  geom_point(alpha = .8, size = 1, color = bcs_colors["dark green"]) + 
  geom_line(data = mean_predicted, aes(x = year, y = mean_mu), color = bcs_colors["dark green"]) +
  geom_ribbon(data = mean_predicted, aes(x = year, ymin = mean_lower, ymax = mean_upper), 
              alpha = 0.2, fill = bcs_colors["dark green"], inherit.aes = FALSE) + 
  labs(title = "Average species richness at NBP sites over time", y = "Average number of species at NBP sites (first-order jackknife estimate)", 
       x = "Year") +
  theme_bcs()

png("results/figures/02_average_species_richness_over_time.png", width = 4, height = 4, units = "in", res = 300)
p.S_trend
dev.off()

beta_per_year <- summary(mod.1)$coefficients["yr", "Estimate"] / year_sd
se_per_year <- summary(mod.1)$coefficients["yr", "Std. Error"] / year_sd


mean_year_effect <- exp(beta_per_year)
upper_year_effect <- exp(beta_per_year + 1.96 * se_per_year)
lower_year_effect <- exp(beta_per_year - 1.96 * se_per_year)

total_change <- mean_year_effect^18
upper_change <- upper_year_effect^18
lower_change <- lower_year_effect^18

# Calculate the percentage change
percentage_change <- (total_change - 1) * 100
upper_percentage_change <- (upper_change - 1) * 100
lower_percentage_change <- (lower_change - 1) * 100

# Print the result
cat("Over the period of 2005 to 2023, we estimate average species richness at study locations to have changed by", 
    round(percentage_change, 2), "%, with a 95% confidence interval of [", 
    round(lower_percentage_change, 2), "%, ", round(upper_percentage_change, 2), "%].\n")



## test early vs late estimates of diversity ###

focal.parks <- c("Carkeek Park", "Discovery Park", "Genesee Park", "Golden Gardens Park", 
                 "Magnuson Park", "Seward Park", "Washington Park Arboretum")

srm.early <- dat %>% 
  filter(!str_detect(species, pattern = " x | sp\\.|Spotted Owl|Domestic"),
         year %in% c(2005:2015),
         park %in% focal.parks,
         station.code %in% covs$station.code) %>%
  mutate(observed = seen + heard + fly) %>% 
  group_by(park, year, station.code, bird.code) %>%
  summarise(count = sum(observed), .groups = "drop") %>%
  pivot_wider(names_from = bird.code, values_from = count, values_fill = 0)


srm.late <- dat %>% 
  filter(!str_detect(species, pattern = " x | sp\\.|Spotted Owl|Domestic"),
         year %in% c(2016:2019, 2021:2023),
         park %in% focal.parks,
         station.code %in% covs$station.code) %>%
  mutate(observed = seen + heard + fly) %>% 
  group_by(park, year, station.code, bird.code) %>%
  summarise(count = sum(observed), .groups = "drop") %>%
  pivot_wider(names_from = bird.code, values_from = count, values_fill = 0)

year.list.early <- list()
for(i in 1:length(focal.parks)){
  year.list.early[[i]] <- sort(unique(srm.early[srm.early$park == focal.parks[i], ]$year))
}

names(year.list.early) <- focal.parks

year.list.late <- list()
for(i in 1:length(focal.parks)){
  year.list.late[[i]] <- sort(unique(srm.late[srm.late$park == focal.parks[i], ]$year))
}

names(year.list.late) <- focal.parks





dat.list.early <- list()

for(i in 1:length(focal.parks)){
  
  for(y in year.list.early[[focal.parks[i]]]) {
    
    subset_data <- srm.early %>% filter(park == focal.parks[i] & year == y) %>% select(-park, -year, -station.code)
    
    dat.list.early[[paste0(focal.parks[i], "-", y)]] <- subset_data
  }
}


dat.list.late <- list()

for(i in 1:length(focal.parks)){
  
  for(y in year.list.late[[focal.parks[i]]]) {
    
    subset_data <- srm.late %>% filter(park == focal.parks[i] & year == y) %>% select(-park, -year, -station.code)
    
    dat.list.late[[paste0(focal.parks[i], "-", y)]] <- subset_data
  }
}


richness.estimates.early <- list()

for(i in 1:length(dat.list.early)){
  est <- specpool(dat.list.early[[i]])
  richness.estimates.early[[names(dat.list[i])]] <- est
}

combined_estimates <- bind_rows(richness.estimates)
head(combined_estimates)

combined_estimates$park.year <- names(dat.list)

combined_estimates.df <- combined_estimates %>%
  separate(park.year, into = c("park", "year"), sep = "-") %>%
  mutate(bootLower = boot - 1.96 * boot.se,
         bootUpper = boot + 1.96 * boot.se,
         year = as.numeric(year))

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

ggplot() +
  geom_col(data = plot.df %>% filter(estimator == "Observed richness"), 
           aes(x = year, y = richness_estimate, group = interaction(park, estimator), 
               fill = estimator)) +
  geom_line(data = plot.df %>% filter(estimator != "Observed richness"), 
            aes(x = year, y = richness_estimate, group = interaction(park, estimator), 
                color = estimator), size = 1, alpha = 0.8) +
  scale_color_manual(name = "Measure",
                     values = c("Chao1 estimate" = "#36BA3A", 
                                "First-order jackknife estimate" = "#FFB98C", 
                                "Second-order jackknife estimate" = "#0A3C23",
                                "Bootstrapped estimate"= "#E6FF55")) +
  scale_fill_manual(name = element_blank(),
                    values = c("Observed richness" = "#0A3C23")) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  facet_wrap(~ park) +
  labs(title = "Estimates of species richness over time", y = "Number of species", x = "Year") +
  theme_bcs() +
  theme(plot.title = element_text(hjust = 0.5))




                     
########## END ###########
