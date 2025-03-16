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
  mutate(species = fct_reorder(species, n))

library(ggbreak)

punr <- ggplot(dunr, aes(x = species, y = n)) + 
  geom_col(fill = bcs_colors["dark green"]) +
  coord_flip() +
  labs(title = "NBP observations not identified to species level", y = "Number of times reported", x = "") +
  #scale_y_break(c(270, 1350), expand = expansion(mult = c(0, 0))) +  # Adjust expansion to remove space
  scale_y_continuous(limits = c(0, 1400), breaks = seq(0, 1400, by = 50), 
                     expand = expansion(mult = c(0, 0)), position = "left") + 
  theme_bcs() +
  theme(
    panel.grid.major.y = element_blank()  # Remove unnecessary grid lines
  )

punr


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

## Details survey history and effort at each location
history <- dat %>% group_by(park, year) %>% 
  summarise(first = min(survey_date), last = max(survey_date), ncirc = n_distinct(station.code), 
          nsurv = n_distinct(survey_id), .groups = "drop") %>%
  mutate(park.year = paste(park, year, sep = "-")) %>%
  select(park.year, first, last, ncirc, nsurv)

write.csv(history, "results/survey_history_by_park.csv", row.names = FALSE)

## NOTE Magnuson has the weird station IDs

unique(dat[dat$park == "Magnuson Park", ]$station.code)

## is species richness correlated with longer survey history?
Sprk <- dat %>% 
  filter(!str_detect(species, pattern = " x | sp\\.|Spotted Owl")) %>%
  group_by(park, year) %>% 
  summarise(nsp = n_distinct(species), .groups = "drop") %>% 
  mutate(park.year = paste(park, year, sep = "-")) %>%
  left_join(., history, join_by("park.year" == "park.year")) %>% 
  mutate(tenure = as.numeric(last - first))


view(history)
Sprk$sncirc <- as.numeric(scale(Sprk$ncirc))
Sprk$snsurv <- as.numeric(scale(Sprk$nsurv))
Sprk$yr <- Sprk$year - min(Sprk$year)

mod.null <- glm(nsp ~ 1, data = Sprk, family = poisson)
mod.1 <- glm(nsp ~ park, data = Sprk, family = poisson)
mod.2 <- glm(nsp ~ park + yr, data = Sprk, family = poisson)
mod.3 <- glm(nsp ~ sncirc, data = Sprk, family = poisson)
mod.4 <- glm(nsp ~ snsurv, data = Sprk, family = poisson)
mod.5 <- glm(nsp ~ sncirc + snsurv, data = Sprk, family = poisson)
mod.6 <- glm(nsp ~ sncirc + snsurv + park, data = Sprk, family = poisson)
mod.7 <- glm(nsp ~ sncirc + snsurv + park + yr, data = Sprk, family = poisson)

mod_list <- list("null" = mod.null, "1" = mod.1, "2" = mod.2, "3" = mod.3, "4" = mod.4, 
                 "5" = mod.5, "6" = mod.6, "7" = mod.7)

aictab(mod_list)

summary(mod.7)


dispersion_test <- simulateResiduals(mod.7)
plot(dispersion_test)


library(glmmTMB)

mod.7.nb1 <- glmmTMB(nsp ~ sncirc + snsurv + park + yr, data = Sprk, family = nbinom1)
mod.7.nb2 <- glmmTMB(nsp ~ sncirc + snsurv + park + yr, data = Sprk, family = nbinom2)

mod_list <- list("nb1" = mod.7.nb1, "nb2" = mod.7.nb2)

aictab(mod_list)
aictab(list(mod.7))

dispersion_test <- simulateResiduals(mod.7.nb1)
plot(dispersion_test)

mod.7.nb1.ran <- glmmTMB(nsp ~ sncirc + snsurv + yr + (1 | park), data = Sprk, family = nbinom1)
mod.7.nb1.ran2 <- glmmTMB(nsp ~ sncirc + snsurv + (yr | park), data = Sprk, family = nbinom1)

mod_list <- list("ran1" = mod.7.nb1.ran, "ran2" = mod.7.nb1.ran2)
aictab(mod_list)

summary(mod.7.nb1)

mod.7.nb1.short <- glmmTMB(nsp ~ snsurv + park + yr, data = Sprk, family = nbinom1)
mod.7.nb1.short.ry <- glmmTMB(nsp ~ snsurv + park + (1 | yr), data = Sprk, family = nbinom1)
mod.7.short.ry <- glmmTMB(nsp ~ snsurv + park + (1 | yr), data = Sprk, family = poisson)

mod_gamma <- glm(nsp ~ snsurv + park + yr, data = Sprk, family = Gamma(link = "log"))

library(tweedie)
library(statmod)

mod_tweedie <- glm(nsp ~ snsurv + park + yr, 
                   family = tweedie(var.power = 1.5), 
                   data = Sprk)

summary(mod_tweedie)

plot(mod_tweedie$fitted.values, residuals(mod_tweedie),
     main = "Residuals vs Fitted", 
     xlab = "Fitted Values", 
     ylab = "Residuals")
abline(h = 0, col = "red")

qqnorm(residuals(mod_tweedie))
qqline(residuals(mod_tweedie), col = "red")

deviance(mod_tweedie) / df.residual(mod_tweedie)

library(mgcv)

AIC(mod.null, mod.7, mod.7.nb1, mod.7.nb1.ran, mod.7.nb1.short, mod.7.nb1.short.ry, mod.7.short.ry, mod_gamma, mod_tweedie)


dispersion_test <- simulateResiduals(mod_tweedie)
plot(dispersion_test)

ggplot(Sprk, aes(x = yr, y = nsp, color = park)) +
  geom_line()

## most support for mod.5. Effect size of scaled tenure = 29.289 per unit increase

# to convert effect size to years, divide by stdv of tenure and multiple by 365 days
29.289 / sd(Sprk$tenure) * 365

ggplot(nsp.predict, aes(x = tenure / 365, y = nsp)) +
  geom_point(color = bcs_colors["dark green"], size = 2) +
  geom_smooth(method = "lm", fill = bcs_colors["peach"], color = bcs_colors["dark green"]) + 
  labs(x = "Years surveyed", y = "Cumulative species richness", title = "Cumulative speices reported", subtitle = "vs. number of years surveyed") +
  theme_bcs()



library(vegan)

focal.parks <- c("Carkeek Park", "Discovery Park", "Genesee Park", "Golden Gardens Park", "Lake Forest Park", 
                 "Magnuson Park", "Seward Park", "Washington Park Arboretum")

srm <- dat %>% 
  filter(!str_detect(species, pattern = " x | sp\\.|Spotted Owl"),
         year %in% c(2005:2019, 2021, 2022, 2023)) %>%
  mutate(observed = seen + heard + fly) %>% 
  group_by(park, year, station.code, bird.code) %>%
  summarise(count = sum(observed), .groups = "drop") %>%
  pivot_wider(names_from = bird.code, values_from = count, values_fill = 0)

srm.focal <- srm %>% filter(park %in% focal.parks)


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


##

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

combined_estimates.df$residuals <- residuals(mod.1, type = "pearson")  # Pearson residuals
combined_estimates.df$fitted <- fitted(mod.1)  # Fitted values

ggplot(combined_estimates.df, aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(x = "Fitted Values", y = "Pearson Residuals", title = "Residual vs. Fitted Plot")


ggplot(combined_estimates.df, aes(x = residuals)) +
  geom_histogram(binwidth = 0.5, fill = "blue", alpha = 0.5, color = "black") +
  theme_minimal() +
  labs(x = "Residuals", y = "Frequency", title = "Histogram of Residuals")

mod.2 <- glmmTMB(round(jack1, 0) ~ yr + park + (yr | park),
                 data = combined_estimates.df, family = poisson)


summary(mod.1)

combined_estimates.df$park <- as.factor(combined_estimates.df$park)

# Create a new data frame for prediction
predict_data <- expand.grid(
  yr = seq(min(combined_estimates.df$yr), max(combined_estimates.df$yr), length.out = 100),
  park = levels(combined_estimates.df$park)
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

ggplot(combined_mean, aes(x = year, y = mean_jack)) +
  geom_point(alpha = .8, size = 2, color = bcs_colors["dark green"]) + 
  geom_line(data = mean_predicted, aes(x = year, y = mean_mu), color = bcs_colors["dark green"]) +
  geom_ribbon(data = mean_predicted, aes(x = year, ymin = mean_lower, ymax = mean_upper), 
              alpha = 0.2, fill = bcs_colors["dark green"], inherit.aes = FALSE) + 
  labs(title = "Trend in species richness over time", y = "Estimated mean number of species at NBP sites", 
       x = "Year") +
  theme_bcs()

beta_per_year <- summary(mod.1)$coefficients["yr", "Estimate"] / year_sd
se_per_year <- summary(mod.1)$coefficients["yr", "Std. Error"] / year_sd


mean_year_effect <- exp(beta_per_year)
upper_year_effect <- exp(beta_per_year + 1.96 * se_per_year)
lower_year_effect <- exp(beta_per_year - 1.96 * se_per_year)


## NOTE: You can run the models on each of the different estimates; always a small but significant
## negative trend

mod.null <- glm(round(boot, 0) ~ 1, data = combined_estimates.df, family = poisson)
mod.1 <- glm(round(boot, 0) ~ yr + park, data = combined_estimates.df, family = poisson) 
mod_qp <- glm(round(boot, 0) ~ yr + park, family = quasipoisson, data = combined_estimates.df)

summary(mod.null)
summary(mod.1)
summary(mod_qp)


dispersion <- sum(residuals(mod.1, type = "pearson")^2) / mod.1$df.residual

#slight UNDERdispersion

combined_estimates.df$residuals <- residuals(mod.1, type = "pearson")  # Pearson residuals
combined_estimates.df$fitted <- fitted(mod.1)  # Fitted values

ggplot(combined_estimates.df, aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(x = "Fitted Values", y = "Pearson Residuals", title = "Residual vs. Fitted Plot")


ggplot(combined_estimates.df, aes(x = residuals)) +
  geom_histogram(binwidth = 0.5, fill = "blue", alpha = 0.5, color = "black") +
  theme_minimal() +
  labs(x = "Residuals", y = "Frequency", title = "Histogram of Residuals")

mod.2 <- glmmTMB(round(jack1, 0) ~ yr + park + (yr | park),
                 data = combined_estimates.df, family = poisson)


summary(mod.1)

combined_estimates.df$park <- as.factor(combined_estimates.df$park)

# Create a new data frame for prediction
predict_data <- expand.grid(
  yr = seq(min(combined_estimates.df$yr), max(combined_estimates.df$yr), length.out = 100),
  park = levels(combined_estimates.df$park)
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


                     
########## END ###########
