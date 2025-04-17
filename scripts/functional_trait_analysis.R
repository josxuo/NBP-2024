functional <- read_xlsx("data/processed/nbp_species_traits.xlsx", sheet = "traits")
nbp <- read_xlsx("data/processed/nbp_tidy_jan_24.xlsx")
covs <- read_csv("data/processed/circ_no_overlap_covariates.csv")

nbp[is.na(nbp$bird.code), ]$species

nbp[nbp$species == "Glaucous-winged x Western Gull", ]$bird.code <- "OLGU"
nbp[nbp$species == "American x Eurasian Wigeon", ]$bird.code <- "AEWI"

nbp[nbp$species == "Glaucous-winged x Western Gull", ]$bird.code <- "GWGU"
nbp[nbp$species == "American x Eurasian Wigeon", ]$bird.code <- "AMWI"

d1 <- nbp %>%
  filter(!str_detect(species, pattern = " sp\\.|Spotted Owl|Domestic"),
         station.code %in% covs$station.code,
         year %in% c(2005:2019, 2021, 2022, 2023)) %>%
        # exclusions == "none") %>%
  mutate(loop2 = as.factor(paste(park, loop, sep = "-")),
         park = as.factor(park),
         month = as.factor(month),
         station.code = as.factor(station.code),
         station.year = paste(station.code, year, sep = "-"),
         park.year = paste(park, year, sep = "-"))



traitsdattraits <- functional %>% select(-common.name, -scientific.name)

dat1 <- d1 %>% left_join(., traits, join_by("bird.code" == "bird.code"))


func1 <- dat1 %>% group_by(year, habitat) %>% reframe(detected = n_distinct(survey_id))



nsurv1 <- dat1 %>% group_by(year) %>% reframe(nsurv = n_distinct(survey_id))
 


func1 <- left_join(func1, nsurv1, join_by("year" == "year")) %>%
  mutate(not_detected = nsurv - detected, prop = detected / nsurv)



ggplot(func1, aes(x = year, y = prop)) +
  geom_line(color = bcs_colors["dark green"]) +
  facet_wrap(~ habitat, scales = "free") +
  labs(x= "Year", y = "Proportion of surveys reporting") +
  theme_bcs()

# Aerial Insectivores

AI <- dat1 %>% group_by(year, AI) %>% reframe(detected = n_distinct(survey_id)) %>%
  left_join(., nsurv1) %>%
  mutate(not_detected = nsurv - detected, prop = detected / nsurv)

ggplot(AI %>% filter(AI == "AI"), aes(x = year, y = prop, color = AI)) +
  geom_line(color = bcs_colors["dark green"]) +
  theme_bcs()

func1$yr <- as.numeric(scale(func1$year))


mod.null <- glmmTMB(cbind(detected, not_detected) ~ 1, data = func1, family = binomial)
mod1 <- glmmTMB(cbind(detected, not_detected) ~ habitat, data = func1, family = binomial)
mod2 <- glmmTMB(cbind(detected, not_detected) ~ yr * habitat, data = func1, family = binomial)


AIC(mod.null, mod1, mod2)

summary(mod2)


functional[functional$AI == "AI", ]$common.name

# Bird Group

group <- dat1 %>% group_by(year, bird.group) %>% reframe(detected = n_distinct(survey_id)) %>%
  left_join(., nsurv1) %>%
  mutate(not_detected = nsurv - detected, prop = detected / nsurv)

ggplot(group, aes(x = year, y = prop)) +
  geom_line(color = bcs_colors["dark green"], linewidth = 1) +
  facet_wrap(~ bird.group, scales = "free") +
  labs(x = "Year", y = "Proportion of surveys reporting") +
  theme_bcs()

functional[functional$bird.group == "other", ]$common.name

# Migratory behavior



Mig <- dat1 %>% group_by(year, migrate) %>% reframe(detected = n_distinct(survey_id)) %>%
  left_join(., nsurv1) %>%
  mutate(not_detected = nsurv - detected, prop = detected / nsurv)

ggplot(mig, aes(x = year, y = prop, color = migrate)) +
  geom_line() +
  #facet_wrap(~ migrate, scales = "free") +
  theme_bcs()

# trophic niche

niche <- dat1 %>% group_by(year, trophic.niche) %>% reframe(detected = n_distinct(survey_id)) %>%
  left_join(., nsurv1) %>%
  mutate(not_detected = nsurv - detected, prop = detected / nsurv)

ggplot(niche %>% filter(trophic.niche != "Scavenger"), aes(x = year, y = prop)) +
  geom_line(color = bcs_colors["dark green"], linewidth = 1) +
  facet_wrap(~ trophic.niche, scales = "free") +
  labs(y = "Proportion of surveys reporting", x = "Year") +
  theme_bcs()




focal.trait <- "breeding.biome"

p.trait <- dat1 %>% group_by(year, winter.biome) %>% reframe(detected = n_distinct(survey_id)) %>%
  left_join(., nsurv1) %>%
  mutate(not_detected = nsurv - detected, prop = detected / nsurv)

ggplot(p.trait, aes(x = year, y = prop)) +
  geom_line(color = bcs_colors["dark green"], linewidth = 1) +
  facet_wrap(~ winter.biome, scales = "free") +
  labs(y = "Proportion of surveys reporting", x = "Year") +
  theme_bcs()



df <- dat1 %>% group_by(year, bird.code) %>%
  reframe(detected = n_distinct(survey_id)) %>% 
  pivot_wider(names_from = bird.code, values_from = detected, values_fill = list(detected = 0)) %>%
  pivot_longer(-1, names_to = "bird.code", values_to = "detected") %>%
  left_join(., nsurv1) %>%
  mutate(not_detected = nsurv - detected, prop = detected / nsurv,
         yr = as.numeric(scale(year)))


df <- left_join(df, traits, join_by("bird.code" == "bird.code")) %>%
  mutate(mass = as.numeric(scale(body.mass)))


view(df)

g <- glm(cbind(detected, not_detected) ~ trophic.niche + habitat + bird.group + AI + yr + mass +
           migrate + native + family + breeding.biome + winter.biome,
                    data = df,
                    family = binomial,
         na.action = "na.fail")

g.zi <- glmmTMB(cbind(detected, not_detected) ~ trophic.niche + habitat + bird.group + AI + yr + mass +
           migrate + breeding.biome,
         data = df,
         zi = ~ 1,
         family = binomial,
         na.action = "na.fail")

g.beta.zi <- glmmTMB(cbind(detected, not_detected) ~ habitat + bird.group + AI + yr + mass,
                data = df,
                zi = ~ 1,
                family = qubinomial,
                na.action = "na.fail")

plot(simulateResiduals(g.zi))
testZeroInflation(g.zi)


summary(g.zi)

mod <- glm(cbind(detected, not_detected) ~ trophic.niche + habitat  + AI + mass + yr + 
             migrate + family + breeding.biome + winter.biome,
           data = df.clean,
           family = binomial,
           na.action = "na.fail")

mod.quasi <- glm(cbind(detected, not_detected) ~ trophic.niche + habitat  + AI + mass + yr + 
             migrate + family + breeding.biome + winter.biome,
           data = df.clean,
           family = quasibinomial,
           na.action = "na.fail")

mod.zi <- glmmTMB(cbind(detected, not_detected) ~  (yr | bird.code),
                  zi = ~1,   
                  data = df.clean,
                    family = binomial)

testZeroInflation(g)

summary(mod.quasi)
plot(resid(mod.quasi), main = "quasi")
plot(resid(mod), main = "binomial")
dispersion <- sum(resid(mod.quasi, type = "pearson")^2) / df.residual(mod.quasi)
dispersion

plot(simulateResiduals(mod.zi))
AIC(mod, mod.quasi)

dredge(g)

rm(global_model)
?dredge

summary(global_model)


quantile(traits$body.mass, probs = c(.33, .67))

traits <- traits %>% mutate(body_size = ifelse(body.mass <= 37.7, "small",
                                               ifelse(body.mass >= 412.4, "large", 
                                                      "medium")))

p.trait <- dat1 %>% group_by(year, AI) %>% reframe(detected = n_distinct(survey_id)) %>%
  left_join(., nsurv1) %>%
  mutate(not_detected = nsurv - detected, prop = detected / nsurv)

ggplot(p.trait %>% filter(AI == "AI"), aes(x = year, y = prop)) +
  geom_line(linewidth = 1) +
  #facet_wrap(~ body_size, scales = "free") +
  labs(y = "Proportion of surveys reporting", x = "Year",
       color = "Aerial Insectivore") +
  theme_bcs()
sort(unique(nbp$park))



park.survs <- d1 %>% group_by(year, park) %>%
  summarise(nsurv = n_distinct(survey_id), .groups = "drop") %>%
  mutate(park.year = paste(park, year, sep = "-"))

park.props <- dat1 %>% group_by(park.year, bird.code) %>%
  reframe(detected = n_distinct(survey_id)) %>% 
  pivot_wider(names_from = bird.code, values_from = detected, values_fill = list(detected = 0)) %>%
  pivot_longer(-1, names_to = "bird.code", values_to = "detected") %>%
  left_join(., park.survs) %>%
  mutate(not_detected = nsurv - detected, prop = detected / nsurv)

# CARKEEK

carkeek <- park.props %>% filter(park == "Carkeek Park") %>% mutate(yr = as.numeric(scale(year)))

spp <- park.props %>% filter(park == "Carkeek Park") %>% group_by(bird.code) %>% 
  summarise(nyrs_det = sum(prop > 0)) %>% filter(nyrs_det >= 5) %>% pull(bird.code)

disp <- vector()
r2 <- vector()
beta_yr <- vector()
se_yr <- vector()
p_yr <- vector()

for(i in 1:length(spp)){
  mod.dat <- carkeek[carkeek$bird.code == spp[i], ]
  mod.null <- glmmTMB(cbind(detected, not_detected) ~ 1, data = mod.dat, family = binomial)
  mod <- glmmTMB(cbind(detected, not_detected) ~ yr, data = mod.dat, family = binomial)
  r2[i] <- 1 - (logLik(mod) / logLik(mod.null))
  disp[i] <- sum(residuals(mod, type = "pearson")^2) / nrow(mod.dat) - length(fixef(mod)$cond)
  beta_yr[i] <- summary(mod)$coefficients$cond["yr", 1]
  se_yr[i] <- summary(mod)$coefficients$cond["yr", 2]
  p_yr[i] <- summary(mod)$coefficients$cond["yr", 4]
}

#for(i in 1:length(spp)){
 # mod.dat <- carkeek[carkeek$bird.code == spp[i], ]
 # mod.null <- glm(cbind(detected, not_detected) ~ 1, data = mod.dat, family = binomial)
 # mod <- glm(cbind(detected, not_detected) ~ yr, data = mod.dat, family = binomial)
  #r2[i] <- 1 - (logLik(mod) / logLik(mod.null))
 # disp[i] <- sum(residuals(mod, type = "pearson")^2) / mod$df.residual
 # beta_yr[i] <- summary(mod)$coefficients["yr", 1]
 # se_yr[i] <- summary(mod)$coefficients["yr", 2]
 # p_yr[i] <- summary(mod)$coefficients["yr", 4]
}

carkeek_results <- as.data.frame(cbind(spp, disp, r2, beta_yr, se_yr, p_yr))


# DISCOVERY

discovery <- park.props %>% filter(park == "Discovery Park") %>% mutate(yr = as.numeric(scale(year)))

spp <- park.props %>% filter(park == "Discovery Park") %>% group_by(bird.code) %>% 
  summarise(nyrs_det = sum(prop > 0)) %>% filter(nyrs_det >= 5) %>% pull(bird.code)

disp <- vector()
r2 <- vector()
beta_yr <- vector()
se_yr <- vector()
p_yr <- vector()

for(i in 1:length(spp)){
  mod.dat <- discovery[discovery$bird.code == spp[i], ]
  mod.null <- glmmTMB(cbind(detected, not_detected) ~ 1, data = mod.dat, family = binomial)
  mod <- glmmTMB(cbind(detected, not_detected) ~ yr, data = mod.dat, family = binomial)
  r2[i] <- 1 - (logLik(mod) / logLik(mod.null))
  disp[i] <- sum(residuals(mod, type = "pearson")^2) / nrow(mod.dat) - length(fixef(mod)$cond)
  beta_yr[i] <- summary(mod)$coefficients$cond["yr", 1]
  se_yr[i] <- summary(mod)$coefficients$cond["yr", 2]
  p_yr[i] <- summary(mod)$coefficients$cond["yr", 4]
}

discovery_results <- as.data.frame(cbind(spp, disp, r2, beta_yr, se_yr, p_yr))
#view(discovery_results %>% filter(p_yr < 0.05))

# GENESEE

genesee <- park.props %>% filter(park == "Genesee Park") %>% mutate(yr = as.numeric(scale(year)))

spp <- park.props %>% filter(park == "Genesee Park") %>% group_by(bird.code) %>% 
  summarise(nyrs_det = sum(prop > 0)) %>% filter(nyrs_det >= 5) %>% pull(bird.code)

disp <- vector()
r2 <- vector()
beta_yr <- vector()
se_yr <- vector()
p_yr <- vector()

for(i in 1:length(spp)){
  mod.dat <- genesee[genesee$bird.code == spp[i], ]
  mod.null <- glmmTMB(cbind(detected, not_detected) ~ 1, data = mod.dat, family = binomial)
  mod <- glmmTMB(cbind(detected, not_detected) ~ yr, data = mod.dat, family = binomial)
  r2[i] <- 1 - (logLik(mod) / logLik(mod.null))
  disp[i] <- sum(residuals(mod, type = "pearson")^2) / nrow(mod.dat) - length(fixef(mod)$cond)
  beta_yr[i] <- summary(mod)$coefficients$cond["yr", 1]
  se_yr[i] <- summary(mod)$coefficients$cond["yr", 2]
  p_yr[i] <- summary(mod)$coefficients$cond["yr", 4]
}

genesee_results <- as.data.frame(cbind(spp, disp, r2, beta_yr, se_yr, p_yr))
#view(genesee_results %>% filter(p_yr < 0.05))

# GOLDEN GARDENS

golden <- park.props %>% filter(park == "Golden Gardens Park") %>% mutate(yr = as.numeric(scale(year)))

spp <- park.props %>% filter(park == "Golden Gardens Park") %>% group_by(bird.code) %>% 
  summarise(nyrs_det = sum(prop > 0)) %>% filter(nyrs_det >= 5) %>% pull(bird.code)

disp <- vector()
r2 <- vector()
beta_yr <- vector()
se_yr <- vector()
p_yr <- vector()

for(i in 1:length(spp)){
  mod.dat <- golden[golden$bird.code == spp[i], ]
  mod.null <- glmmTMB(cbind(detected, not_detected) ~ 1, data = mod.dat, family = binomial)
  mod <- glmmTMB(cbind(detected, not_detected) ~ yr, data = mod.dat, family = binomial)
  r2[i] <- 1 - (logLik(mod) / logLik(mod.null))
  disp[i] <- sum(residuals(mod, type = "pearson")^2) / nrow(mod.dat) - length(fixef(mod)$cond)
  beta_yr[i] <- summary(mod)$coefficients$cond["yr", 1]
  se_yr[i] <- summary(mod)$coefficients$cond["yr", 2]
  p_yr[i] <- summary(mod)$coefficients$cond["yr", 4]
}

golden_results <- as.data.frame(cbind(spp, disp, r2, beta_yr, se_yr, p_yr))
#view(golden_results %>% filter(p_yr < 0.05))


# LAKE FOREST PARK

lfp <- park.props %>% filter(park == "Lake Forest Park") %>% mutate(yr = as.numeric(scale(year)))

spp <- park.props %>% filter(park == "Lake Forest Park") %>% group_by(bird.code) %>% 
  summarise(nyrs_det = sum(prop > 0)) %>% filter(nyrs_det >= 5) %>% pull(bird.code)

disp <- vector()
r2 <- vector()
beta_yr <- vector()
se_yr <- vector()
p_yr <- vector()

for(i in 1:length(spp)){
  mod.dat <- lfp[lfp$bird.code == spp[i], ]
  mod.null <- glmmTMB(cbind(detected, not_detected) ~ 1, data = mod.dat, family = binomial)
  mod <- glmmTMB(cbind(detected, not_detected) ~ yr, data = mod.dat, family = binomial)
  r2[i] <- 1 - (logLik(mod) / logLik(mod.null))
  disp[i] <- sum(residuals(mod, type = "pearson")^2) / nrow(mod.dat) - length(fixef(mod)$cond)
  beta_yr[i] <- summary(mod)$coefficients$cond["yr", 1]
  se_yr[i] <- summary(mod)$coefficients$cond["yr", 2]
  p_yr[i] <- summary(mod)$coefficients$cond["yr", 4]
}

lake_forest_results <- as.data.frame(cbind(spp, disp, r2, beta_yr, se_yr, p_yr))

# LINCOLN PARK

linc <- park.props %>% filter(park == "Lincoln Park") %>% mutate(yr = as.numeric(scale(year)))

spp <- park.props %>% filter(park == "Lincoln Park") %>% group_by(bird.code) %>% 
  summarise(nyrs_det = sum(prop > 0)) %>% filter(nyrs_det >= 5) %>% pull(bird.code)

disp <- vector()
r2 <- vector()
beta_yr <- vector()
se_yr <- vector()
p_yr <- vector()

for(i in 1:length(spp)){
  mod.dat <- linc[linc$bird.code == spp[i], ]
  mod.null <- glmmTMB(cbind(detected, not_detected) ~ 1, data = mod.dat, family = binomial)
  mod <- glmmTMB(cbind(detected, not_detected) ~ yr, data = mod.dat, family = binomial)
  r2[i] <- 1 - (logLik(mod) / logLik(mod.null))
  disp[i] <- sum(residuals(mod, type = "pearson")^2) / nrow(mod.dat) - length(fixef(mod)$cond)
  beta_yr[i] <- summary(mod)$coefficients$cond["yr", 1]
  se_yr[i] <- summary(mod)$coefficients$cond["yr", 2]
  p_yr[i] <- summary(mod)$coefficients$cond["yr", 4]
}

lincoln_results <- as.data.frame(cbind(spp, disp, r2, beta_yr, se_yr, p_yr))

# MAGNUSON PARK

mag <- park.props %>% filter(park == "Magnuson Park") %>% mutate(yr = as.numeric(scale(year)))

spp <- park.props %>% filter(park == "Magnuson Park") %>% group_by(bird.code) %>% 
  summarise(nyrs_det = sum(prop > 0)) %>% filter(nyrs_det >= 5) %>% pull(bird.code)

disp <- vector()
r2 <- vector()
beta_yr <- vector()
se_yr <- vector()
p_yr <- vector()

for(i in 1:length(spp)){
  mod.dat <- mag[mag$bird.code == spp[i], ]
  mod.null <- glmmTMB(cbind(detected, not_detected) ~ 1, data = mod.dat, family = binomial)
  mod <- glmmTMB(cbind(detected, not_detected) ~ yr, data = mod.dat, family = binomial)
  r2[i] <- 1 - (logLik(mod) / logLik(mod.null))
  disp[i] <- sum(residuals(mod, type = "pearson")^2) / nrow(mod.dat) - length(fixef(mod)$cond)
  beta_yr[i] <- summary(mod)$coefficients$cond["yr", 1]
  se_yr[i] <- summary(mod)$coefficients$cond["yr", 2]
  p_yr[i] <- summary(mod)$coefficients$cond["yr", 4]
}

magnuson_results <- as.data.frame(cbind(spp, disp, r2, beta_yr, se_yr, p_yr))

# SEWARD PARK

seward <- park.props %>% filter(park == "Seward Park") %>% mutate(yr = as.numeric(scale(year)))

spp <- park.props %>% filter(park == "Seward Park") %>% group_by(bird.code) %>% 
  summarise(nyrs_det = sum(prop > 0)) %>% filter(nyrs_det >= 5) %>% pull(bird.code)

disp <- vector()
r2 <- vector()
beta_yr <- vector()
se_yr <- vector()
p_yr <- vector()

for(i in 1:length(spp)){
  mod.dat <- seward[seward$bird.code == spp[i], ]
  mod.null <- glmmTMB(cbind(detected, not_detected) ~ 1, data = mod.dat, family = binomial)
  mod <- glmmTMB(cbind(detected, not_detected) ~ yr, data = mod.dat, family = binomial)
  r2[i] <- 1 - (logLik(mod) / logLik(mod.null))
  disp[i] <- sum(residuals(mod, type = "pearson")^2) / nrow(mod.dat) - length(fixef(mod)$cond)
  beta_yr[i] <- summary(mod)$coefficients$cond["yr", 1]
  se_yr[i] <- summary(mod)$coefficients$cond["yr", 2]
  p_yr[i] <- summary(mod)$coefficients$cond["yr", 4]
}

seward_results <- as.data.frame(cbind(spp, disp, r2, beta_yr, se_yr, p_yr))

# WASHINGTON PARK ARBORETUM

arb <- park.props %>% filter(park == "Washington Park Arboretum") %>% mutate(yr = as.numeric(scale(year)))

spp <- park.props %>% filter(park == "Washington Park Arboretum") %>% group_by(bird.code) %>% 
  summarise(nyrs_det = sum(prop > 0)) %>% filter(nyrs_det >= 5) %>% pull(bird.code)

disp <- vector()
r2 <- vector()
beta_yr <- vector()
se_yr <- vector()
p_yr <- vector()

for(i in 1:length(spp)){
  mod.dat <- arb[arb$bird.code == spp[i], ]
  mod.null <- glmmTMB(cbind(detected, not_detected) ~ 1, data = mod.dat, family = binomial)
  mod <- glmmTMB(cbind(detected, not_detected) ~ yr, data = mod.dat, family = binomial)
  r2[i] <- 1 - (logLik(mod) / logLik(mod.null))
  disp[i] <- sum(residuals(mod, type = "pearson")^2) / nrow(mod.dat) - length(fixef(mod)$cond)
  beta_yr[i] <- summary(mod)$coefficients$cond["yr", 1]
  se_yr[i] <- summary(mod)$coefficients$cond["yr", 2]
  p_yr[i] <- summary(mod)$coefficients$cond["yr", 4]
}

arb_results <- as.data.frame(cbind(spp, disp, r2, beta_yr, se_yr, p_yr))

carkeek_results$park <- rep("Carkeek Park", dim(carkeek_results)[1])
discovery_results$park <- rep("Discovery Park", dim(discovery_results)[1])
genesee_results$park <- rep("Genesee Park", dim(genesee_results)[1])
golden_results$park <- rep("Golden Gardens Park", dim(golden_results)[1])
lake_forest_results$park <- rep("Lake Forest Park", dim(lake_forest_results)[1])
lincoln_results$park <- rep("Lincoln Park", dim(lincoln_results)[1])

magnuson_results$park <- rep("Magnuson Park", dim(magnuson_results)[1])
seward_results$park <- rep("Seward Park", dim(seward_results)[1])
arb_results$park <- rep("Washington Park Arboretum", dim(arb_results)[1])

trends_by_park <- bind_rows(carkeek_results, discovery_results, genesee_results, golden_results,
                            lake_forest_results, lincoln_results, magnuson_results, seward_results,
                            arb_results) %>%
  left_join(., traits, join_by("spp" == "bird.code"))

write.csv(trends_by_park, "results/tables/trends_by_park.csv", row.names = FALSE)

park_lists <- park.props %>% group_by(park, bird.code) %>% 
  summarise(nyrs_det = sum(prop > 0)) %>% filter(nyrs_det > 0)

totals <- park_lists %>% group_by(park) %>% summarise(total_assessed = sum(nyrs_det >= 5),
                                                      total_detected = n_distinct(bird.code))


significant <- trends_by_park %>% filter(p_yr < 0.05) %>% group_by(park) %>% 
  summarise(decline = sum(beta_yr < 0), increase = sum(beta_yr > 0))

marginal <- trends_by_park %>% filter(p_yr >= 0.05 & p_yr < 0.1) %>% group_by(park) %>% 
  summarise(decline = sum(beta_yr < 0), increase = sum(beta_yr > 0))

summary(trends_by_park)

ggplot(trends_by_park %>% filter(p_yr < 0.05), 
       aes(x = as.numeric(p_yr), y = as.numeric(beta_yr), color = bird.group)) +
  geom_point() + 
  facet_wrap(~ park, scales = "free")

trend_tbl <- left_join(significant, totals)
