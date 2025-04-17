library(tidyverse)
library(glmmTMB)
library(simr)


trends <- read.csv("results/tables/trends_by_park.csv")

source("functions/theme_bcs.R")

trends %>% filter(park == "Carkeek Park" & p_yr < 0.05) %>% arrange(desc(r2))
nsurv1 <- dat1 %>% group_by(station.year) %>%
  summarise(nsurv = n_distinct(survey_id))


park.props <- dat1 %>% group_by(year, station.code, station.year, bird.code, park) %>%
  reframe(detected = n_distinct(survey_id)) %>% 
  pivot_wider(names_from = bird.code, values_from = detected, values_fill = list(detected = 0)) %>%
  pivot_longer(-c(1,2, 3, 4), names_to = "bird.code", values_to = "detected") %>%
  left_join(., nsurv1) %>%
  #filter(bird.code == "CORA" & park == "Carkeek Park") %>%
  mutate(not_detected = nsurv - detected, prop = detected / nsurv,
         yr = as.numeric(scale(year)))

#mod.null <- glmmTMB(cbind(detected, not_detected) ~ (1 | station.code), data = mod.dat, family = binomial)
mod <- glmmTMB(cbind(detected, not_detected) ~ yr + (1 | station.code), data = mod.dat, family = binomial)
nak <- r2_nakagawa(mod)


# CARKEEK

carkeek <- park.props %>% filter(park == "Carkeek Park") %>% mutate(yr = as.numeric(scale(year)))

spp <- park.props %>% filter(park == "Carkeek Park") %>% group_by(bird.code) %>% 
  summarise(nyrs_det = sum(prop > 0)) %>% filter(nyrs_det >= 5) %>% pull(bird.code)

disp <- numeric()
r2_conditional <- numeric()
r2_marginal <- numeric()
beta_yr <- numeric()
se_yr <- numeric()
p_yr <- numeric()

for(i in 1:length(spp)){
  mod.dat <- carkeek[carkeek$bird.code == spp[i], ]
  #mod.null <- glmmTMB(cbind(detected, not_detected) ~ 1, data = mod.dat, family = binomial)
  mod <- glmmTMB(cbind(detected, not_detected) ~ yr + (1 | station.code), data = mod.dat, family = binomial)
  r2_conditional[i] <- r2_nakagawa(mod)$R2_conditional
  r2_marginal[i] <- r2_nakagawa(mod)$R2_marginal
  disp[i] <- sum(residuals(mod, type = "pearson")^2) / nrow(mod.dat) - length(fixef(mod)$cond)
  beta_yr[i] <- summary(mod)$coefficients$cond["yr", 1]
  se_yr[i] <- summary(mod)$coefficients$cond["yr", 2]
  p_yr[i] <- summary(mod)$coefficients$cond["yr", 4]
}

carkeek_results <- data.frame(spp, disp, r2_conditional, r2_marginal,  beta_yr, se_yr, p_yr)

spp.dd <- park.props %>% filter(park == "Carkeek Park") %>% group_by(bird.code) %>% 
  summarise(nyrs_det = sum(prop > 0)) %>% filter(nyrs_det < 5 & nyrs_det > 0) %>% pull(bird.code)

spp.dd.df <- data.frame(spp = spp.dd,
                        disp = rep(NA, length(spp.dd)),
                        r2_conditional = rep(NA, length(spp.dd)),
                        r2_marginal = rep(NA, length(spp.dd)),
                        beta_yr = rep(NA, length(spp.dd)),
                        se_yr = rep(NA, length(spp.dd)),
                        p_yr = rep(NA, length(spp.dd))
                        )

carkeek_results <- carkeek_results %>% rbind(spp.dd.df) %>% 
  mutate(status = case_when(
    beta_yr >= 0.1 & p_yr < 0.05 & r2_marginal >= 0.2 ~ "Increasing",
    beta_yr < 0.1 & beta_yr > 0 & p_yr < 0.1 & r2_marginal >= 0.1 ~ "Possibly increasing",
    beta_yr <= -0.1 & p_yr < 0.05 & r2_marginal >= 0.2 ~ "Decreasing",
    beta_yr >= -0.1 & beta_yr < 0 & p_yr < 0.1 & r2_marginal >= 0.1 ~ "Possibly decreasing",
    is.na(beta_yr) ~ "Data deficient",
    TRUE ~ "Apparantly stable"
))


view(res_cat)


# DISCOVERY

discovery <- park.props %>% filter(park == "Discovery Park") %>% mutate(yr = as.numeric(scale(year)))

spp <- park.props %>% filter(park == "Discovery Park") %>% group_by(bird.code) %>% 
  summarise(nyrs_det = sum(prop > 0)) %>% filter(nyrs_det >= 5) %>% pull(bird.code)

disp <- numeric()
r2_conditional <- numeric()
r2_marginal <- numeric()
beta_yr <- numeric()
se_yr <- numeric()
p_yr <- numeric()

for(i in 1:length(spp)){
  mod.dat <- discovery[discovery$bird.code == spp[i], ]
  mod <- glmmTMB(cbind(detected, not_detected) ~ yr + (1 | station.code), data = mod.dat, family = binomial)
  r2_conditional[i] <- r2_nakagawa(mod)$R2_conditional
  r2_marginal[i] <- r2_nakagawa(mod)$R2_marginal
  disp[i] <- sum(residuals(mod, type = "pearson")^2) / nrow(mod.dat) - length(fixef(mod)$cond)
  beta_yr[i] <- summary(mod)$coefficients$cond["yr", 1]
  se_yr[i] <- summary(mod)$coefficients$cond["yr", 2]
  p_yr[i] <- summary(mod)$coefficients$cond["yr", 4]
}

discovery_results <- data.frame(spp, disp, r2_conditional, r2_marginal, beta_yr, se_yr, p_yr)

spp.dd <- park.props %>% filter(park == "Discovery Park") %>% group_by(bird.code) %>% 
  summarise(nyrs_det = sum(prop > 0)) %>% filter(nyrs_det < 5 & nyrs_det > 0) %>% pull(bird.code)

spp.dd.df <- data.frame(spp = spp.dd,
                        disp = rep(NA, length(spp.dd)),
                        r2_conditional = rep(NA, length(spp.dd)),
                        r2_marginal = rep(NA, length(spp.dd)),
                        beta_yr = rep(NA, length(spp.dd)),
                        se_yr = rep(NA, length(spp.dd)),
                        p_yr = rep(NA, length(spp.dd))
)

discovery_results <- discovery_results %>% rbind(spp.dd.df) %>% 
  mutate(status = case_when(
    beta_yr >= 0.1 & p_yr < 0.05 & r2_marginal >= 0.2 ~ "Increasing",
    beta_yr < 0.1 & beta_yr > 0 & p_yr < 0.1 & r2_marginal >= 0.1 ~ "Possibly increasing",
    beta_yr <= -0.1 & p_yr < 0.05 & r2_marginal >= 0.2 ~ "Decreasing",
    beta_yr >= -0.1 & beta_yr < 0 & p_yr < 0.1 & r2_marginal >= 0.1 ~ "Possibly decreasing",
    is.na(beta_yr) ~ "Data deficient",
    TRUE ~ "Apparantly stable"
  ))



view(discovery_results)

# GENESEE

genesee <- park.props %>% filter(park == "Genesee Park") %>% mutate(yr = as.numeric(scale(year)))

spp <- park.props %>% filter(park == "Genesee Park") %>% group_by(bird.code) %>% 
  summarise(nyrs_det = sum(prop > 0)) %>% filter(nyrs_det >= 5) %>% pull(bird.code)

disp <- numeric()
r2_conditional <- numeric()
r2_marginal <- numeric()
beta_yr <- numeric()
se_yr <- numeric()
p_yr <- numeric()

for(i in 1:length(spp)){
  mod.dat <- genesee[genesee$bird.code == spp[i], ]
  mod <- glmmTMB(cbind(detected, not_detected) ~ yr + (1 | station.code), data = mod.dat, family = binomial)
  r2_conditional[i] <- r2_nakagawa(mod)$R2_conditional
  r2_marginal[i] <- r2_nakagawa(mod)$R2_marginal
  disp[i] <- sum(residuals(mod, type = "pearson")^2) / nrow(mod.dat) - length(fixef(mod)$cond)
  beta_yr[i] <- summary(mod)$coefficients$cond["yr", 1]
  se_yr[i] <- summary(mod)$coefficients$cond["yr", 2]
  p_yr[i] <- summary(mod)$coefficients$cond["yr", 4]
}

genesee_results <- data.frame(spp, disp, r2_conditional, r2_marginal, beta_yr, se_yr, p_yr)

spp.dd <- park.props %>% filter(park == "Genesee Park") %>% group_by(bird.code) %>% 
  summarise(nyrs_det = sum(prop > 0)) %>% filter(nyrs_det < 5 & nyrs_det > 0) %>% pull(bird.code)

spp.dd.df <- data.frame(spp = spp.dd,
                        disp = rep(NA, length(spp.dd)),
                        r2_conditional = rep(NA, length(spp.dd)),
                        r2_marginal = rep(NA, length(spp.dd)),
                        beta_yr = rep(NA, length(spp.dd)),
                        se_yr = rep(NA, length(spp.dd)),
                        p_yr = rep(NA, length(spp.dd))
)

genesee_results <- genesee_results %>% rbind(spp.dd.df) %>% 
  mutate(status = case_when(
    beta_yr >= 0.1 & p_yr < 0.05 & r2_marginal >= 0.2 ~ "Increasing",
    beta_yr < 0.1 & beta_yr > 0 & p_yr < 0.1 & r2_marginal >= 0.1 ~ "Possibly increasing",
    beta_yr <= -0.1 & p_yr < 0.05 & r2_marginal >= 0.2 ~ "Decreasing",
    beta_yr >= -0.1 & beta_yr < 0 & p_yr < 0.1 & r2_marginal >= 0.1 ~ "Possibly decreasing",
    is.na(beta_yr) ~ "Data deficient",
    TRUE ~ "Apparantly stable"
  ))

view(genesee_results)
#view(genesee_results %>% filter(p_yr < 0.05))

# GOLDEN GARDENS

golden <- park.props %>% filter(park == "Golden Gardens Park") %>% mutate(yr = as.numeric(scale(year)))

spp <- park.props %>% filter(park == "Golden Gardens Park") %>% group_by(bird.code) %>% 
  summarise(nyrs_det = sum(prop > 0)) %>% filter(nyrs_det >= 5) %>% pull(bird.code)

disp <- numeric()
r2_conditional <- numeric()
r2_marginal <- numeric()
beta_yr <- numeric()
se_yr <- numeric()
p_yr <- numeric()

for(i in 1:length(spp)){
  mod.dat <- golden[golden$bird.code == spp[i], ]
  mod <- glmmTMB(cbind(detected, not_detected) ~ yr + (1 | station.code), data = mod.dat, family = binomial)
  r2_conditional[i] <- r2_nakagawa(mod)$R2_conditional
  r2_marginal[i] <- r2_nakagawa(mod)$R2_marginal
  disp[i] <- sum(residuals(mod, type = "pearson")^2) / nrow(mod.dat) - length(fixef(mod)$cond)
  beta_yr[i] <- summary(mod)$coefficients$cond["yr", 1]
  se_yr[i] <- summary(mod)$coefficients$cond["yr", 2]
  p_yr[i] <- summary(mod)$coefficients$cond["yr", 4]
}

golden_results <- data.frame(spp, disp, r2_conditional, r2_marginal, beta_yr, se_yr, p_yr)

spp.dd <- park.props %>% filter(park == "Golden Gardens Park") %>% group_by(bird.code) %>% 
  summarise(nyrs_det = sum(prop > 0)) %>% filter(nyrs_det < 5 & nyrs_det > 0) %>% pull(bird.code)

spp.dd.df <- data.frame(spp = spp.dd,
                        disp = rep(NA, length(spp.dd)),
                        r2_conditional = rep(NA, length(spp.dd)),
                        r2_marginal = rep(NA, length(spp.dd)),
                        beta_yr = rep(NA, length(spp.dd)),
                        se_yr = rep(NA, length(spp.dd)),
                        p_yr = rep(NA, length(spp.dd))
)

golden_results <- golden_results %>% rbind(spp.dd.df) %>% 
  mutate(status = case_when(
    beta_yr >= 0.1 & p_yr < 0.05 & r2_marginal >= 0.2 ~ "Increasing",
    beta_yr < 0.1 & beta_yr > 0 & p_yr < 0.1 & r2_marginal >= 0.1 ~ "Possibly increasing",
    beta_yr <= -0.1 & p_yr < 0.05 & r2_marginal >= 0.2 ~ "Decreasing",
    beta_yr >= -0.1 & beta_yr < 0 & p_yr < 0.1 & r2_marginal >= 0.1 ~ "Possibly decreasing",
    is.na(beta_yr) ~ "Data deficient",
    TRUE ~ "Apparantly stable"
  ))

view(golden_results)

# LAKE FOREST PARK

lfp <- park.props %>% filter(park == "Lake Forest Park") %>% mutate(yr = as.numeric(scale(year)))

spp <- park.props %>% filter(park == "Lake Forest Park") %>% group_by(bird.code) %>% 
  summarise(nyrs_det = sum(prop > 0)) %>% filter(nyrs_det >= 5) %>% pull(bird.code)

disp <- numeric()
r2_conditional <- numeric()
r2_marginal <- numeric()
beta_yr <- numeric()
se_yr <- numeric()
p_yr <- numeric()

for(i in 1:length(spp)){
  mod.dat <- lfp[lfp$bird.code == spp[i], ]
  mod <- glmmTMB(cbind(detected, not_detected) ~ yr, data = mod.dat, family = binomial)
  r2_conditional[i] <- NA
  r2_marginal[i] <- r2_tjur(mod)
  disp[i] <- sum(residuals(mod, type = "pearson")^2) / nrow(mod.dat) - length(fixef(mod)$cond)
  beta_yr[i] <- summary(mod)$coefficients$cond["yr", 1]
  se_yr[i] <- summary(mod)$coefficients$cond["yr", 2]
  p_yr[i] <- summary(mod)$coefficients$cond["yr", 4]
}


lake_forest_results <- data.frame(spp, disp, r2_conditional, r2_marginal, beta_yr, se_yr, p_yr)

spp.dd <- park.props %>% filter(park == "Lake Forest Park") %>% group_by(bird.code) %>% 
  summarise(nyrs_det = sum(prop > 0)) %>% filter(nyrs_det < 5 & nyrs_det > 0) %>% pull(bird.code)

spp.dd.df <- data.frame(spp = spp.dd,
                        disp = rep(NA, length(spp.dd)),
                        r2_conditional = rep(NA, length(spp.dd)),
                        r2_marginal = rep(NA, length(spp.dd)),
                        beta_yr = rep(NA, length(spp.dd)),
                        se_yr = rep(NA, length(spp.dd)),
                        p_yr = rep(NA, length(spp.dd))
)

lake_forest_results <- lake_forest_results %>% rbind(spp.dd.df) %>% 
  mutate(status = case_when(
    beta_yr >= 0.1 & p_yr < 0.05 & r2_marginal >= 0.2 ~ "Increasing",
    beta_yr < 0.1 & beta_yr > 0 & p_yr < 0.1 & r2_marginal >= 0.1 ~ "Possibly increasing",
    beta_yr <= -0.1 & p_yr < 0.05 & r2_marginal >= 0.2 ~ "Decreasing",
    beta_yr >= -0.1 & beta_yr < 0 & p_yr < 0.1 & r2_marginal >= 0.1 ~ "Possibly decreasing",
    is.na(beta_yr) ~ "Data deficient",
    TRUE ~ "Apparantly stable"
  ))

# LINCOLN PARK

linc <- park.props %>% filter(park == "Lincoln Park") %>% mutate(yr = as.numeric(scale(year)))

spp <- park.props %>% filter(park == "Lincoln Park") %>% group_by(bird.code) %>% 
  summarise(nyrs_det = sum(prop > 0)) %>% filter(nyrs_det >= 5) %>% pull(bird.code)

disp <- numeric()
r2_conditional <- numeric()
r2_marginal <- numeric()
beta_yr <- numeric()
se_yr <- numeric()
p_yr <- numeric()

for(i in 1:length(spp)){
  mod.dat <- linc[linc$bird.code == spp[i], ]
  mod <- glmmTMB(cbind(detected, not_detected) ~ yr + (1 | station.code), data = mod.dat, family = binomial)
  r2_conditional[i] <- r2_nakagawa(mod)$R2_conditional
  r2_marginal[i] <- r2_nakagawa(mod)$R2_marginal
  disp[i] <- sum(residuals(mod, type = "pearson")^2) / nrow(mod.dat) - length(fixef(mod)$cond)
  beta_yr[i] <- summary(mod)$coefficients$cond["yr", 1]
  se_yr[i] <- summary(mod)$coefficients$cond["yr", 2]
  p_yr[i] <- summary(mod)$coefficients$cond["yr", 4]
}

lincoln_results <- data.frame(spp, disp, r2_conditional, r2_marginal, beta_yr, se_yr, p_yr)

spp.dd <- park.props %>% filter(park == "Lincoln Park") %>% group_by(bird.code) %>% 
  summarise(nyrs_det = sum(prop > 0)) %>% filter(nyrs_det < 5 & nyrs_det > 0) %>% pull(bird.code)

spp.dd.df <- data.frame(spp = spp.dd,
                        disp = rep(NA, length(spp.dd)),
                        r2_conditional = rep(NA, length(spp.dd)),
                        r2_marginal = rep(NA, length(spp.dd)),
                        beta_yr = rep(NA, length(spp.dd)),
                        se_yr = rep(NA, length(spp.dd)),
                        p_yr = rep(NA, length(spp.dd))
)

lincoln_results <- lincoln_results %>% rbind(spp.dd.df) %>% 
  mutate(status = case_when(
    beta_yr >= 0.1 & p_yr < 0.05 & r2_marginal >= 0.2 ~ "Increasing",
    beta_yr < 0.1 & beta_yr > 0 & p_yr < 0.1 & r2_marginal >= 0.1 ~ "Possibly increasing",
    beta_yr <= -0.1 & p_yr < 0.05 & r2_marginal >= 0.2 ~ "Decreasing",
    beta_yr >= -0.1 & beta_yr < 0 & p_yr < 0.1 & r2_marginal >= 0.1 ~ "Possibly decreasing",
    is.na(beta_yr) ~ "Data deficient",
    TRUE ~ "Apparantly stable"
  ))


# MAGNUSON PARK

mag <- park.props %>% filter(park == "Magnuson Park") %>% mutate(yr = as.numeric(scale(year)))

spp <- park.props %>% filter(park == "Magnuson Park") %>% group_by(bird.code) %>% 
  summarise(nyrs_det = sum(prop > 0)) %>% filter(nyrs_det >= 5) %>% pull(bird.code)

disp <- numeric()
r2_conditional <- numeric()
r2_marginal <- numeric()
beta_yr <- numeric()
se_yr <- numeric()
p_yr <- numeric()

for(i in 1:length(spp)){
  mod.dat <- mag[mag$bird.code == spp[i], ]
  mod <- glmmTMB(cbind(detected, not_detected) ~ yr + (1 | station.code), data = mod.dat, family = binomial)
  r2_conditional[i] <- r2_nakagawa(mod)$R2_conditional
  r2_marginal[i] <- r2_nakagawa(mod)$R2_marginal
  disp[i] <- sum(residuals(mod, type = "pearson")^2) / nrow(mod.dat) - length(fixef(mod)$cond)
  beta_yr[i] <- summary(mod)$coefficients$cond["yr", 1]
  se_yr[i] <- summary(mod)$coefficients$cond["yr", 2]
  p_yr[i] <- summary(mod)$coefficients$cond["yr", 4]
}

magnuson_results <- data.frame(spp, disp, r2_conditional, r2_marginal, beta_yr, se_yr, p_yr)

spp.dd <- park.props %>% filter(park == "Magnuson Park") %>% group_by(bird.code) %>% 
  summarise(nyrs_det = sum(prop > 0)) %>% filter(nyrs_det < 5 & nyrs_det > 0) %>% pull(bird.code)

spp.dd.df <- data.frame(spp = spp.dd,
                        disp = rep(NA, length(spp.dd)),
                        r2_conditional = rep(NA, length(spp.dd)),
                        r2_marginal = rep(NA, length(spp.dd)),
                        beta_yr = rep(NA, length(spp.dd)),
                        se_yr = rep(NA, length(spp.dd)),
                        p_yr = rep(NA, length(spp.dd))
)

magnuson_results <- magnuson_results %>% rbind(spp.dd.df) %>% 
  mutate(status = case_when(
    beta_yr >= 0.1 & p_yr < 0.05 & r2_marginal >= 0.2 ~ "Increasing",
    beta_yr < 0.1 & beta_yr > 0 & p_yr < 0.1 & r2_marginal >= 0.1 ~ "Possibly increasing",
    beta_yr <= -0.1 & p_yr < 0.05 & r2_marginal >= 0.2 ~ "Decreasing",
    beta_yr >= -0.1 & beta_yr < 0 & p_yr < 0.1 & r2_marginal >= 0.1 ~ "Possibly decreasing",
    is.na(beta_yr) ~ "Data deficient",
    TRUE ~ "Apparantly stable"
  ))

# SEWARD PARK

seward <- park.props %>% filter(park == "Seward Park") %>% mutate(yr = as.numeric(scale(year)))

spp <- park.props %>% filter(park == "Seward Park") %>% group_by(bird.code) %>% 
  summarise(nyrs_det = sum(prop > 0)) %>% filter(nyrs_det >= 5) %>% pull(bird.code)

disp <- numeric()
r2_conditional <- numeric()
r2_marginal <- numeric()
beta_yr <- numeric()
se_yr <- numeric()
p_yr <- numeric()

for(i in 1:length(spp)){
  mod.dat <- seward[seward$bird.code == spp[i], ]
  mod <- glmmTMB(cbind(detected, not_detected) ~ yr + (1 | station.code), data = mod.dat, family = binomial)
  r2_conditional[i] <- r2_nakagawa(mod)$R2_conditional
  r2_marginal[i] <- r2_nakagawa(mod)$R2_marginal
  disp[i] <- sum(residuals(mod, type = "pearson")^2) / nrow(mod.dat) - length(fixef(mod)$cond)
  beta_yr[i] <- summary(mod)$coefficients$cond["yr", 1]
  se_yr[i] <- summary(mod)$coefficients$cond["yr", 2]
  p_yr[i] <- summary(mod)$coefficients$cond["yr", 4]
}

seward_results <- data.frame(spp, disp, r2_conditional, r2_marginal, beta_yr, se_yr, p_yr)

spp.dd <- park.props %>% filter(park == "Seward Park") %>% group_by(bird.code) %>% 
  summarise(nyrs_det = sum(prop > 0)) %>% filter(nyrs_det < 5 & nyrs_det > 0) %>% pull(bird.code)

spp.dd.df <- data.frame(spp = spp.dd,
                        disp = rep(NA, length(spp.dd)),
                        r2_conditional = rep(NA, length(spp.dd)),
                        r2_marginal = rep(NA, length(spp.dd)),
                        beta_yr = rep(NA, length(spp.dd)),
                        se_yr = rep(NA, length(spp.dd)),
                        p_yr = rep(NA, length(spp.dd))
)

seward_results <- seward_results %>% rbind(spp.dd.df) %>% 
  mutate(status = case_when(
    beta_yr >= 0.1 & p_yr < 0.05 & r2_marginal >= 0.2 ~ "Increasing",
    beta_yr < 0.1 & beta_yr > 0 & p_yr < 0.1 & r2_marginal >= 0.1 ~ "Possibly increasing",
    beta_yr <= -0.1 & p_yr < 0.05 & r2_marginal >= 0.2 ~ "Decreasing",
    beta_yr >= -0.1 & beta_yr < 0 & p_yr < 0.1 & r2_marginal >= 0.1 ~ "Possibly decreasing",
    is.na(beta_yr) ~ "Data deficient",
    TRUE ~ "Apparantly stable"
  ))

# WASHINGTON PARK ARBORETUM

arb <- park.props %>% filter(park == "Washington Park Arboretum") %>% mutate(yr = as.numeric(scale(year)))

spp <- park.props %>% filter(park == "Washington Park Arboretum") %>% group_by(bird.code) %>% 
  summarise(nyrs_det = sum(prop > 0)) %>% filter(nyrs_det >= 5) %>% pull(bird.code)

disp <- numeric()
r2_conditional <- numeric()
r2_marginal <- numeric()
beta_yr <- numeric()
se_yr <- numeric()
p_yr <- numeric()

for(i in 1:length(spp)){
  mod.dat <- arb[arb$bird.code == spp[i], ]
  mod <- glmmTMB(cbind(detected, not_detected) ~ yr + (1 | station.code), data = mod.dat, family = binomial)
  r2_conditional[i] <- r2_nakagawa(mod)$R2_conditional
  r2_marginal[i] <- r2_nakagawa(mod)$R2_marginal
  disp[i] <- sum(residuals(mod, type = "pearson")^2) / nrow(mod.dat) - length(fixef(mod)$cond)
  beta_yr[i] <- summary(mod)$coefficients$cond["yr", 1]
  se_yr[i] <- summary(mod)$coefficients$cond["yr", 2]
  p_yr[i] <- summary(mod)$coefficients$cond["yr", 4]
}

arb_results <- data.frame(spp, disp, r2_conditional, r2_marginal, beta_yr, se_yr, p_yr)

spp.dd <- park.props %>% filter(park == "Washington Park Arboretum") %>% group_by(bird.code) %>% 
  summarise(nyrs_det = sum(prop > 0)) %>% filter(nyrs_det < 5 & nyrs_det > 0) %>% pull(bird.code)

spp.dd.df <- data.frame(spp = spp.dd,
                        disp = rep(NA, length(spp.dd)),
                        r2_conditional = rep(NA, length(spp.dd)),
                        r2_marginal = rep(NA, length(spp.dd)),
                        beta_yr = rep(NA, length(spp.dd)),
                        se_yr = rep(NA, length(spp.dd)),
                        p_yr = rep(NA, length(spp.dd))
)

arb_results <- arb_results %>% rbind(spp.dd.df) %>% 
  mutate(status = case_when(
    beta_yr >= 0.1 & p_yr < 0.05 & r2_marginal >= 0.2 ~ "Increasing",
    beta_yr < 0.1 & beta_yr > 0 & p_yr < 0.1 & r2_marginal >= 0.1 ~ "Possibly increasing",
    beta_yr <= -0.1 & p_yr < 0.05 & r2_marginal >= 0.2 ~ "Decreasing",
    beta_yr >= -0.1 & beta_yr < 0 & p_yr < 0.1 & r2_marginal >= 0.1 ~ "Possibly decreasing",
    is.na(beta_yr) ~ "Data deficient",
    TRUE ~ "Apparantly stable"
  ))

res_cat$park <- rep("Carkeek Park", dim(res_cat)[1])
discovery_results$park <- rep("Discovery Park", dim(discovery_results)[1])
genesee_results$park <- rep("Genesee Park", dim(genesee_results)[1])
golden_results$park <- rep("Golden Gardens Park", dim(golden_results)[1])
lake_forest_results$park <- rep("Lake Forest Park", dim(lake_forest_results)[1])
lincoln_results$park <- rep("Lincoln Park", dim(lincoln_results)[1])

magnuson_results$park <- rep("Magnuson Park", dim(magnuson_results)[1])
seward_results$park <- rep("Seward Park", dim(seward_results)[1])
arb_results$park <- rep("Washington Park Arboretum", dim(arb_results)[1])

trends_by_park <- bind_rows(res_cat, discovery_results, genesee_results, golden_results,
                            lake_forest_results, lincoln_results, magnuson_results, seward_results,
                            arb_results) %>%
  left_join(., traits, join_by("spp" == "bird.code"))

write.csv(trends_by_park, "results/tables/trends_by_park.csv", row.names = FALSE)

filter(trends_by_park, p_yr < 0.1, beta_yr < -0.05, r2_marginal > 0.15)


trends_by_park <- trends_by_park %>% mutate(status2 = case_when(
  beta_yr > 0.05 & p_yr < 0.1 & r2_marginal > 0.1 ~ "Increasing",
  beta_yr < 0.05 & beta_yr > 0 & p_yr < 0.1 & r2_marginal >= 0.1 ~ "Possibly increasing",
  beta_yr < -0.05 & p_yr < 0.1 & r2_marginal >= 0.1 ~ "Decreasing",
  beta_yr > -0.05 & beta_yr < 0 & p_yr < 0.1 & r2_marginal >= 0.1 ~ "Possibly decreasing",
  is.na(beta_yr) ~ "Data deficient",
  TRUE ~ "Apparantly stable"
))




donut.dat <- trends_by_park %>% group_by(park, status2) %>% summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = status2, values_from = n, values_fill = list(n = 0)) %>%
  pivot_longer(-1, names_to = "status2", values_to = "n") %>%
  mutate(status2 = factor(status2, levels = c( "Data deficient", "Decreasing", "Apparantly stable", "Increasing")))

donut_pal <- c(
  "Data deficient" = "lightgrey",
  "Decreasing" = "#FFB98C",
  "Apparantly stable" = "#B2E632",
  "Increasing" = "#36BA3A"
)

ggplot(donut.dat %>% filter(park == "Carkeek Park"), aes(x = "", y = n, fill = status2)) + 
  geom_bar(stat ="identity", width = 0.2) +
  scale_fill_manual(values = donut_pal) +
  coord_flip() +
  theme_void() +
  theme(legend.position = 'none')

donut.dat %>% filter(park == "Carkeek Park")
