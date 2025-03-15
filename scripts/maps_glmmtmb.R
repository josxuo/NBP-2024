## Modeling p / a 

# libraries
library(rjags)
load.module("glm")
library(tidyverse)
library(readxl)
library(glmmTMB)
library(DHARMa)


# functions

summarise<-function(x){
  if(!is.matrix(x)){ x<-as.matrix(x) }
  low<-apply(x,2,quantile,probs=0.025)
  upp<-apply(x,2,quantile,probs=0.975)
  mean<-apply(x,2,mean)
  out<-cbind(mean,low,upp)
  colnames(out)<-c("mean","lower","upper")
  signif(out,digits=4)
}


# data
nbp <- read_xlsx("data/b_intermediate_data/nbp_tidy_jan_24.xlsx")
covs <- read_csv("data/c_analysis_data/circ_no_overlap_covariates.csv")

covs[is.na(covs)] <- 0

covs$scaleDowntown <- as.numeric(scale(covs$mDowntown))
covs$scaleWater <- as.numeric(scale(covs$mMajorWater))
covs$scaleCanopyProp2016 <- as.numeric(scale(covs$CanopyProp2016))
covs$scaleCanopyProp2021 <- as.numeric(scale(covs$CanopyProp2021))

# station / park table
station.table <- nbp %>%
  select(station.code, park) %>%
  unique()


d <- nbp %>%
  filter(!str_detect(species, pattern = " sp\\."),
         station.code %in% covs$Station, 
         year %in% c(2005:2019, 2022, 2023)) %>%
         #exclusions == "none") %>%
  mutate(loop2 = as.factor(paste(park, loop, sep = "-")),
         park = as.factor(park),
         month = as.factor(month),
         station.code = as.factor(station.code),
         station.year = paste(station.code, year, sep = "-"),
         park.year = paste(park, year, sep = "-"))


# inspect data
## sort(unique(d$species))
## sort(unique(d$year))
## sort(unique(d$code))


# Set up tables/calculations useful later
## number of surveys conducted each year
nsurvYearsp <- d %>%
  group_by(year, park, station.code) %>%
  reframe(nsurv = n_distinct(survey_id)) %>%
  mutate(station.year = paste(station.code, year, sep = "-"),
         park.year = paste(park, year, sep = "-")) %>%
  select(station.year, park.year, nsurv)
nsurvYearsp



## species p/a matrix
presence.site <- d %>%
  # filter(station.code %in% sus, year %in% year, month %in% months) %>%
  group_by(year, bird.code, park, loop2, station.code, station.year) %>%
  reframe(detected = n()) %>%
  pivot_wider(names_from = bird.code, values_from = detected, values_fill = 0) %>%
  pivot_longer(-grep("year|park|loop|station", colnames(.)), names_to = "bird.code", values_to = "detected") %>%
  left_join(., nsurvYearsp, join_by("station.year" == "station.year")) %>%
  mutate(not_detected = nsurv - detected)


focal <- "DEJU"
focal.dat.site <- presence.site %>% filter(bird.code == focal) %>% 
  left_join(., covs, join_by("station.code" == "Station")) %>%
  mutate(station.code = factor(station.code),
         yr = year - min(year))

view(focal.dat.site)

mod.null <- glmmTMB(cbind(detected, not_detected) ~ 1,
                family = binomial(link = "logit"), 
                data = focal.dat)

mod1 <- glmmTMB(cbind(detected, not_detected) ~ yr,
                family = binomial(link = "logit"), 
                data = focal.dat)

mod2 <- glmmTMB(cbind(detected, not_detected) ~ yr + park,
                family = binomial(link = "logit"), 
                data = focal.dat)

mod3 <- glmmTMB(cbind(detected, not_detected) ~ yr + park + (1 | station.code),
                family = binomial(link = "logit"),
                data = focal.dat)

mod4 <- glmmTMB(cbind(detected, not_detected) ~ yr + park + (1 | station.code) + (1| yr),
        family = binomial(link = "logit"),
        data = focal.dat)

mod5 <- glmmTMB(cbind(detected, not_detected) ~ yr + park + (yr | station.code),
        family = binomial(link = "logit"),
        data = focal.dat)

mod6 <- glmmTMB(cbind(detected, not_detected) ~ yr + (yr | station.code),
                family = binomial(link = "logit"),
                data = focal.dat)

summary(mod5)

AIC(mod.null, mod1, mod2, mod3, mod4, mod5, mod6)  ## most support for model 5. check diagnostics

dispersion_test <- simulateResiduals(mod5)
plot(dispersion_test)

testZeroInflation(dispersion_test)

## OK this model works; diagnostics are OK; overdispersion not an issue; zero inflation not an issue

summary(mod5)

# Coefficient for yr
intercept <- -4.278955
coef_yr <- 0.150302

# Function to convert log-odds to probability
log_odds_to_prob <- function(log_odds) {
  return(1 / (1 + exp(-log_odds)))
}

# Probabilities for different years
probability_yr0 <- log_odds_to_prob(intercept + coef_yr * 0)  # yr = 0
probability_yr1 <- log_odds_to_prob(intercept + coef_yr * 1)  # yr = 1
probability_yr2 <- log_odds_to_prob(intercept + coef_yr * 2)  # yr = 2
probability_yr3 <- log_odds_to_prob(intercept + coef_yr * 3)

probability_yr20 <- log_odds_to_prob(intercept + coef_yr * 20)
probability_yr21 <- log_odds_to_prob(intercept + coef_yr * 21)


# Results
probability_yr0
probability_yr1
probability_yr2


(probability_yr3 - probability_yr2) / probability_yr2
(probability_yr2 - probability_yr1) / probability_yr1
(probability_yr1 - probability_yr0) / probability_yr0

probability_yr21/probability_yr20
# Create a new data frame for prediction
predict_data <- expand.grid(
  yr = seq(min(focal.dat$yr), max(focal.dat$yr), length.out = 100),  # time range for prediction
  station.code = levels(focal.dat$station.code)) %>%
  left_join(., station.table, join_by("station.code" == "station.code"))

# Create predictions (including fixed and random effects)
predicted <- predict(mod5, newdata = predict_data, type = "response", se.fit = TRUE)

predicted_df <- cbind(predict_data, predicted[[1]], predicted[[2]]) %>%
  rename(mu = `predicted[[1]]`, se = `predicted[[2]]`) %>%
  mutate(lower = mu - (1.96 * se),
         upper = mu + (1.96 * se))

mean_predicted <- predicted_df %>%
  group_by(yr) %>%
  reframe(mean_mu = mean(mu),
          mean_lower = mean(lower),
          mean_upper = mean(upper))

focal_mean <- focal.dat %>%
  mutate(prop = detected / nsurv) %>%
  group_by(yr) %>%
  reframe(mean_prop = mean(prop))


ggplot(focal_mean, aes(x = (yr + 2005), y = mean_prop)) +
  geom_point() + 
  geom_line(data = mean_predicted, aes(x = (yr + 2005), y = mean_mu)) +
  geom_line(data = mean_predicted, aes(x = (yr + 2005), y = mean_lower), linetype = "dashed") +
  geom_line(data = mean_predicted, aes(x = (yr + 2005), y = mean_upper), linetype = "dashed") + 
  labs(title = "Proportion of surveys reporting Dark-eyed Junco over time", x = "Year", y = "Proportion of surveys reporting") +
  theme_minimal()


## check for autocorrelation
n <- dim(focal.dat)[1]
res1<-residuals(simulateResiduals(mod5,n=nsim,seed=F))
plot(res1[1:(n-1)],res1[2:n],xlab="residual at year t",ylab="residual at year t+1")


## No apparent temporal autocorrelation

######################## TRY FOR ANNA'S HUMMINGBIRD ########################
focal <- "ANHU"
focal.dat <- presence %>% filter(bird.code == focal) %>%
                                 #!station.year %in% outliers_station.year) %>% 
  left_join(., covs, join_by("station.code" == "Station")) %>%
  mutate(station.code = factor(station.code),
         yr = year - min(year))


mod.null <- glmmTMB(cbind(detected, not_detected) ~ 1,
                    family = binomial(link = "logit"), 
                    data = focal.dat)

mod1 <- glmmTMB(cbind(detected, not_detected) ~ yr,
                family = binomial(link = "logit"), 
                data = focal.dat)

mod2 <- glmmTMB(cbind(detected, not_detected) ~ yr + park,
                family = binomial(link = "logit"), 
                data = focal.dat)

mod3 <- glmmTMB(cbind(detected, not_detected) ~ yr + park + (1 | station.code),
                family = binomial(link = "logit"),
                data = focal.dat)

mod4 <- glmmTMB(cbind(detected, not_detected) ~ yr + park + (1 | station.code) + (1| yr),
                family = binomial(link = "logit"),
                data = focal.dat)

mod5 <- glmmTMB(cbind(detected, not_detected) ~ yr + park + (yr | station.code),
                family = binomial(link = "logit"),
                data = focal.dat)

mod6 <- glmmTMB(cbind(detected, not_detected) ~ yr + (yr | station.code),
                family = binomial(link = "logit"),
                data = focal.dat)

AIC(mod.null, mod1, mod2, mod3, mod4, mod5, mod6)  ## most support for model 5. check diagnostics

dispersion_test <- simulateResiduals(mod5)
plot(dispersion_test)

testZeroInflation(dispersion_test)

testOutliers(dispersion_test)

outliers_test <- outliers(dispersion_test)
outliers_data <- focal.dat[outliers_test, ]
outliers_station.year <- focal.dat[outliers_test, ]$station.year


## overdispersion not a problem; excluded a few outliers to avoid dragging estimates 
## high from a few high high count years

summary(mod5)


# Coefficient for yr
intercept <- -1.877202
coef_yr <- 0.067079

# Function to convert log-odds to probability
log_odds_to_prob <- function(log_odds) {
  return(1 / (1 + exp(-log_odds)))
}

# Probabilities for different years
probability_yr0 <- log_odds_to_prob(intercept + coef_yr * 0)  # yr = 0
probability_yr1 <- log_odds_to_prob(intercept + coef_yr * 1)  # yr = 1
probability_yr2 <- log_odds_to_prob(intercept + coef_yr * 2)  # yr = 2
probability_yr3 <- log_odds_to_prob(intercept + coef_yr * 3)

probability_yr20 <- log_odds_to_prob(intercept + coef_yr * 20)
probability_yr21 <- log_odds_to_prob(intercept + coef_yr * 21)


# Results
probability_yr0
probability_yr1
probability_yr2


(probability_yr3 - probability_yr2) / probability_yr2
(probability_yr2 - probability_yr1) / probability_yr1
(probability_yr1 - probability_yr0) / probability_yr0

probability_yr21/probability_yr20
# Create a new data frame for prediction
predict_data <- expand.grid(
  yr = seq(min(focal.dat$yr), max(focal.dat$yr), length.out = 100),  # time range for prediction
  station.code = levels(focal.dat$station.code)) %>%
  left_join(., station.table, join_by("station.code" == "station.code"))

# Create predictions (including fixed and random effects)
predicted <- predict(mod5, newdata = predict_data, type = "response", se.fit = TRUE)

predicted_df <- cbind(predict_data, predicted[[1]], predicted[[2]]) %>%
  rename(mu = `predicted[[1]]`, se = `predicted[[2]]`) %>%
  mutate(lower = mu - (1.96 * se),
         upper = mu + (1.96 * se))

mean_predicted <- predicted_df %>%
  group_by(yr) %>%
  reframe(mean_mu = mean(mu),
          mean_lower = mean(lower),
          mean_upper = mean(upper))

focal_mean <- focal.dat %>%
  mutate(prop = detected / nsurv) %>%
  group_by(yr) %>%
  reframe(mean_prop = mean(prop))


ggplot(focal_mean, aes(x = (yr + 2005), y = mean_prop)) +
  geom_point() + 
  geom_line(data = mean_predicted, aes(x = (yr + 2005), y = mean_mu)) +
  geom_line(data = mean_predicted, aes(x = (yr + 2005), y = mean_lower), linetype = "dashed") +
  geom_line(data = mean_predicted, aes(x = (yr + 2005), y = mean_upper), linetype = "dashed") + 
  labs(title = "Proportion of surveys reporting Anna's Hummingbird over time", x = "Year", y = "Proportion of surveys reporting") +
  theme_minimal()


## check for autocorrelation
n <- dim(focal.dat)[1]
res1<-residuals(simulateResiduals(mod5,n=nsim,seed=F))
plot(res1[1:(n-1)],res1[2:n],xlab="residual at year t",ylab="residual at year t+1")


######################## TRY FOR AMERICAN CROW ########################
focal <- "AMCR"
focal.dat <- presence %>% filter(bird.code == focal) %>%
  #!station.year %in% outliers_station.year) %>% 
  left_join(., covs, join_by("station.code" == "Station")) %>%
  mutate(station.code = factor(station.code),
         yr = year - min(year))


mod.null <- glmmTMB(cbind(detected, not_detected) ~ 1,
                    family = binomial(link = "logit"), 
                    data = focal.dat)

mod1 <- glmmTMB(cbind(detected, not_detected) ~ yr,
                family = binomial(link = "logit"), 
                data = focal.dat)

mod2 <- glmmTMB(cbind(detected, not_detected) ~ yr + park,
                family = binomial(link = "logit"), 
                data = focal.dat)

mod3 <- glmmTMB(cbind(detected, not_detected) ~ yr + park + (1 | station.code),
                family = binomial(link = "logit"),
                data = focal.dat)

mod4 <- glmmTMB(cbind(detected, not_detected) ~ yr + park + (1 | station.code) + (1| yr),
                family = binomial(link = "logit"),
                data = focal.dat)

mod5 <- glmmTMB(cbind(detected, not_detected) ~ yr + park + (yr | station.code),
                family = binomial(link = "logit"),
                data = focal.dat)

mod6 <- glmmTMB(cbind(detected, not_detected) ~ yr + (yr | station.code),
                family = binomial(link = "logit"),
                data = focal.dat)

AIC(mod.null, mod1, mod2, mod3, mod4, mod5, mod6)  ## most support for model 5. check diagnostics

dispersion_test <- simulateResiduals(mod5)
plot(dispersion_test)

testZeroInflation(dispersion_test)

testOutliers(dispersion_test)

outliers_test <- outliers(dispersion_test)
outliers_data <- focal.dat[outliers_test, ]
outliers_station.year <- focal.dat[outliers_test, ]$station.year


## overdispersion not a problem; excluded a few outliers to avoid dragging estimates 
## high from a few high high count years

summary(mod5)


# Coefficient for yr
intercept <- -1.877202
coef_yr <- 0.067079

# Function to convert log-odds to probability
log_odds_to_prob <- function(log_odds) {
  return(1 / (1 + exp(-log_odds)))
}

# Probabilities for different years
probability_yr0 <- log_odds_to_prob(intercept + coef_yr * 0)  # yr = 0
probability_yr1 <- log_odds_to_prob(intercept + coef_yr * 1)  # yr = 1
probability_yr2 <- log_odds_to_prob(intercept + coef_yr * 2)  # yr = 2
probability_yr3 <- log_odds_to_prob(intercept + coef_yr * 3)

probability_yr20 <- log_odds_to_prob(intercept + coef_yr * 20)
probability_yr21 <- log_odds_to_prob(intercept + coef_yr * 21)


# Results
probability_yr0
probability_yr1
probability_yr2


(probability_yr3 - probability_yr2) / probability_yr2
(probability_yr2 - probability_yr1) / probability_yr1
(probability_yr1 - probability_yr0) / probability_yr0

probability_yr21/probability_yr20
# Create a new data frame for prediction
predict_data <- expand.grid(
  yr = seq(min(focal.dat$yr), max(focal.dat$yr), length.out = 100),  # time range for prediction
  station.code = levels(focal.dat$station.code)) %>%
  left_join(., station.table, join_by("station.code" == "station.code"))

# Create predictions (including fixed and random effects)
predicted <- predict(mod5, newdata = predict_data, type = "response", se.fit = TRUE)

predicted_df <- cbind(predict_data, predicted[[1]], predicted[[2]]) %>%
  rename(mu = `predicted[[1]]`, se = `predicted[[2]]`) %>%
  mutate(lower = mu - (1.96 * se),
         upper = mu + (1.96 * se))

mean_predicted <- predicted_df %>%
  group_by(yr) %>%
  reframe(mean_mu = mean(mu),
          mean_lower = mean(lower),
          mean_upper = mean(upper))

focal_mean <- focal.dat %>%
  mutate(prop = detected / nsurv) %>%
  group_by(yr) %>%
  reframe(mean_prop = mean(prop))


ggplot(focal_mean, aes(x = (yr + 2005), y = mean_prop)) +
  geom_point() + 
  geom_line(data = mean_predicted, aes(x = (yr + 2005), y = mean_mu)) +
  geom_line(data = mean_predicted, aes(x = (yr + 2005), y = mean_lower), linetype = "dashed") +
  geom_line(data = mean_predicted, aes(x = (yr + 2005), y = mean_upper), linetype = "dashed") + 
  labs(title = "Proportion of surveys reporting Anna's Hummingbird over time", x = "Year", y = "Proportion of surveys reporting") +
  theme_minimal()


## check for autocorrelation
n <- dim(focal.dat)[1]
res1<-residuals(simulateResiduals(mod5,n=nsim,seed=F))
plot(res1[1:(n-1)],res1[2:n],xlab="residual at year t",ylab="residual at year t+1")

######################## TRY FOR BALD EAGLE ########################
focal <- "BAEA"
focal.dat <- presence %>% filter(bird.code == focal) %>%
  #!station.year %in% outliers_station.year) %>% 
  left_join(., covs, join_by("station.code" == "Station")) %>%
  mutate(station.code = factor(station.code),
         yr = year - min(year))


mod.null <- glmmTMB(cbind(detected, not_detected) ~ 1,
                    family = binomial(link = "logit"), 
                    data = focal.dat)

mod1 <- glmmTMB(cbind(detected, not_detected) ~ yr,
                family = binomial(link = "logit"), 
                data = focal.dat)

mod2 <- glmmTMB(cbind(detected, not_detected) ~ yr + park,
                family = binomial(link = "logit"), 
                data = focal.dat)

mod3 <- glmmTMB(cbind(detected, not_detected) ~ yr + park + (1 | station.code),
                family = binomial(link = "logit"),
                data = focal.dat)

mod4 <- glmmTMB(cbind(detected, not_detected) ~ yr + park + (1 | station.code) + (1| yr),
                family = binomial(link = "logit"),
                data = focal.dat)

mod5 <- glmmTMB(cbind(detected, not_detected) ~ yr + park + (yr | station.code),
                family = binomial(link = "logit"),
                data = focal.dat)

mod6 <- glmmTMB(cbind(detected, not_detected) ~ yr + (yr | station.code),
                family = binomial(link = "logit"),
                data = focal.dat)

AIC(mod.null, mod1, mod2, mod3, mod4, mod5, mod6)  ## most support for model 5. check diagnostics

dispersion_test <- simulateResiduals(mod5)
plot(dispersion_test)

testZeroInflation(dispersion_test)

testOutliers(dispersion_test)

outliers_test <- outliers(dispersion_test)
outliers_data <- focal.dat[outliers_test, ]
outliers_station.year <- focal.dat[outliers_test, ]$station.year


## overdispersion not a problem; excluded a few outliers to avoid dragging estimates 
## high from a few high high count years

summary(mod5)


# Coefficient for yr
intercept <- -1.877202
coef_yr <- 0.067079

# Function to convert log-odds to probability
log_odds_to_prob <- function(log_odds) {
  return(1 / (1 + exp(-log_odds)))
}

# Probabilities for different years
probability_yr0 <- log_odds_to_prob(intercept + coef_yr * 0)  # yr = 0
probability_yr1 <- log_odds_to_prob(intercept + coef_yr * 1)  # yr = 1
probability_yr2 <- log_odds_to_prob(intercept + coef_yr * 2)  # yr = 2
probability_yr3 <- log_odds_to_prob(intercept + coef_yr * 3)

probability_yr20 <- log_odds_to_prob(intercept + coef_yr * 20)
probability_yr21 <- log_odds_to_prob(intercept + coef_yr * 21)


# Results
probability_yr0
probability_yr1
probability_yr2


(probability_yr3 - probability_yr2) / probability_yr2
(probability_yr2 - probability_yr1) / probability_yr1
(probability_yr1 - probability_yr0) / probability_yr0

probability_yr21/probability_yr20
# Create a new data frame for prediction
predict_data <- expand.grid(
  yr = seq(min(focal.dat$yr), max(focal.dat$yr), length.out = 100),  # time range for prediction
  station.code = levels(focal.dat$station.code)) %>%
  left_join(., station.table, join_by("station.code" == "station.code"))

# Create predictions (including fixed and random effects)
predicted <- predict(mod5, newdata = predict_data, type = "response", se.fit = TRUE)

predicted_df <- cbind(predict_data, predicted[[1]], predicted[[2]]) %>%
  rename(mu = `predicted[[1]]`, se = `predicted[[2]]`) %>%
  mutate(lower = mu - (1.96 * se),
         upper = mu + (1.96 * se))

mean_predicted <- predicted_df %>%
  group_by(yr) %>%
  reframe(mean_mu = mean(mu),
          mean_lower = mean(lower),
          mean_upper = mean(upper))

focal_mean <- focal.dat %>%
  mutate(prop = detected / nsurv) %>%
  group_by(yr) %>%
  reframe(mean_prop = mean(prop))


ggplot(focal_mean, aes(x = (yr + 2005), y = mean_prop)) +
  geom_point() + 
  geom_line(data = mean_predicted, aes(x = (yr + 2005), y = mean_mu)) +
  geom_line(data = mean_predicted, aes(x = (yr + 2005), y = mean_lower), linetype = "dashed") +
  geom_line(data = mean_predicted, aes(x = (yr + 2005), y = mean_upper), linetype = "dashed") + 
  labs(title = "Proportion of surveys reporting Bald Eagle over time", x = "Year", y = "Proportion of surveys reporting") +
  theme_minimal()


## check for autocorrelation
n <- dim(focal.dat)[1]
res1<-residuals(simulateResiduals(mod5,n=nsim,seed=F))
plot(res1[1:(n-1)],res1[2:n],xlab="residual at year t",ylab="residual at year t+1")

######################## TRY FOR SPOTTED TOWHEE ########################
focal <- "STJA"
focal.dat <- presence %>% filter(bird.code == focal) %>%
  #!station.year %in% outliers_station.year) %>% 
  left_join(., covs, join_by("station.code" == "Station")) %>%
  mutate(station.code = factor(station.code),
         yr = year - min(year))


mod.null <- glmmTMB(cbind(detected, not_detected) ~ 1,
                    family = binomial(link = "logit"), 
                    data = focal.dat)

mod1 <- glmmTMB(cbind(detected, not_detected) ~ yr,
                family = binomial(link = "logit"), 
                data = focal.dat)

mod2 <- glmmTMB(cbind(detected, not_detected) ~ yr + park,
                family = binomial(link = "logit"), 
                data = focal.dat)

mod3 <- glmmTMB(cbind(detected, not_detected) ~ yr + park + (1 | station.code),
                family = binomial(link = "logit"),
                data = focal.dat)

mod4 <- glmmTMB(cbind(detected, not_detected) ~ yr + park + (1 | station.code) + (1| yr),
                family = binomial(link = "logit"),
                data = focal.dat)

mod5 <- glmmTMB(cbind(detected, not_detected) ~ yr + park + (yr | station.code),
                family = binomial(link = "logit"),
                data = focal.dat)

mod6 <- glmmTMB(cbind(detected, not_detected) ~ yr + (yr | station.code),
                family = binomial(link = "logit"),
                data = focal.dat)

AIC(mod.null, mod1, mod2, mod3, mod4, mod5, mod6)  ## most support for model 5. check diagnostics

dispersion_test <- simulateResiduals(mod4)
plot(dispersion_test)

testZeroInflation(dispersion_test)

testOutliers(dispersion_test)

outliers_test <- outliers(dispersion_test)
outliers_data <- focal.dat[outliers_test, ]
outliers_station.year <- focal.dat[outliers_test, ]$station.year


## overdispersion not a problem; 


summary(mod5)

# Create a new data frame for prediction
predict_data <- expand.grid(
  yr = seq(min(focal.dat$yr), max(focal.dat$yr), length.out = 100),  # time range for prediction
  station.code = levels(focal.dat$station.code)) %>%
  left_join(., station.table, join_by("station.code" == "station.code"))

# Create predictions (including fixed and random effects)
predicted <- predict(mod4, newdata = predict_data, type = "response", se.fit = TRUE, allow.new.levels = TRUE)

predicted_df <- cbind(predict_data, predicted[[1]], predicted[[2]]) %>%
  rename(mu = `predicted[[1]]`, se = `predicted[[2]]`) %>%
  mutate(lower = mu - (1.96 * se),
         upper = mu + (1.96 * se))

mean_predicted <- predicted_df %>%
  group_by(yr) %>%
  reframe(mean_mu = mean(mu),
          mean_lower = mean(lower),
          mean_upper = mean(upper))

focal_mean <- focal.dat %>%
  mutate(prop = detected / nsurv) %>%
  group_by(yr) %>%
  reframe(mean_prop = mean(prop))


ggplot(focal_mean, aes(x = (yr + 2005), y = mean_prop)) +
  geom_point() + 
  geom_line(data = mean_predicted, aes(x = (yr + 2005), y = mean_mu)) +
  geom_line(data = mean_predicted, aes(x = (yr + 2005), y = mean_lower), linetype = "dashed") +
  geom_line(data = mean_predicted, aes(x = (yr + 2005), y = mean_upper), linetype = "dashed") + 
  labs(title = "Proportion of surveys reporting Steller's Jay over time", x = "Year", y = "Proportion of surveys reporting") +
  theme_minimal()


## check for autocorrelation
n <- dim(focal.dat)[1]
res1<-residuals(simulateResiduals(mod5,n=nsim,seed=F))
plot(res1[1:(n-1)],res1[2:n],xlab="residual at year t",ylab="residual at year t+1")