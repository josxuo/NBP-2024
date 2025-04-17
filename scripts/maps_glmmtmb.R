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
nbp <- read_xlsx("data/processed/nbp_tidy_jan_24.xlsx")
covs <- read_csv("data/processed/circ_no_overlap_covariates.csv")

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
         station.code %in% covs$station.code, 
         year %in% c(2005:2019, 2021, 2022, 2023),
         month %in% c(5, 6)) %>%
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

focal.species.name <- "Rufous Hummingbird" 
focal <- "RUHU"
focal.dat.site <- presence.site %>% filter(bird.code == focal) %>% 
  #left_join(., covs, join_by("station.code" == "station.code")) %>%
  mutate(station.code = factor(station.code),
         yr = as.numeric(scale(year)))

view(focal.dat.site)

mod.null <- glmmTMB(cbind(detected, not_detected) ~ 1,
                family = binomial(link = "logit"), 
                data = focal.dat.site)

mod1 <- glmmTMB(cbind(detected, not_detected) ~ yr,
                family = binomial(link = "logit"), 
                data = focal.dat.site)

mod2 <- glmmTMB(cbind(detected, not_detected) ~ yr + park,
                family = binomial(link = "logit"), 
                data = focal.dat.site)

mod3 <- glmmTMB(cbind(detected, not_detected) ~ yr + park + (1 | station.code),
                family = binomial(link = "logit"),
                data = focal.dat.site)

mod4 <- glmmTMB(cbind(detected, not_detected) ~ yr + park + (1 | station.code) + (1| yr),
        family = binomial(link = "logit"),
        data = focal.dat.site)

mod5 <- glmmTMB(cbind(detected, not_detected) ~ yr + park + (yr | station.code),
        family = binomial(link = "logit"),
        data = focal.dat.site)

mod6 <- glmmTMB(cbind(detected, not_detected) ~ yr + (yr | station.code),
                family = binomial(link = "logit"),
                data = focal.dat.site)


AIC(mod.null, mod1, mod2, mod3, mod4, mod5, mod6)  ## most support for model 5. check diagnostics


mod5.pr <- glmmTMB(cbind(detected, not_detected) ~ yr + park + (yr | station.code),
                family = binomial(link = "probit"),
                data = focal.dat.site)


AIC(mod5, mod5.pr) # probit link does not improve fit

summary(mod5)

dispersion_test <- simulateResiduals(mod5)
plot(dispersion_test)

testZeroInflation(dispersion_test)
testOutliers(dispersion_test, type = "bootstrap")

## Overdispersion not an issue; zero inflation not an issue, outliers not an issue


## OK this model works; diagnostics are OK; overdispersion not an issue; zero inflation not an issue
intercept <- summary(mod5)$coefficients$cond[1, 1]
sd_year_effect <- summary(mod5)$coefficients$cond[2, 1]
sd_year_error <- summary(mod5)$coefficients$cond[2, 2]

mean_year <- mean(focal.dat.site$year)
sd_year <- sd(focal.dat.site$year)

# Function to convert log-odds to probability
source("functions/custom_functions.R")

yr_2005 <- (2005 - mean_year) / sd_year
yr_2023 <- (2023 - mean_year) / sd_year

logit_2005 <- intercept + (sd_year_effect * yr_2005)
logit_2023 <- intercept + (sd_year_effect * yr_2023)


logit_2005_lower <- intercept + (sd_year_effect - 1.96 * sd_year_error) * yr_2005
logit_2005_upper <- intercept + (sd_year_effect + 1.96 * sd_year_error) * yr_2005

logit_2023_lower <- intercept + (sd_year_effect - 1.96 * sd_year_error) * yr_2023
logit_2023_upper <- intercept + (sd_year_effect + 1.96 * sd_year_error) * yr_2023

p_2005 <- plogis(logit_2005)
p_2023 <- plogis(logit_2023)

p_2005_lower <- plogis(logit_2005_lower)
p_2005_upper <- plogis(logit_2005_upper)

p_2023_lower <- plogis(logit_2023_lower)
p_2023_upper <- plogis(logit_2023_upper)

(p_2023 - p_2005) / p_2005

(p_2023_lower - p_2005_upper) / p_2005_upper
(p_2023_upper - p_2005_lower) / p_2005_lower

# Probabilities for different years
probability_yr0 <- plogis(intercept + sd_year_effect * 0)  # yr = 0
probability_yr1 <- plogis(intercept + (sd_year_effect / sd_year) * 1)  # yr = 1
probability_yr2 <- plogis(intercept + (sd_year_effect / sd_year) * 2)  # yr = 2
probability_yr3 <- plogis(intercept + (sd_year_effect / sd_year) * 3)

probability_yr17 <- plogis(intercept + (sd_year_effect / sd_year) * 17)
probability_yr18 <- plogis(intercept + (sd_year_effect / sd_year) * 18)


# Results
probability_yr0
probability_yr1
probability_yr2


lower

(probability_yr3 - probability_yr2) / probability_yr2
(probability_yr2 - probability_yr1) / probability_yr1
(probability_yr1 - probability_yr0) / probability_yr0

(probability_yr18 - probability_yr0) / probability_yr0

# Create a new data frame for prediction
predict_data <- expand.grid(
  yr = seq(min(focal.dat.site$yr), max(focal.dat.site$yr), length.out = 100),  # time range for prediction
  station.code = levels(focal.dat.site$station.code)) %>%
  left_join(., station.table, join_by("station.code" == "station.code"))

# Create predictions (including fixed and random effects)
predicted <- predict(mod5, newdata = predict_data, type = "link", se.fit = TRUE, re.form = NULL)

predicted_df <- cbind(predict_data, predicted[[1]], predicted[[2]]) %>%
  rename(mu_logit = `predicted[[1]]`, se_logit = `predicted[[2]]`) %>%
  mutate(
    lower_logit = mu_logit - (1.96 * se_logit),
    upper_logit = mu_logit + (1.96 * se_logit),
    
    # Back-transforming from log-odds to probability
    mu = plogis(mu_logit),
    lower = plogis(lower_logit),
    upper = plogis(upper_logit)
  )


lower_logit <- predicted_df$mu_logit - 1.96 * predicted_df$se_logit
upper_logit<- predicted_df$mu_logit + 1.96 * predicted_df$se_logit

prob_lower <- plogis(lower_logit)
prob_upper <- plogis(upper_logit)


mean_predicted <- predicted_df %>%
  group_by(yr) %>%
  dplyr::summarise(mean_mu = mean(mu),
          mean_lower = mean(lower),
          mean_upper = mean(upper), .groups = "drop") %>%
  mutate(year = yr * sd_year + mean_year)

ggplot(mean_predicted, aes(x = yr, y = mean_mu)) +
  geom_line()


focal_mean <- focal.dat.site %>%
  mutate(prop = detected / nsurv) %>%
  group_by(yr) %>%
  dplyr::summarise(mean_prop = mean(prop), .groups = "drop") %>%
  mutate(year = yr * sd_year + mean_year)


ggplot(focal_mean, aes(x = year, y = mean_prop)) +
  geom_point() + 
  geom_line(data = mean_predicted, aes(x = year, y = mean_mu)) +
  #geom_ribbon(data = mean_predicted, aes(x = year, ymin = mean_lower, ymax = mean_upper), 
  #           alpha = 0.2, fill = bcs_colors["dark green"], inherit.aes = FALSE) + 
  labs(title = paste("Proportion of surveys reporting", focal.species.name, "over time"), x = "Year", y = "Proportion of surveys reporting") +
  theme_bcs()




## check for autocorrelation
n <- dim(focal.dat)[1]
res1<-residuals(simulateResiduals(mod5,n=nsim,seed=F))
plot(res1[1:(n-1)],res1[2:n],xlab="residual at year t",ylab="residual at year t+1")


y <- focal.dat.site$detected
n <- focal.dat.site$nsurv
npark <- length(unique(focal.dat.site$park))
yr <- focal.dat.site$yr
N <- dim(focal.dat.site)[1]
parks <- data.frame(park = sort(unique(focal.dat.site$park))) %>%
  mutate(park.code = row_number())
focal.dat.site <- left_join(focal.dat.site, parks)
park <- focal.dat.site$park.code
year_z <- data.frame(yr = sort(unique(focal.dat.site$yr))) %>%
  mutate(year.code = row_number())
focal.dat.site <- left_join(focal.dat.site, year_z)
year.code <- focal.dat.site$year.code




nchains <- 2
niters <- 10000
nburn <- niters / 2

data = list(y = y, n = n, N = N, yr = yr, park = park, npark = npark, year_z = year_z, year.code = year.code)
mod1 <- jags.model(textConnection(mod_string), data = data, n.chains = nchains)
monitor <-c("beta_0", "beta_year", "beta_park", "beta_park_year", "tau_park", "tau_year")
samp <- coda.samples(mod1, monitor, n.iter = niters)

colnames(samp[[1]])
ind <- grep("beta_year|tau_park|tau_year", colnames(samp[[1]]))
gelman.diag(samp[, ind], multivariate = FALSE)

dim(samp[[1]])

samp <- samp[(nburn + 1):niters, ]
res1 <- as.matrix(samp[[1]])
res2 <- as.matrix(samp[[2]])

par(mfrow = c(1,1))
for(i in ind){
  ylims <- range(res1[, i], res2[, i])
  plot(res1[, i], type = "l", ylim = ylims, xlab = "", ylab = "", main = colnames(res1)[i])
  lines(res2[, i], col = "red")
}

par(mfrow = c(2, 2))
for(i in ind){
  den1 <- density(res1[, i])
  den2 <- density(res2[, i])
  ylims <- range(den1$y, den2$y)
  xlims <- range(den1$x, den2$x)
  plot(den1$x, den1$y, type = "l", ylim = ylims, xlim = xlims, xlab = "", ylab = "", main = colnames(res1)[i])
  lines(den2$x, den2$y, col = "red")
}

res <- rbind(res1, res2)

par(mfrow = c(2, 2))
for(i in ind){
  den <- density(res[, i])
  plot(den$x, den$y, type = "l", xlab = "", ylab = "", main = colnames(res1)[i])
}


pairs(res[1:10^3,ind])
corr<-cor(res[,ind])
corr<-corr-diag(diag(corr))
round(range(corr,na.rm=T),2)

summarise(res[, ind])
rm(mod1)
parks