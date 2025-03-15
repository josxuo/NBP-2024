## Exploring Bayesian models for NBP data ##

# libraries
library(rjags)
load.module("glm")
library(tidyverse)
library(readxl)

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

#sort(unique(paste(d$park, d$loop, d$station)))

# Exclude Magnuson Park Back Fence Loop, Main Drag Loop, South End Loop due to 
# errors in data entry that make bird counts at any given circle unreliable.
# filter our sp records.
d <- nbp %>%
  filter(!str_detect(species, pattern = " sp\\."),
         station.code %in% covs$Station, 
         year %in% c(2012:2019, 2022, 2023),
         exclusions == "none")

sus <- unique(d$station.code)
sus.samp <- sus[sample(length(sus), 40, FALSE)]

d <- d %>%
  filter(station.code %in% sus.samp)

# inspect data
## sort(unique(d$species))
## sort(unique(d$year))
## sort(unique(d$code))


# Set up tables/calculations useful later
## number of surveys conducted each year
nsurvYearsp <- d %>%
  group_by(year, station.code) %>%
  reframe(nsurv = n_distinct(survey_id))
nsurvYearsp



## species response matrix -- mean abundance per survey by year
srmAbund <- d %>%
  # filter(station.code %in% sus, year %in% year, month %in% months) %>%
  group_by(year, species, park, loop, station.code) %>%
  reframe(count = sum(seen, heard, fly)) %>%
  left_join(., nsurvYearsp) %>%
  mutate(maps = round(count / nsurv, 0)) %>%
  pivot_wider(names_from = species, values_from = maps, values_fill = 0) %>%
  pivot_longer(-grep("year|park|loop|station|count|nsurv", colnames(.)), names_to = "species", values_to = "maps")

focal <- "Anna's Hummingbird"
focal.dat <- srmAbund %>% filter(species == focal) %>% mutate(yr = (year - min(year))) %>% 
  select(park, loop, station.code, maps, year, yr)


y <- focal.dat$maps
ns <- length(unique(focal.dat$station.code))
yr <- focal.dat$yr
n <- dim(focal.dat)[1]
np <- n_distinct(focal.dat$park)
sites <- data.frame(station.code = sort(unique(focal.dat$station.code))) %>%
  mutate(num.code = row_number())
parks <- data.frame(park = sort(unique(focal.dat$park))) %>%
  mutate(park.code = row_number())
focal.dat <- left_join(focal.dat, parks)
park <- focal.dat$park.code
focal.dat <- left_join(focal.dat, sites)
site <- focal.dat$num.code

nchains <- 2
niters <- 10000
nburn <- niters / 2

data = list(y = y, n = n, ns = ns, site = site, yr = yr, np = np, park = park)
mod1 <- jags.model("code/models/trend_w_rand_site_effect_negbin.R", data = data, n.chains = nchains)
monitor <- c("tr", "sigmas")
samp <- coda.samples(mod1, monitor, n.iter = niters)

colnames(samp[[1]])
ind <- grep("tr|sigmas", colnames(samp[[1]]))
gelman.diag(samp[, ind], multivariate = FALSE)


samp <- samp[(nburn + 1):niters, ]
res1 <- as.matrix(samp[[1]])
res2 <- as.matrix(samp[[2]])

par(mfrow = c(5,2))
for(i in ind){
  ylims <- range(res1[, i], res2[, i])
  plot(res1[, i], type = "l", ylim = ylims, xlab = "", ylab = "", main = colnames(res1)[i])
  lines(res2[, i], col = "red")
}

par(mfrow = c(5, 2))
for(i in ind){
  den1 <- density(res1[, i])
  den2 <- density(res2[, i])
  ylims <- range(den1$y, den2$y)
  xlims <- range(den1$x, den2$x)
  plot(den1$x, den1$y, type = "l", ylim = ylims, xlim = xlims, xlab = "", ylab = "", main = colnames(res1)[i])
  lines(den2$x, den2$y, col = "red")
}

res <- rbind(res1, res2)

par(mfrow = c(5, 2))
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
# interpreting the regression coefficients
# probability of detection has to be constant over time and space
# any changes we observe could be due to changes in observers, etc. 

tr0 <- res[ , grep("tr", colnames(res), fixed=T)]

mean.trend <- exp(tr0 / sd(focal.dat$year))
summarise(mean.trend)

par(mfrow = c(1, 1))
plot(focal.dat$year, focal.dat$maps)

colnames(res)
b0 <- res[, grep("b0", colnames(res))]
tr <- res[, grep("tr", colnames(res))] 


yrs <- unique(focal.dat$yr)

pop.lambda <- array(data = NA, dim = c(niters, length(yrs)))

for(i in 1:length(yrs)){
  pop.lambda[, i] <- exp(b0 + tr * yrs[i])
}

mean.pop.lambda <- apply(pop.lambda, 2, mean)
low.pop.lambda <- apply(pop.lambda, 2, quantile, probs = 0.025)
upp.pop.lambda <- apply(pop.lambda, 2, quantile, probs = 0.975)

par(mfrow=c(1,1))

ylims<-range(low.pop.lambda, upp.pop.lambda)
plot(yrs, mean.pop.lambda, ylab = "population index", ylim=ylims, xlab="year",
     pch=19, col="blue")
segments(yrs, low.pop.lambda, yrs, upp.pop.lambda, col="blue")




## But now I need to check other indicators of model fit--qqplot, for example.
## Posterior predictive check using overall discrepancy, posterior predictive check
## on observations, autocorrelation, distribution of error terms.

# Posterior predictive check
ysim <- res[, grep("ysim", colnames(res))]

hist(ysim)
hist(y)

T_obs_mean <- mean(y)
T_obs_var <- var(y)

T_sim_mean <- apply(ysim, 2, mean)
T_sim_var <- apply(ysim, 2, var)

ppp_mean <- mean(T_sim_mean >= T_obs_mean)
ppp_var <- mean(T_sim_var >= T_obs_var)

hist(T_sim_mean)
abline(v = T_obs_mean, col = "red")

lambda <- res[, grep("lambda", colnames(res))]
head(colnames(lambda))
tail(colnames(lambda))

mulambda <- apply(lambda, 1, mean)

summary(mulambda)

ysim <- res[, grep("ysim", colnames(res))]
head(colnames(ysim))
tail(colnames(ysim))

nsim <- 2000
ind <- sample(1:niters, nsim)
Dobs <- Dsim <- vector()

for(i in 1:n){
  y_hat_i <- mean(lambda[, i])
  Dobs[i] <- mean((y[i] - y_hat_i)^2 / y_hat_i)
  Dsim[i] <- mean((ysim[ind[i], ] - y_hat_i)^2 / y_hat_i)
}

mean((ysim[ind[i], i] - y_hat_i)^2 / y_hat_i)

summary(ysim[, 1859])

for(i in 1:nsim){
  Dobs[i] <- mean((y - mu[ind[i]])^2 / mu[ind[i]])  # pearson chi squared statistic for poisson data
  Dsim[i] <- mean((ysim[ind[i], ] - mu[ind[i]])^2 / mu[ind[i]])
}

summary(Dobs)
summary(Dsim)
su

par(mfrow = c(1, 2))
xylims <- range(Dsim, Dobs)

summary(xylims)
plot(Dsim, Dobs, xlim = xylims, ylim = xylims)
abline(0, 1)
hist(Dobs / Dsim, main = "")
mean(Dobs / Dsim > 1)  ## Bayesian p-value

# FROM AI

# Number of observations and posterior samples
N <- length(y)
S <- nrow(lambda)

# Initialize vectors to store chi-squared discrepancies
chi_sq_obs <- numeric(N)     # Chi-squared for observed data
chi_sq_rep <- matrix(NA, nrow = S, ncol = N)  # Chi-squared for replicated data

# Compute chi-squared discrepancy for each observation
for (i in 1:N) {
  # Calculate mean of posterior predictive for each observation (expected count)
  y_hat_i <- mean(lambda[, i])
  
  # Chi-squared for observed data
  chi_sq_obs[i] <- (y[i] - y_hat_i)^2 / y_hat_i
  
  # Chi-squared for replicated data (across posterior samples)
  for (s in 1:S) {
    chi_sq_rep[s, i] <- (ysim[s, i] - lambda[s, i])^2 / lambda[s, i]
  }
}

T_obs_chi_sq <- sum(chi_sq_obs)
T_sim_chi_sq <- apply(chi_sq_rep, 2, mean)

ppp_chi_sq <- mean(T_sim_chi_sq >= chi_sq_obs)

## predictive check showing possible overdisperion. Try with zero-inflated model?

focal <- "Dark-eyed Junco"
focal.dat <- srmAbund %>% filter(species == focal) %>% mutate(yr = as.numeric(scale(year))) %>% 
  select(park, loop, station.code, maps, year, yr)


y <- focal.dat$maps
ns <- length(unique(focal.dat$station.code))
yr <- focal.dat$yr
n <- dim(focal.dat)[1]
sites <- data.frame(station.code = sort(unique(focal.dat$station.code))) %>%
  mutate(num.code = row_number())
focal.dat <- left_join(focal.dat, sites)
site <- focal.dat$num.code

nchains <- 2
niters <- 10000
nburn <- niters / 2

data = list(y = y, n = n, ns = ns, site = site, yr = yr)
mod1 <- jags.model("code/models/trend_w_rand_site_effect_zi-negbin.R", data = data, n.chains = nchains)
monitor <- c("b0", "tr", "sigmas", "ysim", "r", "y", "lambda", "psi")
samp <- coda.samples(mod1, monitor, n.iter = niters)

zero_inflated[i] * dnegbin(p[i], r) + (1 - zero_inflated[i]) * 0

(1*.4) + (1-1) * 0
