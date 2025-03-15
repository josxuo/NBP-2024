## Exploring baysian models for NBP data ##

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

complete <- function(years, months, data) {
  int <- data %>%
    filter(year %in% years, month %in% months, exclusions == "none")
  
  comp <- int %>%
    group_by(station.code, year) %>%
    reframe(nsurv = n_distinct(survey_id)) %>%
    pivot_wider(names_from = year, values_from = nsurv) %>%
    replace(is.na(.), 0) %>%
    mutate(row_sum = apply(.[, -1], 1, sum)) %>%
    filter(row_sum == length(years) * length(months)) %>%
    pull(station.code)
}

# data
d <- read_xlsx("data/b_intermediate_data/nbp_tidy_jan_24.xlsx")
sort(unique(paste(d$park, d$loop, d$station)))

# Exclude Magnuson Park Back Fence Loop, Main Drag Loop, South End Loop due to 
# errors in data entry that make bird counts at any given circle unreliable
nbp <- d %>% filter(exclusions == "none")


## Single month trend analysis

# subset data for circles with consistent data
years <-  c(2006:2019, 2022)
months <- 2
sus <- complete(years, months, data = nbp)

counts <- nbp %>%
#  filter(station.code %in% sus, 
   #      year %in% years, month %in% months) %>%
  group_by(station.code, year, month, bird.code) %>%
  reframe(count = sum(seen, heard, fly)) %>%
  pivot_wider(names_from = bird.code, values_from = count) %>%
  replace(is.na(.), 0) %>%
  pivot_longer(-c(1, 2, 3), names_to = "bird.code", values_to = "count")

counts$bird.code

bird <- filter(counts, bird.code == "DEJU")

y <- bird$count
ns <- length(unique(bird$station.code))
yr <- as.numeric(scale(bird$year))
n <- dim(bird)[1]
sites <- data.frame(station.code = sort(unique(bird$station.code))) %>%
  mutate(num.code = row_number())
bird <- left_join(bird, sites)
site <- bird$num.code

nchains <- 2
niters <- 10000
nburn <- niters / 2

data = list(y = y, n = n, ns = ns, site = site, yr = yr)
mod1 <- jags.model("code/models/trend_w_rand_site_effect.R", data = data, n.chains = nchains)
monitor <- c("b0", "bs", "tr0", "sigmas", "sigmatr")
samp <- coda.samples(mod1, monitor, n.iter = niters)

colnames(samp[[1]])
ind <- grep("b0|tr0|sigma", colnames(samp[[1]]))
gelman.diag(samp[, ind], multivariate = FALSE)

samp <- samp[(nburn + 1):niters, ]
res1 <- as.matrix(samp[[1]])
res2 <- as.matrix(samp[[2]])

par(mfrow = c(5, 5))
for(i in c(1:49)){
  ylims <- range(res1[, i], res2[, i])
  plot(res1[, i], type = "l", ylim = ylims, xlab = "", ylab = "", main = colnames(res1)[i])
  lines(res2[, i], col = "red")
}

par(mfrow = c(5, 5))
for(i in 1:49){
  den1 <- density(res1[, i])
  den2 <- density(res2[, i])
  ylims <- range(den1$y, den2$y)
  xlims <- range(den1$x, den2$x)
  plot(den1$x, den1$y, type = "l", ylim = ylims, xlim = xlims, xlab = "", ylab = "", main = colnames(res1)[i])
  lines(den2$x, den2$y, col = "red")
}

res <- rbind(res1, res2)

par(mfrow = c(2, 3))
for(i in ind){
  den <- density(res[, i])
  plot(den$x, den$y, type = "l", xlab = "", ylab = "", main = colnames(res1)[i])
}


pairs(res[1:10^3,ind])
corr<-cor(res[,ind])
corr<-corr-diag(diag(corr))
round(range(corr,na.rm=T),2)

summarise(res[, ind])

## FOR DEJU, SD between sites on linear predictor scale is ~1.9
## FOR DEJU, SD between trends on linear predictor scale is ~0.5

# interpreting the regression coefficients
# probability of detection has to be constant over time and space
# any changes we observe could be due to changes in observers, etc. 

tr0 <- res[ , grep("tr0", colnames(res), fixed=T)]

mean.trend <- exp(tr0 / sd(bird$year))
summarise(mean.trend)

plot(bird$year, bird$count)

colnames(res)
b0<-res[, grep("b0", colnames(res))]
tr0 <- res[, grep("tr0", colnames(res))] 

yrmean <- mean(bird$year)
yrsd <- sd(bird$year)
yrs <- sort(unique(bird$year))
scyr <- (yrs - yrmean) / yrsd

pop.mu <- array(data = NA, dim = c(niters, length(scyr)))

for(i in 1:length(scyr)){
  pop.mu[, i] <- exp(b0 + tr0*scyr[i])
}

mean.pop.mu <- apply(pop.mu, 2, mean)
low.pop.mu <- apply(pop.mu, 2, quantile, probs = 0.025)
upp.pop.mu <- apply(pop.mu, 2, quantile, probs = 0.975)

par(mfrow=c(1,1))

ylims<-range(low.pop.mu,upp.pop.mu)
plot(years,mean.pop.mu,ylab="population index",ylim=ylims,xlab="year",pch=19,col="blue")
segments(years,low.pop.mu,years,upp.pop.mu,col="blue")


## TRY WITH ANOTHER SPECIES

bird <- filter(counts, bird.code == "RUHU")

y <- bird$count
ns <- length(unique(bird$station.code))
yr <- as.numeric(scale(bird$year))
n <- dim(bird)[1]
sites <- data.frame(station.code = sort(unique(bird$station.code))) %>%
  mutate(num.code = row_number())
bird <- left_join(bird, sites)
site <- bird$num.code

nchains <- 2
niters <- 10000
nburn <- niters / 2

data = list(y = y, n = n, ns = ns, site = site, yr = yr)
mod1 <- jags.model("code/models/trend_w_rand_site_effect.R", data = data, n.chains = nchains)
monitor <- c("b0", "bs", "tr0", "sigmas", "sigmatr")
samp <- coda.samples(mod1, monitor, n.iter = niters)

colnames(samp[[1]])
ind <- grep("b0|tr0|sigma", colnames(samp[[1]]))
gelman.diag(samp[, ind], multivariate = FALSE)

samp <- samp[(nburn + 1):niters, ]
res1 <- as.matrix(samp[[1]])
res2 <- as.matrix(samp[[2]])

par(mfrow = c(5, 5))
for(i in c(1:49)){
  ylims <- range(res1[, i], res2[, i])
  plot(res1[, i], type = "l", ylim = ylims, xlab = "", ylab = "", main = colnames(res1)[i])
  lines(res2[, i], col = "red")
}

par(mfrow = c(5, 5))
for(i in 1:49){
  den1 <- density(res1[, i])
  den2 <- density(res2[, i])
  ylims <- range(den1$y, den2$y)
  xlims <- range(den1$x, den2$x)
  plot(den1$x, den1$y, type = "l", ylim = ylims, xlim = xlims, xlab = "", ylab = "", main = colnames(res1)[i])
  lines(den2$x, den2$y, col = "red")
}

res <- rbind(res1, res2)

par(mfrow = c(2, 3))
for(i in ind){
  den <- density(res[, i])
  plot(den$x, den$y, type = "l", xlab = "", ylab = "", main = colnames(res1)[i])
}


pairs(res[1:10^3,ind])
corr<-cor(res[,ind])
corr<-corr-diag(diag(corr))
round(range(corr,na.rm=T),2)

summarise(res[, ind])

## FOR RCKI, SD between sites on linear predictor scale is ~2
## FOR RCKI, SD between trends on linear predictor scale is ~-0.3

# interpreting the regression coefficients
# probability of detection has to be constant over time and space
# any changes we observe could be due to changes in observers, etc. 

tr0 <- res[ , grep("tr0", colnames(res), fixed=T)]

mean.trend <- exp(tr0 / sd(bird$year))
summarise(mean.trend)

plot(bird$year, bird$count)

colnames(res)
b0<-res[, grep("b0", colnames(res))]
tr0 <- res[, grep("tr0", colnames(res))] 

yrmean <- mean(bird$year)
yrsd <- sd(bird$year)
yrs <- sort(unique(bird$year))
scyr <- (yrs - yrmean) / yrsd

pop.mu <- array(data = NA, dim = c(niters, length(scyr)))

for(i in 1:length(scyr)){
  pop.mu[, i] <- exp(b0 + tr0*scyr[i])
}

mean.pop.mu <- apply(pop.mu, 2, mean)
low.pop.mu <- apply(pop.mu, 2, quantile, probs = 0.025)
upp.pop.mu <- apply(pop.mu, 2, quantile, probs = 0.975)

par(mfrow=c(1,1))

ylims<-range(low.pop.mu,upp.pop.mu)
plot(years,mean.pop.mu,ylab="population index",ylim=ylims,xlab="year",pch=19,col="blue")
segments(years,low.pop.mu,years,upp.pop.mu,col="blue")


focal <- c("AMRO", "STJA", "AMCR", "NOFL", "GBHE", "ANHU")

yrmean <- mean(bird$year)
yrsd <- sd(bird$year)
yrs <- sort(unique(bird$year))
scyr <- (yrs - yrmean) / yrsd

nchains <- 2
niters <- 30000
nburn <- niters / 2
x <- list(array(NA, dim = c(3, length(scyr))))
pop.mus <- rep(x, times = length(focal))

monitor <- c("b0", "bs", "tr0", "sigmas", "sigmatr")

for(i in 1:length(focal)){
  
  bird <- filter(counts, bird.code == focal[i])
  
  y <- bird$count
  ns <- length(unique(bird$station.code))
  yr <- as.numeric(scale(bird$year))
  n <- dim(bird)[1]
  sites <- data.frame(station.code = sort(unique(bird$station.code))) %>%
    mutate(num.code = row_number())
  bird <- left_join(bird, sites)
  site <- bird$num.code
  
  data = list(y = y, n = n, ns = ns, site = site, yr = yr)

  mod1 <- jags.model("code/models/trend_w_rand_site_effect.R", data = data, n.chains = nchains)
  samp <- coda.samples(mod1, monitor, n.iter = niters)
  
  ind <- grep("b0|tr0|sigma", colnames(samp[[1]]))
  
  print(focal[i])
  print(gelman.diag(samp[, ind], multivariate = FALSE))
  
  samp <- samp[(nburn + 1):niters, ]
  res1 <- as.matrix(samp[[1]])
  res2 <- as.matrix(samp[[2]])
  
  par(mfrow = c(2, 3))
  for(j in ind){
    ylims <- range(res1[, j], res2[, j])
    plot(res1[, j], type = "l", ylim = ylims, xlab = "", ylab = "", main = paste(focal[i], colnames(res1)[j]))
    lines(res2[, j], col = "red")
  }
  
  par(mfrow = c(2, 3))
  for(j in ind){
    den1 <- density(res1[, j])
    den2 <- density(res2[, j])
    ylims <- range(den1$y, den2$y)
    xlims <- range(den1$x, den2$x)
    plot(den1$x, den1$y, type = "l", ylim = ylims, xlim = xlims, xlab = "", ylab = "", main = paste(focal[i], colnames(res1)[j]))
    lines(den2$x, den2$y, col = "red")
  }
  
  res <- rbind(res1, res2)
  
  par(mfrow = c(2, 3))
  for(j in ind){
    den <- density(res[, j])
    plot(den$x, den$y, type = "l", xlab = "", ylab = "", main = paste(focal[i], colnames(res1)[j]))
  }
  
  
  pairs(res[1:10^3,ind])
  corr<-cor(res[,ind])
  corr<-corr-diag(diag(corr))
  round(range(corr,na.rm=T),2)
  
  print(focal[i])
  print(summarise(res[, ind]))
  
  b0<-res[, grep("b0", colnames(res))]
  tr0 <- res[, grep("tr0", colnames(res))] 

  pop.mu <- array(data = NA, dim = c(niters, length(scyr)))
  
  for(j in 1:length(scyr)){
    pop.mu[, j] <- exp(b0 + tr0*scyr[j])
  }
  
  pop.mus[[i]][1, ] <- apply(pop.mu, 2, mean)
  pop.mus[[i]][2, ] <- apply(pop.mu, 2, quantile, probs = 0.025)
  pop.mus[[i]][3, ] <- apply(pop.mu, 2, quantile, probs = 0.975)
  
  par(mfrow=c(1,1))
  
  ylims<-range(pop.mus[[i]][3, ], pop.mus[[i]][2, ])
  plot(years,pop.mus[[i]][1, ], ylab = "population index", ylim = ylims, 
       xlab="year", pch=19, col="blue",
       main = focal[i])
  segments(years, pop.mus[[i]][2, ], years, pop.mus[[i]][3, ], col="blue")
}

pop.mus
mean.trend <- exp(-.1872/ sd(bird$year))
summarise(mean.trend)



#### TRY ADDING RANDOM MONTH EFFECT

bird <- filter(counts, year > 2015 & bird.code == "DEJU")

y <- bird$count
ns <- length(unique(bird$station.code))
yr <- as.numeric(scale(bird$year))
n <- dim(bird)[1]
ny <- n_distinct(bird$year)
sites <- data.frame(station.code = sort(unique(bird$station.code))) %>%
  mutate(num.code = row_number())
bird <- left_join(bird, sites)
site <- bird$num.code
mon <- bird$month
year <- bird$year - min(bird$year) + 1

nchains <- 2
niters <- 5000
nburn <- niters / 2

data = list(y = y, n = n, ns = ns, ny = ny, site = site, yr = yr, year = year, mon = mon)
mod1 <- jags.model("code/models/trend_w_rand_site_effect_and_error.R", data = data, n.chains = nchains)
monitor <- c("b0", "tr10", "bs", "mu", "sigmas", "sigmatr1", "ysim", "sigmae")
samp <- coda.samples(mod1, monitor, n.iter = niters)

colnames(samp[[1]])
ind <- grep("tr|sigma", colnames(samp[[1]]))
gelman.diag(samp[, ind], multivariate = FALSE)

samp <- samp[(nburn + 1):niters, ]
res1 <- as.matrix(samp[[1]])
res2 <- as.matrix(samp[[2]])

par(mfrow = c(2, 3))
for(i in ind){
  ylims <- range(res1[, i], res2[, i])
  plot(res1[, i], type = "l", ylim = ylims, xlab = "", ylab = "", main = colnames(res1)[i])
  lines(res2[, i], col = "red")
}

par(mfrow = c(2, 3))
for(i in ind){
  den1 <- density(res1[, i])
  den2 <- density(res2[, i])
  ylims <- range(den1$y, den2$y)
  xlims <- range(den1$x, den2$x)
  plot(den1$x, den1$y, type = "l", ylim = ylims, xlim = xlims, xlab = "", ylab = "", main = colnames(res1)[i])
  lines(den2$x, den2$y, col = "red")
}

res <- rbind(res1, res2)

par(mfrow = c(2,3))
for(i in ind){
  den <- density(res[, i])
  plot(den$x, den$y, type = "l", xlab = "", ylab = "", main = colnames(res1)[i])
}


pairs(res[1:10^3,ind])
corr<-cor(res[,ind])
corr<-corr-diag(diag(corr))
round(range(corr,na.rm=T),2)

summarise(res[, ind])

# interpreting the regression coefficients
# probability of detection has to be constant over time and space
# any changes we observe could be due to changes in observers, etc. 

tr0 <- res[ , grep("tr0", colnames(res), fixed=T)]

mean.trend <- exp(tr0 / sd(bird$year))
summarise(mean.trend)

plot(bird$year, bird$count)

colnames(res)
b0<-res[, grep("b0", colnames(res))]
tr0 <- res[, grep("tr0", colnames(res))] 

yrmean <- mean(bird$year)
yrsd <- sd(bird$year)
yrs <- sort(unique(bird$year))
scyr <- (yrs - yrmean) / yrsd

pop.mu <- array(data = NA, dim = c(niters, length(scyr)))

for(i in 1:length(scyr)){
  pop.mu[, i] <- exp(b0 + tr0*scyr[i])
}

mean.pop.mu <- apply(pop.mu, 2, mean)
low.pop.mu <- apply(pop.mu, 2, quantile, probs = 0.025)
upp.pop.mu <- apply(pop.mu, 2, quantile, probs = 0.975)

par(mfrow=c(1,1))

years <- sort(unique(bird$year))

ylims<-range(low.pop.mu,upp.pop.mu)
plot(years,mean.pop.mu,ylab="population index",ylim=ylims,xlab="year",pch=19,col="blue")
segments(years,low.pop.mu,years,upp.pop.mu,col="blue")
points(bird$year, bird$count)

sum(mean.trend > 1.)

## But now I need to check other indicators of model fit--qqplot, for example.
## Posterior predictive check using overall discrepancy, posterior predictive check
## on observations, autocorrelation, distribution of error terms.

# Posterior predictive check
ysim <- res[, grep("ysim", colnames(res))]
mu <- res[, grep("mu", colnames(res))]

nsim <- 1000
ind <- sample(1:niters, nsim)
Dobs <- Dsim <- vector()
for(i in 1:nsim){
  Dobs[i] <- mean((y - mu[ind[i], ])^2 / mu[ind[i]])  # pearson chi squared statistic for poisson data
  Dsim[i] <- mean((ysim[ind[i], ] - mu[ind[i], ])^2 / mu[ind[i], ])
}

par(mfrow = c(1, 2))
xylims <- range(Dsim, Dobs)
plot(Dsim, Dobs, xlim = xylims, ylim = xylims)
abline(0, 1)
hist(Dobs / Dsim, main = "")
mean(Dobs / Dsim > 1)  ## Bayesian p-value


## plot data and model
mean.mu <- apply(mu, 2, mean)
low.mu <- apply(mu, 2, quantile, probs = 0.025)
upp.mu <- apply(mu, 2, quantile, probs = 0.975)

mean.mu <- as.data.frame(cbind(bird$year, unname(mean.mu)))
mean.mu <- mean.mu %>% group_by(V1) %>% reframe(mean = mean(V2))

low.mu <- as.data.frame(cbind(bird$year, unname(low.mu))) %>% group_by(V1) %>% reframe(mean = mean(V2))
upp.mu <- as.data.frame(cbind(bird$year, unname(upp.mu))) %>% group_by(V1) %>% reframe(mean = mean(V2))


ylims <- range(bird$count, low.ysim, upp.ysim)
par(mfrow = c(1, 1))
plot(bird$year, bird$count)
lines(mean.mu$V1, mean.mu$mean, col = "blue")

lines(low.mu$V1, low.mu$mean, lty = 3)
lines(upp.mu$V1, upp.mu$mean, lty = 3)

bird.means <- bird %>% group_by(year) %>% reframe(mean.count = mean(count))

ylims <- range(bird.means$mean.count, low.mu$mean, upp.mu$mean)
plot(bird.means$year, bird.means$mean.count, ylim = ylims)
lines(mean.mu$V1, mean.mu$mean, col = "blue")

lines(low.mu$V1, low.mu$mean, lty = 3)
lines(upp.mu$V1, upp.mu$mean, lty = 3)

rpois(100, .5999)


d <- array(data = NA, dim = c(12, 16, 99))
d
