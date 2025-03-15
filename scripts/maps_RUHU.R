## MODELING BARN SWALLOW MEAN ABUNDANCE PER SURVEY ##

# clear environment

rm(list = ls())

# libraries
library(tidyverse)
library(readxl)
library(tidyverse)
library(lme4)

# load data

nbp <- read_xlsx("data/b_intermediate_data/nbp_tidy_jan_24.xlsx")
covs <- read_csv("data/c_analysis_data/circ_no_overlap_covariates.csv")


ilogit <- function(x){
  exp(x) / (1 + exp(x))
}

qqunif<-function(res,nsim=10^3,lab="",xlims=c(0,1),ylims=c(0,1)){
  res<-sort(res)
  n<-length(res)
  usim<-array(NA,c(nsim,n))
  for(i in 1:nsim){ usim[i,]<-sort(runif(n)) }
  u.med<-apply(usim,2,median)
  u.low<-apply(usim,2,quantile,probs=0.025)
  u.upp<-apply(usim,2,quantile,probs=0.975)
  plot(u.med,res,ylab="observed",xlab="expected",xlim=xlims,ylim=ylims,main=lab)
  lines(u.med,u.med)
  lines(u.med,u.low)
  lines(u.med,u.upp)
  out<-which((res>u.upp)|(res<u.low))
  p.out<-length(out)/n
  if(p.out>0){
    message(paste("observations outside envelope (",round(100*p.out,1),"%)",sep=""))
    out
  }
  else{ message("no observations outside envelope") }
}

phi<-function(m){ sum(residuals(m,type="pearson")^2)/df.residual(m) }


d <- nbp %>%
  # filter for years after 2003, excluding incomplete years 2020 and 2024.
  # filter out sp. records
  filter(!str_detect(species, pattern = " sp\\."),
         # filter out introduced species
         # !species %in% intros,
         # filter for just non-overlapping count circles
         station.code %in% covs$Station, 
         exclusions == "none")

# inspect data
## sort(unique(d$species))
## sort(unique(d$year))
## sort(unique(d$code))


# Set up tables/calculations useful later
## number of surveys conducted each year
nsurvYearsp <- d %>%
  #  filter(year %in% years, month %in% months, station.code %in% sus) %>%
  group_by(year, species, station.code) %>%
  summarise(nsurv = n_distinct(survey_id))
nsurvYearsp



## species response matrix -- mean abundance per survey by year
srmAbund <- d %>%
  # filter(station.code %in% sus, year %in% year, month %in% months) %>%
  group_by(year, species, park, loop, station.code) %>%
  summarise(count = sum(seen, heard, fly)) %>%
  ungroup() %>%
  left_join(., nsurvYearsp) %>%
  mutate(maps = round(count / nsurv, 0)) %>%
  pivot_wider(names_from = species, values_from = maps) %>%
  replace(is.na(.), 0) %>%
  pivot_longer(-grep("year|park|loop|station|count|nsurv", colnames(.)), names_to = "species", values_to = "maps")


focal <- "Steller's Jay"
focal.dat <- srmAbund %>% filter(species == focal) %>% mutate(yr = as.numeric(scale(year))) %>% 
  select(park, loop, station.code, maps, year, yr)

yr <- focal.dat$yr
n.yr <- length(unique(yr))

plot(focal.dat$yr, focal.dat$maps)

m1 <- glmmTMB(maps ~ yr + (1 | station.code), data = focal.dat, family = poisson)
m2 <- glmmTMB(maps ~ yr + (1 | station.code), data = focal.dat, family = nbinom1)
m3 <- glmmTMB(maps ~ yr + (1 | station.code), data = focal.dat, family = nbinom2)

AICtab(m1, m2, m3, mnames = c("m1", "m2", "m3"))

## most support for nbiom1. Compare agains null models
m01 <- glmmTMB(maps ~ 1, data = focal.dat, family = nbinom1)
m02 <- glmmTMB(maps ~ 1 + (1 | station.code), data = focal.dat, family = nbinom1)

AICtab(m2, m01, m02, mnames = c("m2", "m01", "m02")) ## M2 still supported

summary(m2)

mod <- m2

### qq-plot for pearson residuals

par(mfrow=c(1,1))
res<-residuals(mod,type="pearson")
qqnorm(res,xlab="expected",ylab="observed",main="")
abline(0,1,lwd=1.5)

### diagnostic plots using PIT residuals
nsim <- 1000
n <- dim(focal.dat)[1]

res<-residuals(simulateResiduals(mod,n=nsim,seed=F))
qqunif(res,nsim)
plot(rank(fitted(mod)) / n, res, xlab = "fitted values (rank-transformed)", ylab = "residuals")
boxplot(res ~ yr, data= focal.dat, ylab= "residuals", xlab="")
stripchart(res~focal.dat$yr,vertical=T,xlab="",ylab="residuals",main="",pch=21,xlim=c(0,n.yr)+0.5)

par(mfrow=c(1,2))

### check for overdispersion

phi.obs<-phi(mod)
phi.sim<-vector()
for(i in 1:nsim){
  print(i)
  ysim <- simulate(mod)$sim_1
  phi.sim[i] <- phi(glmmTMB(ysim ~ yr + (1 | station.code), family = nbinom1 ,data=focal.dat))
}
phi.p<-mean(phi.sim>=phi.obs)
hmax<-max(hist(phi.sim,plot=F)$density)
hist(phi.sim,main=paste("Dispersion (p=",round(phi.p,2),")",sep=""),freq=F,xlab="",yaxt="n",ylab="",xlim=range(phi.sim,phi.obs))
segments(phi.obs,0,phi.obs,hmax,col="red",lwd=2)

### Overdispersion does not appear to be an issue here.

### check for zero-inflation

pzero.obs<-mean(focal.dat$maps==0)
pzero.sim<-vector()
for(i in 1:nsim){ pzero.sim[i]<-mean(simulate(mod)$sim_1==0) }
pzero.p<-mean(pzero.sim>=pzero.obs)
hmax<-max(hist(pzero.sim,plot=F)$density)
hist(pzero.sim,main=paste("Proportion of zeros (p=",round(pzero.p,2),")",sep=""),freq=F,xlab="",yaxt="n",ylab="",xlim=range(pzero.sim,pzero.obs))
segments(pzero.obs,0,pzero.obs,hmax,col="red",lwd=2)

## Zero inflation also ok

summary(mod)

## plot results

yrs <- unique(focal.dat$yr)
stations <- unique(focal.dat$station.code)

newdata <- data.frame(yr = rep(yrs, length(stations)),
                      station.code = rep(stations, length(yrs)))

fit <- predict(mod, newdata = newdata, type = "response", se.fit = TRUE)

str(fit)

head(fit)

preds <- cbind(newdata, fit$fit, fit$se.fit)

mean.preds <- preds %>% group_by(yr) %>% reframe(mean.fit = mean(`fit$fit`), mean.se = mean(`fit$se.fit`))
upp <- mean.preds$mean.fit + 1.96*mean.preds$mean.se
low <- mean.preds$mean.fit - 1.96*mean.preds$mean.se

mean.maps <- focal.dat %>% group_by(year) %>% reframe(mean = mean(maps))

par(mfrow = c(1, 1))

ylims <- range(upp, low)

plot(x = c(1996:2024), y = mean.preds$mean.fit, ylim = ylims, type = "l")
points(mean.maps$year, mean.maps$mean)
lines(x = c(1996:2024), y = upp, lty = 2)
lines(x = c(1996:2024), y = low, lty = 2)

effect <- (exp(summary(mod)$coefficients[1]$cond[2, 1]/sd(focal.dat$year)) - 1) * 100
summary(mod)

(exp(.0477 / sd(focal.dat$year)) - 1) * 100

ggplot(mean.maps, aes(x = year, y = mean)) +
  geom_point(color = bcs_colors["dark green"]) +
  geom_line(aes(x = c(1996:2024), y = mean.preds$mean.fit), color = bcs_colors["bright green"]) +
  geom_line(aes(x = c(1996:2024), y = upp), color = bcs_colors["bright green"], linetype = 2) +
  geom_line(aes(x = c(1996:2024), y = low), color = bcs_colors["bright green"], linetype = 2) +
  labs(x = "Year", y = "Mean abundance per survey", title = "Mean abundance per survey", subtitle = focal) +
  annotate("text", x =2002, y = 0.31, label = paste0("Mean abundance per survey decreasing by ", 
                                                     round(effect, 1), "% per year"),
           size = 5, color = bcs_colors["dark green"]) +
  theme_bcs()

