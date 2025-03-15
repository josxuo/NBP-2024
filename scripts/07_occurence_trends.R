#############################
####### 00. occupancy #######
#############################

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

### ANALYSIS OBJECTIVE: IDENTIFY SPECIES EXPERIENCING FASTEST LOCAL DECLINES ###
# Seek: 1) Declines in reported abundance; 2) Declines in proportion of surveys reporting;
# 3) Declines in occupancy.

## t1 = 2004
## t2 = 2023
## exclude 2020

## aggregate data across all parks

# functions
## extract overall p value from linear models
overall_p <- function(my_model){
  f <- summary(my_model)$fstatistic
  p <- pf(f[1], f[2], f[3], lower.tail = FALSE)
  attributes(p) <- NULL
  return(p)
}

# function for returning count circle codes with a complete set of survey data for a given set of years and months
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

# will want to exclude introduced species
intros <- c("European Starling", "House Sparrow", "Rock Pigeon", "Eurasian Collared-Dove",
            "California Quail")


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

## Select circles with complete data for a given set of years, months to analyze
# years <- c(2005:2019, 2022, 2023)
# months <- c(1:12)
# sus <- complete(years, months, data = d)

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
  



## Species mean abundance per survey trends
d.maps <- srmAbund %>% filter(year %in% years) %>%
  mutate(yr = factor(year - min(year) + 1, levels= c(1:19)),
         station.code = as.factor(station.code),
         park = as.factor(park), 
         loop = as.factor(loop))


focal <- "Dark-eyed Junco"
focal.dat <- d.maps %>% filter(species == focal) %>% select(park, loop, station.code, maps, year, yr)
yr <- focal.dat$yr
n.yr <- length(unique(yr))

plot(focal.dat$yr, focal.dat$maps)

test <- glm(maps ~ yr, data = focal.dat, family = poisson(link = "log"))
test <- glmmTMB(maps ~ yr + (1 | station.code), data = focal.dat, family = nbinom2)
summary(test)

(emm1<-emmeans(test,~yr,type="lp"))
(emm2<-emmeans(test,~yr,type="response"))

(confint(pairs(emm1,adjust="none")))
(confint(pairs(emm2,adjust="none")))

(confint(contrast(emm1,method="del.eff",adjust="none")))
(confint(contrast(emm2,method="del.eff",adjust="none")))

(contrast(emm1,method="del.eff",adjust="none"))
(contrast(emm2,method="del.eff",adjust="none"))

### plot estimates and confidence intervals

(df.emm1<-as.data.frame(emm1))
(df.emm2<-as.data.frame(emm2))

str(df.emm2)

par(mfrow=c(1,1))
ylims<-range(df.emm2$asymp.LCL,df.emm2$asymp.UCL)
stripchart(df.emm2$response~unique(yr),vertical=T,ylim=ylims,xlab="",ylab="",main="",pch=19,xlim=c(0,n.yr)+0.5)
segments(1:n.yr,df.emm2$asymp.LCL,1:n.yr,df.emm2$asymp.UCL)

### plot confidence densities (see reference on confidence densities)

par(mfrow=c(2,4))
for(i in 1:n.yr){
  r<-exp(rnorm(10^6,df.emm1$emmean[i],df.emm1$SE[i]))
  rlim<-quantile(r,probs=c(0.001,0.999))
  plot(density(r),main=unique(yr)[i],xlab="",ylab="",xlim=rlim)
}

### qq-plot for pearson residuals

par(mfrow=c(1,1))
res<-residuals(test,type="pearson")
qqnorm(res,xlab="expected",ylab="observed",main="")
abline(0,1,lwd=1.5)

### diagnostic plots using PIT residuals

res<-residuals(simulateResiduals(test,n=nsim,seed=F))
qqunif(res,nsim)
plot(rank(fitted(test))/n,res,xlab="fitted values (rank-transformed)",ylab="residuals")
boxplot(res~yr,data=focal.dat,ylab="residuals",xlab="")
stripchart(res~focal.dat$yr,vertical=T,xlab="",ylab="residuals",main="",pch=21,xlim=c(0,n.yr)+0.5)

par(mfrow=c(1,2))

### check for overdispersion

phi.obs<-phi(test)
phi.sim<-vector()
for(i in 1:nsim){
  print(i)
  ysim<-simulate(test)$sim_1
  phi.sim[i]<-phi(glmmTMB(ysim~ yr + (1 | station.code), family=nbinom2 ,data=focal.dat))
}
phi.p<-mean(phi.sim>=phi.obs)
hmax<-max(hist(phi.sim,plot=F)$density)
hist(phi.sim,main=paste("Dispersion (p=",round(phi.p,2),")",sep=""),freq=F,xlab="",yaxt="n",ylab="",xlim=range(phi.sim,phi.obs))
segments(phi.obs,0,phi.obs,hmax,col="red",lwd=2)

### Overdispersion does not appear to be an issue here.

### check for zero-inflation

pzero.obs<-mean(focal.dat$maps==0)
pzero.sim<-vector()
for(i in 1:nsim){ pzero.sim[i]<-mean(simulate(test)$sim_1==0) }
pzero.p<-mean(pzero.sim>=pzero.obs)
hmax<-max(hist(pzero.sim,plot=F)$density)
hist(pzero.sim,main=paste("Proportion of zeros (p=",round(pzero.p,2),")",sep=""),freq=F,xlab="",yaxt="n",ylab="",xlim=range(pzero.sim,pzero.obs))
segments(pzero.obs,0,pzero.obs,hmax,col="red",lwd=2)

## Zero inflation also ok

# Compare various model options

m1 <- glmmTMB(maps ~ yr + (1 | station.code), data = focal.dat, family = poisson)
m2 <- glmmTMB(maps ~ yr + (1 | station.code), data = focal.dat, family = nbinom1)
m3 <- glmmTMB(maps ~ yr + (1 | station.code), data = focal.dat, family = nbinom2)
m0 <- glmmTMB(maps ~ 1 + (1 | station.code), data = focal.dat, family = nbinom2)

AICtab(m1,m2,m3,m0,logLik=T,base=T,mnames=c("Poisson","NB1","NB2", "nullNB2"))

# Switch to nb1
mod <- m2

res<-residuals(simulateResiduals(mod,n=nsim,seed=F))
qqunif(res,nsim)
plot(rank(fitted(mod))/n,res,xlab="fitted values (rank-transformed)",ylab="residuals")
boxplot(res~yr,data=focal.dat,ylab="residuals",xlab="")

par(mfrow=c(1,2))

phi.obs<-phi(mod)
phi.sim<-vector()
for(i in 1:nsim){
  print(i)
  ysim<-simulate(mod)$sim_1
  phi.sim[i]<-phi(glmmTMB(ysim ~ yr + (1 | station.code), family=nbinom1, data=focal.dat))
}
phi.p<-mean(phi.sim>=phi.obs)
hmax<-max(hist(phi.sim,plot=F)$density)
hist(phi.sim,main=paste("Dispersion (p=",round(phi.p,2),")",sep=""),freq=F,xlab="",yaxt="n",ylab="",xlim=range(phi.sim,phi.obs))
segments(phi.obs,0,phi.obs,hmax,col="red",lwd=2)

pzero.obs<-mean(focal.dat$maps==0)
pzero.sim<-vector()
for(i in 1:nsim){ pzero.sim[i]<-mean(simulate(m3)$sim_1==0) }
pzero.p<-mean(pzero.sim>=pzero.obs)
hmax<-max(hist(pzero.sim,plot=F)$density)
hist(pzero.sim,main=paste("Proportion of zeros (p=",round(pzero.p,2),")",sep=""),freq=F,xlab="",yaxt="n",ylab="",xlim=range(pzero.sim,pzero.obs))
segments(pzero.obs,0,pzero.obs,hmax,col="red",lwd=2)

## This looks like a fine model

## compare to null with nbinom1
m0 <- glmmTMB(maps ~ 1 + (1 | station.code), data = focal.dat, family = nbinom1)


AICtab(mod, m0, logLik = TRUE, base = TRUE, mnames = c("mod", "m0"))

## More support for model with year as factor then null model

summary(mod)

# What if yr is not a factor but scaled continuous?

focal.sp <- "Dark-eyed Junco"
focal.dat <- d.maps %>% filter(species == focal.sp) %>% mutate(yr = as.numeric(scale(year)))

mod <- glmmTMB(maps ~ yr + (1 | station.code), data = focal.dat, family = nbinom1)


summary(mod)

## qq-plot with pearson diagnositics
par(mfrow=c(1,1))
res<-residuals(mod,type="pearson")
qqnorm(res,xlab="expected",ylab="observed",main="")
abline(0,1,lwd=1.5)

### diagnostic plots using PIT residuals

res<-residuals(simulateResiduals(mod,n=nsim,seed=F))
qqunif(res,nsim)
plot(rank(fitted(test))/n,res,xlab="fitted values (rank-transformed)",ylab="residuals")
boxplot(res~yr,data=focal.dat,ylab="residuals",xlab="")
stripchart(res~focal.dat$yr,vertical=T,xlab="",ylab="residuals",main="",pch=21,xlim=c(0,n.yr)+0.5)


phi.obs<-phi(mod)
phi.sim<-vector()
for(i in 1:nsim){
  print(i)
  ysim<-simulate(mod)$sim_1
  phi.sim[i]<-phi(glmmTMB(ysim ~ yr + (1 | station.code), family=nbinom1, data=focal.dat))
}
phi.p<-mean(phi.sim>=phi.obs)
hmax<-max(hist(phi.sim,plot=F)$density)
hist(phi.sim,main=paste("Dispersion (p=",round(phi.p,2),")",sep=""),freq=F,xlab="",yaxt="n",ylab="",xlim=range(phi.sim,phi.obs))
segments(phi.obs,0,phi.obs,hmax,col="red",lwd=2)

pzero.obs<-mean(focal.dat$maps==0)
pzero.sim<-vector()
for(i in 1:nsim){ pzero.sim[i]<-mean(simulate(m3)$sim_1==0) }
pzero.p<-mean(pzero.sim>=pzero.obs)
hmax<-max(hist(pzero.sim,plot=F)$density)
hist(pzero.sim,main=paste("Proportion of zeros (p=",round(pzero.p,2),")",sep=""),freq=F,xlab="",yaxt="n",ylab="",xlim=range(pzero.sim,pzero.obs))
segments(pzero.obs,0,pzero.obs,hmax,col="red",lwd=2)



summary(mod)

(emm1<-emmeans(mod,~yr,type="lp"))
(emm2<-emmeans(mod,~yr,type="response"))

### plot estimates and confidence intervals

(df.emm1<-as.data.frame(emm1))
(df.emm2<-as.data.frame(emm2))

par(mfrow=c(1,1))
ylims<-range(na.omit(df.emm2$asymp.LCL,df.emm2$asymp.UCL))
stripchart(df.emm2$response~c(1:19),vertical=T,ylim=ylims,xlab="",ylab="",main="",pch=19,xlim=c(0,19)+0.5)
segments(1:n.yr,df.emm2$asymp.LCL,1:n.yr,df.emm2$asymp.UCL)

newdata <- data.frame(yr = unique(focal.dat$yr),
                      station.code = rep("CV6", length(unique(focal.dat$year))))

focal.dat$station.code

fit <- predict(mod, newdata ,type = "response", se.fit = TRUE)

?predict

mean.maps <- focal.dat %>% group_by(year) %>% reframe(mean = mean(maps))

upp <- fit$fit + 1.96*fit$se.fit
low <- fit$fit - 1.96 * fit$se.fit

ylims <- range(upp, low)
plot(x = c(2005:2019, 2022, 2023), y = fit$fit, type = "l", ylim = ylims)
points(x = c(2005:2019, 2022, 2023), y = mean.maps$mean, add = TRUE)
lines(x = c(2005:2019, 2022, 2023), y = fit$fit + 1.96 * fit$se.fit, lty = 2)
lines(x = c(2005:2019, 2022, 2023), y = fit$fit - 1.96 * fit$se.fit, lty = 2)


## What if I include all years? 
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

d.maps <- srmAbund %>%
  mutate(station.code = as.factor(station.code),
         park = as.factor(park), 
         loop = as.factor(loop))

focal <- "Dark-eyed Junco"
focal.dat <- d.maps %>% filter(species == focal) %>% mutate(yr = as.numeric(scale(year))) %>% select(park, loop, station.code, maps, year, yr)
yr <- focal.dat$yr
n.yr <- length(unique(yr))

modall <- glmmTMB(maps ~ yr + (1 | station.code), data = focal.dat, family = nbinom1)
summary(mod)
summary(modall)


res<-residuals(simulateResiduals(modall,n=nsim,seed=F))
qqunif(res,nsim)

(emm1<-emmeans(modall,~yr,type="lp"))
(emm2<-emmeans(modall,~yr,type="response"))

summary(modall)

coef(modall)
exp(.57)
exp(-1.7)
exp(-3.64)


yrs <- unique(focal.dat$yr)
stations <- unique(focal.dat$station.code)

newdata <- data.frame(yr = rep(yrs, length(stations)),
                      station.code = rep(stations, length(yrs)))



fit <- predict(modall, newdata = newdata, type = "response", se.fit = TRUE)
str(fit)
head(fit)

preds <- cbind(newdata, fit$fit, fit$se.fit)

mean.preds <- preds %>% group_by(yr) %>% reframe(mean.fit = mean(`fit$fit`), mean.se = mean(`fit$se.fit`))


upp <- mean.preds$mean.fit + 1.96*mean.preds$mean.se
low <- mean.preds$mean.fit - 1.96*mean.preds$mean.se

ylims <- range(upp, low)

plot(x = c(1996:2024), y = mean.preds$mean.fit, ylim = ylims, type = "l")
points(mean.maps$year, mean.maps$mean)
lines(x = c(1996:2024), y = upp, lty = 2)
lines(x = c(1996:2024), y = low, lty = 2)

effect <- (exp(summary(modall)$coefficients[1]$cond[2, 1]/sd(focal.dat$year)) - 1) * 100


ggplot(mean.maps, aes(x = year, y = mean)) +
  geom_point(color = bcs_colors["dark green"]) +
  geom_line(aes(x = c(1996:2024), y = mean.preds$mean.fit), color = bcs_colors["bright green"]) +
  geom_line(aes(x = c(1996:2024), y = upp), color = bcs_colors["bright green"], linetype = 2) +
  geom_line(aes(x = c(1996:2024), y = low), color = bcs_colors["bright green"], linetype = 2) +
  labs(x = "Year", y = "Mean abundance per survey", title = "Mean abundance per survey", subtitle = "Dark-eyed Junco") +
  annotate("text", x =2002, y = 0.31, label = paste0("Mean abundance per survey increasing by ", 
                                                   round(effect, 1), "% per year"),
           size = 5, color = bcs_colors["dark green"]) +
  theme_bcs()

exp(summary(modall)$coefficients[1]$cond[2, 1]/sd(focal.dat$year))
.0811/sd(focal.dat$year)

emm2
.01*29

for(i in 1:length(spp)){
  df <- d.maps[d.maps$species == spp[i], ]
  mod <- glmer(maps ~ yr + (1 | station.code) + (1 | yr) + (1 | station.code:yr), data = df, family = "poisson")
  rate[i] <- exp(summary(mod)$coefficients[2, 1])
  pval[i] <- summary(mod)$coefficients[2, 4]
}


## combine species, rates, and pvalues into data table
maps.rates <- data.frame(species = spp,
                         maps.rate = rate,
                         maps.pval = pval)


## identify statistically significant and top 10 fastest declining species by mean abund per suv
maps.decline <- maps.rates %>%
  filter(maps.pval < 0.05, maps.rate < 1)

## Species proportion of surveys reporting trends
spp <- sort(unique(srmDet$species)) # list of unique species in maps data set to use for indexing
rate <- numeric(length = length(spp))
pval <- numeric(length = length(spp))

for(i in 1:length(spp)){
  dets <- cbind(srmDet[srmDet$species == spp[i], ]$det, srmDet[srmDet$species == spp[i], ]$notdet)  
  yrs <- srmDet[srmDet$species == spp[i], ]$year - min(srmDet[srmDet$species == spp[i], ]$year) + 1
  mod <- glm(dets ~ yrs, family = binomial(link = "logit"))
  rate[i] <- ilogit(coef(mod)[2])
  pval[i] <- coef(summary(mod))[2, 4]
}


prop.rates <- data.frame(species = spp,
                         prop.rate = rate,
                         prop.pval = pval)

prop.decline.20 <- prop.rates %>%
  filter(prop.pval < 0.05) %>%
  slice_min(., prop.rate, n = 20)

topDecliners <- inner_join(maps.decline.20, prop.decline.20)
sort(unique(prop.decline.20$species))
sort(unique(maps.decline.20$species))

# we now have our declining species list
concern <- sort(unique(topDecliners$species)) 

ggplot(filter(d.maps, species == "Violet-green Swallow"), aes(x = year, y = maps)) +
  geom_smooth(method = "loess") 

ggplot(filter(srmDet, species == "Greater Scaup"), aes(x = year, y = det / (det + notdet))) +
  geom_smooth(method = "loess", se = FALSE, aes(color = species))# +
  #theme_bcs()


ggplot(filter(d.maps, species == "American Coot"), aes(x = year, y = maps)) +
  geom_point() +
  geom_smooth(method = lm)


amro <- lm(maps ~ year, data = filter(d.maps, species == "American Robin"))
summary(amro)
overall_p(coot)
coef(coot)[2]

## Maybe we could also look at maximum flock size per year?

maxFlock <- d %>%
  mutate(count = seen + fly + heard) %>%
  group_by(year, species) %>%
  summarise(max = max(count)) %>%
  ungroup() %>%
  pivot_wider(names_from = species, values_from = max) %>%
  replace(is.na(.), 0)

d.max <- maxFlock %>%
  pivot_longer(2:length(maxFlock), names_to = "species", values_to = "max")

spp <- sort(unique(d.max$species)) # list of unique species in maps data set to use for indexing
rate <- numeric(length = length(spp))
pval <- numeric(length = length(spp))

for(i in 1:length(spp)){
  df <- d.max[d.max$species == spp[i], ]
  mod <- lm(max ~ year, data = df)
  rate[i] <- coef(mod)[2]
  pval[i] <- overall_p(mod)
}

max.rates <- data.frame(species = spp,
                         max.rate = rate,
                         max.pval = pval)

max.decline.20 <- max.rates %>%
  filter(max.pval < 0.05) %>%
  slice_min(., max.rate, n = 20)

max.decline.20

ggplot(filter(d.max, species %in% max.decline.20$species), aes(x = year, y = max)) +
  geom_smooth(method = "loess", se = FALSE) +#, color = bcs_colors['dg']) + 
  facet_wrap(~species, scales = "free") #+
#  theme_bcs()

