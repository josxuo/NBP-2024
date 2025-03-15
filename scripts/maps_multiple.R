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
data <- nbp %>%
  filter(!str_detect(species, pattern = " sp\\."),
         station.code %in% covs$Station, 
         year %in% c(2012:2019, 2022, 2023),
         exclusions == "none") %>%
  mutate(count = sum(seen, heard, fly))

# inspect data
## sort(unique(d$species))
## sort(unique(d$year))
## sort(unique(d$code))



focal <- c("DEJU", "ANHU", "NOFL")


dat <- data %>%
  group_by(survey_id, station.code, park, bird.code, year, month, weather, wind, precipitation) %>%
  reframe(count = sum(seen, heard, fly)) %>%  # Use sum() directly
  pivot_wider(names_from = bird.code, values_from = count, values_fill = list(count = 0)) %>%
  pivot_longer(cols = -c(survey_id:precipitation), names_to = "bird.code", values_to = "count") %>%
  mutate(year = year - min(year, na.rm = TRUE)) %>%  # Ensure no NA issues %>%
  filter(bird.code %in% focal)

head(dat)

N <- nrow(dat)
S <- length(unique(dat$bird.code))
L <- length(unique(dat$station.code))
Nmonths <- length(unique(dat$month))
counts <- dat$count
species <- as.numeric(as.factor(dat$bird.code))
site <- as.numeric(as.factor(dat$station.code))
month <- dat$month
year <- dat$year

data_list <- list(N = N, S = S, L = L, Nmonths = Nmonths, counts = counts, species = species, 
             site = site, month = month, year = year)

nchains <- 2
niters <- 5000
nburn <- niters / 2

mod <- jags.model("code/models/trend_all_species.R", data = data_list, n.chains = nchains)
monitor <- c("beta1", "sigma_site", "r", "sigma_month", "month_effect", "ysim")
samp <- coda.samples(mod, monitor, n.iter = niters)

colnames(samp[[1]])
ind <- grep("beta1|sigma_site|month_effect", colnames(samp[[1]]))
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
for(i in 1:12){
  den1 <- density(res1[, i])
  den2 <- density(res2[, i])
  ylims <- range(den1$y, den2$y)
  xlims <- range(den1$x, den2$x)
  plot(den1$x, den1$y, type = "l", ylim = ylims, xlim = xlims, xlab = "", ylab = "", main = colnames(res1)[i])
  lines(den2$x, den2$y, col = "red")
}

res <- rbind(res1, res2)

par(mfrow = c(5, 2))
for(i in 1:12){
  den <- density(res[, i])
  plot(den$x, den$y, type = "l", xlab = "", ylab = "", main = colnames(res1)[i])
}


pairs(res[1:10^3,ind])
corr<-cor(res[,ind])
corr<-corr-diag(diag(corr))
round(range(corr,na.rm=T),2)

summarise((exp(res[, ind]) - 1) * 100)
rm(mod1)
parks



posterior_samples <- res[, grep("ysim", colnames(res))]
coun

hist(posterior_samples, col = "lightblue", main = "Posterior Predictive Check")
abline(v = mean(observed_values), col = "red", lwd = 2)  # Compare with observed mean
library(stats)

# Compute mean predictions from simulated data
y_pred <- apply(posterior_samples, 2, mean)  # Expected y from the posterior
phi <- mean(res[, "r"])

neg_bin_deviance <- function(y, mu, phi) {
  -2 * sum(dnbinom(y, size = phi, mu = mu, log = TRUE))
}


# Compute discrepancy for observed data
D_obs <- neg_bin_deviance(counts, mu = y_pred, phi)

# Compute discrepancy for replicated data
D_rep <- apply(posterior_samples, 1, function(y_rep) neg_bin_deviance(y_rep, mu = y_pred, phi))

# Compute Bayesian p-value
bayesian_p_value <- mean(D_rep > D_obs)
print(bayesian_p_value)


saveRDS(mod, file = "eight_species_trends.rds")
mod <- readRDS("eight_species_trends.rds")

summary(mod)

model_samples <- as.mcmc(mod)
summary(model_samples)

str(mod)

#1 - AMCR: Decrease; -3.4%
#2 - AMRO: Decrease; -3.3%
#3 - ANHU: Increase; +2.1%
#4 - BAEA: Decrease; -4.3%
#5 - DEJU: Increase; +8.2%
#6 - NOFL: Decrease; -5.1%
#7 - SOSP: No trend
#8 - STJA: No trend




exp()

library(MASS)

model_DEJU <- glm.nb(count ~ year + as.factor(park) + as.factor(month) + as.factor(weather) + as.factor(wind) + as.factor(precipitation), data = subset(dat, bird.code == "DEJU"))



# Subset data
dat_DEJU <- subset(dat, bird.code == "DEJU") 




# Add the fitted values (predictions) to the summary data
dat_DEJU$predicted_count <- predict(model_DEJU, newdata = dat_DEJU, type = "response")

dat_DEJU_summary <- dat_DEJU %>%
  group_by(year) %>%
  reframe(mean_count = mean(count, na.rm = TRUE),
          mean_predicted = mean(predicted_count, na.rm = TRUE))  # Calculate mean count per year

# Plot the actual vs. predicted trend
ggplot(dat_DEJU_summary, aes(x = year)) +
  geom_point(aes(y = mean_count), size = 3, color = "blue") +  # Observed data (mean per year)
  geom_line(aes(y = mean_predicted), color = "red", size = 1) +  # Predicted trend line
  labs(title = "Trend of Dark-eyed Junco Mean Detections Over Time (Model Fit)",
       x = "Year", y = "Mean Count per Survey") +
  theme_minimal()
