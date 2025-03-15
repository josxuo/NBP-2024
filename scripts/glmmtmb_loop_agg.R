d <- nbp %>%
  filter(!str_detect(species, pattern = " sp\\."),
         station.code %in% covs$Station, 
         year %in% c(2005:2019, 2022, 2023)) %>%
  mutate(loop2 = as.factor(paste(park, loop, sep = "-")),
         park = as.factor(park),
         month = as.factor(month),
         station.code = as.factor(station.code),
         park.year = paste(park, year, sep = "-"))


# inspect data
## sort(unique(d$species))
## sort(unique(d$year))
## sort(unique(d$code))

loop.table <- d %>%
  select(loop2, park) %>%
  unique()

# Set up tables/calculations useful later
## number of surveys conducted each year
nsurvYearsp <- d %>%
  group_by(year, park) %>%
  reframe(nsurv = n_distinct(survey_id)) %>%
  mutate(park.year = paste(park, year, sep = "-")) %>%
  select(park.year, nsurv)
nsurvYearsp



## species p/a matrix
presence <- d %>%
  # filter(station.code %in% sus, year %in% year, month %in% months) %>%
  group_by(year, bird.code, park, park.year) %>%
  reframe(detected = n()) %>%
  pivot_wider(names_from = bird.code, values_from = detected, values_fill = 0) %>%
  pivot_longer(-grep("year|park", colnames(.)), names_to = "bird.code", values_to = "detected") %>%
  left_join(., nsurvYearsp, join_by("park.year" == "park.year")) %>%
  mutate(not_detected = nsurv - detected)


focal <- "NOFL"
focal.dat <- presence %>% filter(bird.code == focal) %>%
                                 #!park %in% c("Cheasty Greenspace", "Washington Park Arboretum")) %>% #,
                                 #!park.year %in% outliers_park.year) %>%
  mutate(yr = year - min(year), 
         park = factor(park))


mod.null <- glmmTMB(cbind(detected, not_detected) ~ 1,
                    family = binomial(link = "logit"), 
                    data = focal.dat)

mod1 <- glmmTMB(cbind(detected, not_detected) ~ yr,
                family = binomial(link = "logit"), 
                data = focal.dat)

mod2 <- glmmTMB(cbind(detected, not_detected) ~ yr + park,
                family = binomial(link = "logit"), 
                data = focal.dat)

mod2.beta <- glmmTMB(cbind(detected, not_detected) ~ yr + park,
                     family = betabinomial(link = "logit"), 
                     data = focal.dat)

mod3 <- glmmTMB(cbind(detected, not_detected) ~ yr + (1 | park),
                family = binomial(link = "logit"),
                data = focal.dat)

mod4 <- glmmTMB(cbind(detected, not_detected) ~ yr + (yr | park),
                family = binomial(link = "logit"),
                data = focal.dat)

mod5 <- glmmTMB(cbind(detected, not_detected) ~ yr + park + (1 | yr),
                family = betabinomial(link = "logit"), 
                data = focal.dat)


AIC(mod.null, mod1, mod2, mod3, mod4, mod5, mod2.beta)  ## most support for model 2. check diagnostics
AIC(mod2, mod2.beta)

dispersion_test <- simulateResiduals(mod5)
plot(dispersion_test)

testZeroInflation(dispersion_test)

## check for autocorrelation
n <- dim(focal.dat)[1]
res1<-residuals(simulateResiduals(mod5,n=1000,seed=F))
plot(res1[1:(n-1)],res1[2:n],xlab="residual at year t",ylab="residual at year t+1")

summary(mod5)

log_odds_to_prob <- function(log_odds) {
  return(1 / (1 + exp(-log_odds)))
}


# Create a new data frame for prediction
predict_data <- expand.grid(
  yr = seq(min(focal.dat$yr), max(focal.dat$yr), length.out = 100),
  park = levels(focal.dat$park)
)

# Generate predictions on the link scale (logit)
predicted <- predict(mod5, newdata = predict_data, type = "link", se.fit = TRUE)

# Convert predictions to a data frame and transform back to response scale
predicted_df <- cbind(predict_data, mu = predicted$fit, se = predicted$se.fit) %>%
  mutate(
    lower = plogis(mu - 1.96 * se),
    upper = plogis(mu + 1.96 * se),
    mu = plogis(mu)
  )

# Average predictions across parks for a general trend
mean_predicted <- predicted_df %>%
  group_by(yr) %>%
  summarize(mean_mu = mean(mu),
            mean_lower = mean(lower),
            mean_upper = mean(upper))

# Compute mean observed proportion per year
focal_mean <- focal.dat %>%
  mutate(prop = detected / nsurv) %>%
  group_by(yr) %>%
  summarize(mean_prop = mean(prop, na.rm = TRUE))

# Plot observed data and model predictions
ggplot(focal_mean, aes(x = (yr + 2005), y = mean_prop)) +
  geom_point(alpha = .3, size = 2) + 
  geom_line(data = mean_predicted, aes(x = (yr + 2005), y = mean_mu), color = "blue") +
  geom_ribbon(data = mean_predicted, aes(x = (yr + 2005), ymin = mean_lower, ymax = mean_upper), 
              alpha = 0.2, fill = "blue", inherit.aes = FALSE) + 
  labs(title = paste("Proportion of surveys reporting", focal, "over time"),
       x = "Year", y = "Proportion of surveys reporting") +
  theme_minimal()
