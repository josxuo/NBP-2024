d <- nbp %>%
  group_by(survey_id, park, year, bird.code) %>%
  summarize(count = sum(seen, heard, fly), .groups = "drop") %>%
  pivot_wider(names_from = bird.code, values_from = count, values_fill = 0) %>%
  pivot_longer(-c("year", "park", "survey_id"), names_to = "bird.code", values_to = "count") %>%
  mutate(park.year = paste(park, year, sep = "-"))


survey_totals <- nbp %>%
  group_by(year, park) %>%
  summarize(nsurv = n_distinct(survey_id), .groups = "drop") %>%
  mutate(park.year = paste(park, year, sep = "-")) %>%
  select(park.year, nsurv)


agg <- d %>% group_by(bird.code, year, park, park.year) %>%
  summarise(total_count = sum(count), .groups = "drop") %>%
  left_join(., survey_totals, join_by("park.year" == "park.year"))


focal = "ANHU"
focal.dat <- filter(agg, bird.code == focal) %>%
  mutate(yr = factor(year - min(year) + 1, levels = 1:28)) %>%
  arrange(park, year)

max(focal.dat$year) - min(focal.dat$year)
  
  
mod.null <- glmmTMB(total_count ~ 1, 
                    family = poisson, 
                    offset = log(nsurv),
                    data = focal.dat)

mod1 <- glmmTMB(total_count ~ yr, 
                    family = poisson, 
                    offset = log(nsurv),
                    data = focal.dat)

mod2 <- glmmTMB(total_count ~ yr + park, 
                family = poisson, 
                offset = log(nsurv),
                data = focal.dat)

mod3 <- glmmTMB(total_count ~ yr + (1 | park), 
                        family = poisson, 
                        offset = log(nsurv),
                        data = focal.dat)

mod4 <- glmmTMB(total_count ~ yr + park + (1 | park), 
                family = poisson, 
                offset = log(nsurv),
                data = focal.dat)

AICs <- AIC(mod.null, mod1, mod2, mod3, mod4)
AICs


best_mod <- get(rownames(AICs)[which.min(AICs$AIC)])

dispersion_test <- simulateResiduals(mod2)
plot(dispersion_test)

mod2.nb1 <- glmmTMB(total_count ~ yr + park, 
                family = nbinom1(link = "log"), 
                offset = log(nsurv),
                data = focal.dat)

mod2.nb2 <- glmmTMB(total_count ~ yr + park, 
                    family = nbinom2(link = "log"), 
                    offset = log(nsurv),
                    data = focal.dat)


AIC(mod2, mod2.nb1, mod2.nb2)




library(nlme)

mod <- glmmTMB(total_count ~ yr + park + (1 | park) + ar1((yr - 1) | park), 
               family = nbinom2(link = "log"), 
               offset = log(nsurv),
               data = focal.dat)
lme4::glmer.nb(total_count ~ yr + )


dispersion_test <- simulateResiduals(mod4.nb2)
plot(dispersion_test)

testZeroInflation(dispersion_test)

temp <- focal.dat %>%
  arrange(park, year)

## check for temporal autocorrelation

n <- dim(focal.dat)[1]
res1<-residuals(simulateResiduals(mod4.nb2, n=1000, seed=F))
plot(res1[1:(n-1)],res1[2:n],xlab="residual at year t",ylab="residual at year t+1")

acf(residuals(mod4.nb2), main = "Autocorrelation of residuals")

acf(residuals(mod4), main = "Autocorrelation of residuals")

library(mgcv)
mod.gam <- gam(total_count ~ s(yr, k = 5) + park,
               family = nb(),
               offset = log(nsurv),
               data = focal.dat)

acf(residuals(mod.gam), main = "Autocorrelation of residuals")
gam.check(mod.gav)
gam.check(mod.gam)

mod4.nb2.ar <- glmmTMB(total_count ~ yr + park + (1 | park),
                       family = nbinom2(link = "log"), 
                       offset = log(nsurv),
                       dispformula = ~1,
                       data = focal.dat,
                       control = glmmTMBControl(optCtrl = list(iter.max = 1e5, eval.max = 1e5)),
                       REML = FALSE)  

# Add AR1 structure manually by including yr as a random slope (approximate AR1)
mod4.nb2.ar <- update(mod4.nb2.ar, . ~ . + (1 | park/yr))

acf(residuals(mod4.nb2.ar), main = "Autocorrelation of residuals")

mod.gam.ar <- gam(total_count ~ s(yr, k=4) + park,
                  family = nb(),
                  offset = log(nsurv),
                  method = "REML",
                  correlation = corAR1(form = ~ yr | park),
                  data = focal.dat)
acf(residuals(mod.gam.ar), main = "Autocorrelation of residuals")

# Calculate autocorrelation per park
acf_results <- focal.dat %>%
  group_by(park) %>%
  do(acf_result = acf(.$total_count, plot = FALSE)) %>%
  ungroup()

# Plot ACF for each park

# Plot ACF for each park
focal.dat %>%
  group_by(park) %>%
  do({
    acf_res <- acf(.$total_count, plot = FALSE)
    
    # Create a data frame with lags and ACF values
    acf_df <- data.frame(
      lag = acf_res$lag,
      acf = acf_res$acf
    )
    
    # Print the plot instead of returning it
    print(ggplot(acf_df, aes(x = lag, y = acf)) +
            geom_line() +
            labs(title = paste("ACF for", unique(.$park)), x = "Lag", y = "Autocorrelation") +
            theme_minimal())
  })
