nbp <- read_xlsx("data/processed/nbp_tidy_jan_24.xlsx")
covs <- read.csv("data/processed/circ_no_overlap_covariates.csv")

d2 <- nbp %>%
  filter(!str_detect(species, pattern = " sp\\.|Spotted Owl|Domestic"),
         station.code %in% covs$station.code,
         #!park %in% c("Cheasty Greenspace", "Lincoln Park"),
         year %in% c(2005:2019, 2021, 2022, 2023)) %>%
  # exclusions == "none") %>%
  mutate(loop2 = as.factor(paste(park, loop, sep = "-")),
         park = as.factor(park),
         month = as.factor(month),
         #station.code = as.factor(station.code),
         station.year = paste(station.code, year, sep = "-"),
         park.year = paste(park, year, sep = "-"),
         count = seen + heard + fly)


nurv <- d2 %>% group_by(park.year, month) %>%
  summarise(nsurv = n_distinct(survey_id))


monthly_obs <- d2 %>% group_by(park, year, park.year, month, bird.code) %>%
  summarise(nobs = sum(count), .groups = "drop") %>%
  pivot_wider(names_from = bird.code, values_from = nobs, values_fill = list(nobs = 0)) %>%
  pivot_longer(-c(1, 2, 3, 4), names_to = "bird.code", values_to = "nobs")

focal <- "DEJU"
focal.dat <- monthly_obs %>% 
  filter(bird.code == focal) %>% 
  mutate(yr = as.numeric(scale(year)))

mod <- glm(nobs ~ yr + month + park, data = focal.dat, family = poisson)
plot(simulateResiduals(mod))
testZeroInflation(mod)


mod.pois <- mod.pln <- glmmTMB(nobs ~ month + (yr | park), 
                               data = focal.dat, 
                               family = poisson, 
                               dispformula = ~1)

sqrt(focal.dat$nobs)

mod.nb <- glmmTMB(nobs ~ yr * park + month, data = focal.dat, family = nbinom1())

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  dispersion_ratio <- sum(rp^2) / rdf
  p_value <- pchisq(sum(rp^2), df = rdf, lower.tail = FALSE)
  cat("Dispersion ratio:", dispersion_ratio, "\n")
  cat("p-value:", p_value, "\n")
}
overdisp_fun(mod.nb)

testZeroInflation(mod.nb)
plot(simulateResiduals(mod.nb))
testOutliers(mod.nb)

summary(mod.nb)
summary(mod.nb.re)

mod.nb.zi <- glmmTMB(nobs ~ yr + month + park, zi = ~1, data = focal.dat, family = nbinom1())


mod.nb2 <- glmmTMB(nobs ~ yr + month + park, data = focal.dat, family = nbinom2())

testZeroInflation(mod.nb2)
plot(simulateResiduals(mod.nb2))
summary(mod.nb2)
                       

predict_data <- expand.grid(
  yr = seq(min(focal.dat$yr), max(focal.dat$yr), length.out = 100),  # time range for prediction
  month = levels(focal.dat$month),
  park = levels(focal.dat$park))

predicted <- predict(mod.nb, newdata = predict_data, type = "link", se.fit = TRUE, re.form = NULL)

predicted_df <- cbind(predict_data, predicted$fit, predicted$se.fit) %>%
  rename(mu_log = `predicted$fit`, se_log = `predicted$se.fit`) %>%
  mutate(
    lower_log = mu_log - (1.96 * se_log),
    upper_log = mu_log + (1.96 * se_log),
    mu = exp(mu_log),
    lower = exp(lower_log),
    upper = exp(upper_log)
  )

mean_year <- mean(focal.dat$year)
sd_year <- sd(focal.dat$year)

mean_predicted <- predicted_df %>%
  group_by(yr, park) %>%
  dplyr::summarise(mean_mu = mean(mu),
                   mean_lower = mean(lower),
                   mean_upper = mean(upper), .groups = "drop") %>%
  mutate(year = yr * sd_year + mean_year)

focal_mean <- focal.dat %>% group_by(yr, park) %>%
  mutate(mean_count = mean(nobs)) %>%
  mutate(year = yr * sd_year + mean_year)

focal.species <- "Dark-eyed Junco"
p.sp <- ggplot(focal_mean, aes(x = year, y = mean_count)) +
  geom_point(fill = bcs_colors["dark green"], alpha = 0.5) + 
  geom_line(data = mean_predicted, aes(x = year, y = mean_mu)) +
  geom_ribbon(data = mean_predicted, aes(x = year, ymin = mean_lower, ymax = mean_upper), 
             alpha = 0.2, fill = bcs_colors["dark green"], inherit.aes = FALSE) + 
  labs(title = paste("Average monthly counts of", focal.species, "over time"), x = "Year", y = "Average monthly count") +
  theme_bcs()


png(filename = "results/figures/avg_monthly_deju_counts.png", 
    height = 4, width = 5, units = "in", res = 300)
p.sp
dev.off()

p.sp <- ggplot(mean_predicted, aes(x = (yr * sd_year + mean_year), y = mean_mu, color = park)) + 
  geom_line() +
  geom_point(data = focal_mean,  size = .8, aes(x = year, y = mean_count, color = park)) +
  scale_color_manual(values = unname(chart_colors)) +
  labs(title = "Average count of Dark-eyed Junco observations", 
       x = "Year", y = "Average count",
       color = "Park") + 
  theme_bcs()


summary(mod.nb)

## group by stationID 
nsurv <- d2 %>% group_by(station.year) %>%
  summarise(nsurv = n_distinct(survey_id))

cnt.dat <- d2 %>% group_by(station.year, station.code, park, year, bird.code) %>%
  summarise(ann_count = sum(count), .groups = "drop") %>% 
  pivot_wider(names_from = bird.code, values_from = ann_count, values_fill = list(ann_count = 0)) %>%
  pivot_longer(-c(1:4), names_to = "bird.code", values_to = "ann_count") %>%
  left_join(., nsurv)

view(cnt.dat[1:100, ])

focal.park <- "Carkeek Park"
focal.species <- "CORA"
focal.dat <- cnt.dat %>% filter(park == focal.park, bird.code == focal.species) %>%
  mutate(yr = as.numeric(scale(year))) %>% left_join(., nsurv)

mod.null <- glmmTMB(ann_count ~ (1 | station.code), offset = log(nsurv), data = focal.dat, family = poisson)
mod <- glmmTMB(ann_count ~ yr + (1 | station.code), offset = log(nsurv), 
               zi = ~ 1, data = focal.dat, family = nbinom1())

summary(mod)

library(DHARMa)

plot(simulateResiduals(mod))
testZeroInflation(mod)


# Create a new data frame for prediction
predict_data <- expand.grid(
    yr = seq(min(focal.dat$yr), max(focal.dat$yr), length.out = 100),
    station.code = levels(as.factor(focal.dat$station.code))
  )

predict_data$nsurv <- round(mean(focal.dat$nsurv),0)

# Generate predictions on the link scale (log)
predicted <- predict(mod, newdata = predict_data, type = "link", se.fit = TRUE)



# Convert predictions to a data frame and transform back to response scale
predicted_df <- cbind(predict_data, mu = predicted$fit, se = predicted$se.fit) %>%
  mutate(
    lower = exp(mu - 1.96 * se),
    upper = exp(mu + 1.96 * se),
    mu = exp(mu)
  )

mean_predicted <- predicted_df %>% group_by(yr) %>%
  summarise(mean_mu = mean(mu),
            mean_lower = mean(lower),
            mean_upp = mean(upper),
            .groups = "drop")

mean_year <- mean(focal.dat$year)
sd_year <- sd(focal.dat$year)

focal.mean <- focal.dat %>% group_by(year) %>%
  summarise(mean_cnt = mean(ann_count))

ggplot() +
  geom_point(data = focal.mean, aes(x = year, y = mean_cnt)) +
  geom_line(data = mean_predicted, aes(x = (yr * sd_year + mean_year), y = mean_mu)) +
  geom_ribbon(data = mean_predicted, aes(x = (yr * sd_year + mean_year), ymin = mean_lower, ymax = mean_upp), 
              alpha = 0.2, fill = bcs_colors["dark green"], inherit.aes = FALSE) +
  theme_bcs() 


focal.parks <- c("Carkeek Park", "Discovery Park", "Genesee Park", "Golden Gardens Park", 
                 "Lincoln Park", "Magnuson Park", "Seward Park", "Washington Park Arboretum")

r2 <- numeric()
#r2_marginal <- numeric()
disp <- numeric()
beta_yr <- numeric()
se_yr <- numeric()
park <- vector()

res_list <- list()

for(i in 1:length(focal.parks)) {
  
  park.dat <- cnt.dat[cnt.dat$park == focal.parks[i], ]
  
  spp <- park.dat %>% group_by(year, bird.code) %>% summarise(dets = sum(ann_count > 0), .groups = "drop") %>%
    group_by(bird.code) %>% summarise(total_dets = sum(dets), .groups = "drop") %>% filter(total_dets >= 10) %>% 
    select(bird.code) %>% unique() %>% arrange(bird.code) %>% pull()
  
  r2 <- numeric(length(spp))
  disp <- numeric(length(spp))
  beta_yr <- numeric(length(spp))
  se_yr <- numeric(length(spp))
  park <- character(length(spp))
  species <- character(length(spp))
  
  
  for(j in 1:length(spp)){
    
    
    
    mod.dat <- park.dat[park.dat$bird.code == spp[j],] %>% mutate(yr = as.numeric(scale(year)))
    mod <- glmmTMB(ann_count ~ yr + (1 | station.code), offset = log(nsurv),
                   zi = ~ 1, data = mod.dat, family = nbinom1(link = "log"))
    r2[j] <- r2(mod)
    disp[j] <- sum(residuals(mod, type = "pearson")^2) / df.residual(mod)
    beta_yr[j] <- summary(mod)$coefficients$cond["yr", 1]
    se_yr[j] <- summary(mod)$coefficients$cond["yr", 2]
    p_yr[j] <- summary(mod)$coefficients$cond["yr", 4]
    park[j] <- focal.parks[i]
    species[j] <- spp[j]
    print(paste(focal.parks[i], spp[j], "complete"))
    
  }

  res_list[[focal.parks[i]]] <- data.frame(park, species, beta_yr, se_yr, p_yr, r2, disp)
  print(paste(focal.parks[i], "complete"))

}

## UPDATED CODE
res_list <- list()

for(i in 1:length(focal.parks)) {
  
  park.dat <- cnt.dat[cnt.dat$park == focal.parks[i], ]
  
  spp <- park.dat %>% 
    group_by(year, bird.code) %>% 
    summarise(dets = sum(ann_count > 0), .groups = "drop") %>%
    group_by(bird.code) %>% 
    summarise(total_dets = sum(dets), .groups = "drop") %>% 
    filter(total_dets >= 10) %>% 
    arrange(bird.code) %>% 
    pull(bird.code)
  
  # Preallocate storage
  r2 <- numeric(length(spp))
  disp <- numeric(length(spp))
  beta_yr <- numeric(length(spp))
  se_yr <- numeric(length(spp))
  p_yr <- numeric(length(spp))
  park <- character(length(spp))
  species <- character(length(spp))
  
  for(j in seq_along(spp)){
    
    mod.dat <- park.dat %>% filter(bird.code == spp[j]) %>% mutate(yr = as.numeric(scale(year)))
    
    mod <- glmmTMB(ann_count ~ yr + (1 | station.code), offset = log(nsurv),
                   zi = ~ 1, data = mod.dat, family = nbinom1(link = "log"))
    
    # Compute R2
    r2_vals <- r2_nakagawa(mod)
    if (is.list(r2_vals)) {
      r2[j] <- r2_vals$R2_marginal
    } else {
      r2[j] <- NA  # Fallback if the function fails
    }
    
    # Dispersion calculation
    disp[j] <- sum(residuals(mod, type = "pearson")^2) / df.residual(mod)
    
    # Extract model coefficients
    beta_yr[j] <- summary(mod)$coefficients$cond["yr", "Estimate"]
    se_yr[j] <- summary(mod)$coefficients$cond["yr", "Std. Error"]
    p_yr[j] <- summary(mod)$coefficients$cond["yr", "Pr(>|z|)"]
    
    # Store identifiers
    park[j] <- focal.parks[i]
    species[j] <- spp[j]
    
    print(paste(focal.parks[i], spp[j], "complete"))
  }
  
  # Save results
  res_list[[focal.parks[i]]] <- data.frame(park, species, beta_yr, se_yr, p_yr, r2, disp)
  print(paste(focal.parks[i], "complete"))
}


count_trends <- bind_rows(res_list, .id = "source")

write.csv(count_trends, "results/tables/count_trends_by_species_by_park.csv")

view(count_trends)

count_trends <- read.csv("results/tables/count_trends_by_species_by_park.csv")

spp.list <- list()

for(i in 1:length(focal.parks)){
  park.dat <- cnt.dat[cnt.dat$park == focal.parks[i], ]
  
  spp.list[[focal.parks[i]]] <- park.dat %>% 
    group_by(year, bird.code) %>% 
    summarise(dets = sum(ann_count > 0), .groups = "drop") %>%
    group_by(bird.code) %>% 
    summarise(total_dets = sum(dets), .groups = "drop") %>% 
    filter(total_dets >= 1) %>% 
    arrange(bird.code) %>% 
    pull(bird.code)
  
  print(paste(focal.parks[i], "complete"))
  
}

combined_spp.list <- bind_rows(lapply(names(spp.list), function(park_name) {
  data.frame(park = park_name, species = spp.list[[park_name]]) %>%
    mutate(park.species = paste(park, species, sep = "-")) %>%
    select(park.species)
}))


count_trends$park.species <- paste(count_trends$park, count_trends$species, sep = "-")

combined_trends <- left_join(combined_spp.list, count_trends) %>% 
  mutate(status = factor(case_when(
    beta_yr >= 0.1 & p_yr < 0.05 & r2 >= 0.2 & disp < 2 ~ "Increasing",
    beta_yr > 0 & p_yr < 0.1 & r2 >= 0.1 & disp < 2 ~ "Possibly increasing",
    beta_yr <= -0.1 & p_yr < 0.05 & r2 >= 0.2 & disp < 2 ~ "Decreasing",
    beta_yr < 0 & p_yr < 0.1 & r2 >= 0.1 & disp < 2 ~ "Possibly decreasing",
    is.na(beta_yr) ~ "Data deficient",
    TRUE ~ "Trend uncertain"), levels = c("Increasing", "Possibly increasing",
                                            "Trend uncertain", "Possibly decreasing",
                                            "Decreasing")))
    


z_mults <- data.frame(park = character(),
                      mean_yr = numeric(),
                      sd_yr = numeric())

for(i in 1:length(focal.parks)){
  z_mults[i, 1] <- focal.parks[i]
  z_mults[i, 2] <- mean(cnt.dat[cnt.dat$park == focal.parks[i], ]$year)
  z_mults[i, 3] <- sd(cnt.dat[cnt.dat$park == focal.parks[i], ]$year)
}

plots <- list()
for(i in 1:length(focal.parks)){
  plot.dat <- combined_trends %>% filter(park == focal.parks[i]) %>%
    mutate(sp2 = fct_reorder(species, beta_yr))
  
  plots[[i]] <- ggplot(plot.dat, 
                      aes(y = sp2, 
                          x = (exp(beta_yr / z_mults[z_mults$park == focal.parks[i], "sd_yr"][[1]]) - 1),
                  fill = status)) + 
    geom_col(color = "black", linetype = "solid", linewidth = 0.3) +
    scale_fill_manual(values = c(
      "Trend uncertain" = "#FAF5F0",
      "Decreasing" = "coral",
      "Increasing" = "#36BA3A",
      "Possibly decreasing" = "#FFB98C",
      "Possibly increasing" = "#72D672"
    )) +
    theme_bcs() +
    labs(x = "Estimated annual percent change", y = "Species code", fill = "Trend description",
         title = "Suggestive species abundance trends", subtitle = focal.parks[i]) +
    theme(#axis.text.x = element_text(angle = 45, hjust = 1),
          #panel.grid.major.y = element_blank(),
          legend.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"),
          legend.key = element_rect(color = NA, fill = NA),
          #legend.key.size = unit(0.5, "cm"),
          legend.position = "inside",
          legend.justification = c(1, 0),
          legend.position.inside = c(1, 0),
          legend.background = element_rect(fill = "white", color = "#0A3C23")
          )
}


print(plots[[1]])


for(i in 1:length(focal.parks)){
  png(filename = paste0("results/figures/species_abundance_trends_", 
                      focal.parks[i],".png"),
    width = 4, height = 10, units = "in", res = 300)
  print(plots[[i]])
  dev.off()
}


plots[[2]]
head(plots)

str(z_mults)


park.dat <- cnt.dat[cnt.dat$park == focal.parks[1], ] 

spp <- park.dat %>% group_by(year, bird.code) %>% summarise(dets = sum(ann_count > 0), .groups = "drop") %>%
  group_by(bird.code) %>% summarise(total_dets = sum(dets)) %>% filter(total_dets >= 10) %>% select(bird.code) %>% unique() %>% arrange() %>% pull()
  
  mod.dat <- park.dat[park.dat$bird.code == spp[1],] %>% mutate(yr = as.numeric(scale(year)))
  
  mod <- glmmTMB(ann_count ~ yr + (1 | station.code), offset = log(nsurv),
                 zi = ~ 1, data = mod.dat, family = nbinom1(link = "log"))
  
  r2_conditional[i] <- r2_nakagawa(mod)$R2_conditional
  r2_marginal[i] <- r2_nakagawa(mod)$R2_marginal
  disp[i] <- sum(residuals(mod, type = "pearson")^2) / nrow(mod.dat) - length(fixef(mod)$cond)
  beta_yr[i] <- summary(mod)$coefficients$cond["yr", 1]
  se_yr[i] <- summary(mod)$coefficients$cond["yr", 2]
  p_yr[i] <- summary(mod)$coefficients$cond["yr", 4]
  park[i] <- focal.parks[i]



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



