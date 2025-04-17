#### Common Raven at Discovery Park


nbp <- read_xlsx("data/processed/nbp_tidy_jan_24.xlsx")
covs <- read_csv("data/processed/circ_no_overlap_covariates.csv")

nsurv <- dat %>% mutate(station.year = paste(station.code, year, sep = "-")) %>%
  group_by(station.year) %>%
  summarise(nsurv = n_distinct(survey_id))

disco.cora <- dat %>% filter(park == "Discovery Park",
                             year %in% c(2005:2019, 2021:2023)) %>%
  group_by(year, station.code, bird.code) %>%
  summarise(count = sum(seen, heard, fly), .groups = "drop") %>%
  pivot_wider(names_from = bird.code, values_from = count, values_fill = list(count = 0)) %>%
  pivot_longer(-c(1, 2), names_to = "bird.code", values_to = "count") %>%
  filter(bird.code == "CORA") %>%
  mutate(station.year = paste(station.code, year, sep = "-"),
         yr = as.numeric(scale(year))) %>%
  left_join(., nsurv)

view(disco.cora)

library(glmmTMB)

mod.null <- glmmTMB(count ~ (1 | station.code), offset = log(nsurv), 
                    zi = ~ 1, family = nbinom1, data = disco.cora)
mod <- glmmTMB(count ~ yr + (1 | station.code), offset = log(nsurv), 
               zi = ~ yr + (1 | station.code), family = nbinom1, data = disco.cora)

library(DHARMa)

AIC(mod.null, mod)

plot(simulateResiduals(mod))
testZeroInflation(mod)

## diagnostics look OK

# Create a new data frame for prediction
predict_data <- expand.grid(
  yr = seq(min(disco.cora$yr), max(disco.cora$yr), length.out = 100),
  station.code = levels(as.factor(disco.cora$station.code)),
  nsurv = mean(disco.cora$nsurv)
)

# Generate predictions on the link scale (logit)
predicted <- predict(mod, newdata = predict_data, type = "link", se.fit = TRUE)

# Convert predictions to a data frame and transform back to response scale
predicted_df <- cbind(predict_data, mu = predicted$fit, se = predicted$se.fit) %>%
  mutate(
    lower = exp(mu - 1.96 * se),
    upper = exp(mu + 1.96 * se),
    mu = exp(mu)
  )

sd_year <- sd(disco.cora$year)
mean_year <- mean(disco.cora$year)

predicted_df$year <- predicted_df$yr * sd_year + mean_year

mean_predicted <- predicted_df %>%
  group_by(year) %>%
  mutate(mean_mu = mean(mu),
         mean_upper = mean(upper), 
         mean_lower = mean(lower))

mean_observed <- disco.cora %>% group_by(year) %>%
  mutate(mean_count = mean(count))


ggplot(data = mean_observed, aes(x = year, y = mean_count)) +
  geom_point() +
  geom_smooth(method = "loess", color = bcs_colors["peach"], se = FALSE) +
  #geom_line(data = mean_predicted, aes(x = year, y = mean_mu)) +
  #geom_ribbon(data = mean_predicted, aes(x = year, ymin = mean_lower, ymax = mean_upper), 
  #            alpha = 0.2, fill = bcs_colors["dark green"], inherit.aes = FALSE) +
  theme_bcs()


min(disco.cora[disco.cora$count > 0, ]$year)
min(dat[dat$bird.code == "CORA", ]$survey_date)

filter(dat, bird.code == "CORA") %>% arrange(survey_date)

summary(mod)

exp(0.5560 / sd_year)

library(mgcv)

disco.cora$station.code <- as.factor(disco.cora$station.code)

gam(count ~ yr + (1 | station.code), offset = log(nsurv), 
    zi = ~ yr + (1 | station.code), family = nbinom1, data = disco.cora)

mod <- gam(count ~ s(yr) + s(station.code, bs = "re") + offset(log(nsurv)), 
           family = nb(), # Negative Binomial (for count data)
           data = disco.cora)

summary(mod)

# Create a new data frame for prediction
predict_data <- expand.grid(
  yr = seq(min(disco.cora$yr), max(disco.cora$yr), length.out = 100),
  station.code = levels(as.factor(disco.cora$station.code)),
  nsurv = mean(disco.cora$nsurv))

predicted <- predict(mod, newdata = predict_data, type = "link", se.fit = TRUE)

# Convert predictions to a data frame and transform back to response scale
predicted_df <- cbind(predict_data, mu = predicted$fit, se = predicted$se.fit) %>%
  mutate(
    lower = exp(mu - 1.96 * se),
    upper = exp(mu + 1.96 * se),
    mu = exp(mu)
  )

sd_year <- sd(disco.cora$year)
mean_year <- mean(disco.cora$year)

predicted_df$year <- predicted_df$yr * sd_year + mean_year

mean_predicted <- predicted_df %>%
  group_by(year) %>%
  mutate(mean_mu = mean(mu),
         mean_upper = mean(upper), 
         mean_lower = mean(lower))

mean_observed <- disco.cora %>% group_by(year) %>%
  mutate(mean_count = mean(count))


disco.cora.p <- ggplot(data = mean_observed, aes(x = year, y = mean_count)) +
  geom_line(data = mean_predicted, aes(x = year, y = mean_mu, color = "Modeled trend"),
            linetype = "dashed") +
  geom_ribbon(data = mean_predicted, aes(x = year, ymin = mean_lower, ymax = mean_upper, fill = "95% Confidence interval"), 
              alpha = 0.2, inherit.aes = FALSE) +
  geom_point(aes(color = "Observed counts")) +
  labs(title = "AVERAGE ANNUAL COUNT OF COMMON RAVEN", subtitle = "Per count station at Discovery Park", 
       x = "Year", y = "Average annual count") + 
  
  scale_color_manual(name = "Key",
                    values = c("Observed counts" = "#0A3C23",
                               "Modeled trend" = "#FFB98C"),
                    limits = c("Observed counts",
                               "Modeled trend")) + 
  scale_fill_manual(name = element_blank(),
                     values = c("95% Confidence interval" = "#FFB98C"),
                     limits = c("95% Confidence interval")) +

  guides(color = guide_legend(order = 1),
         fill = guide_legend(order = 2)) +
  theme_bcs() + 
  theme(legend.spacing = unit(-1.17, "lines"),
        legend.key.size = unit(0.25, "cm"), 
        #plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "in"),
        panel.background = element_rect(fill = "transparent", color = NA), # Transparent panel
        plot.background = element_rect(fill = "transparent", color = NA),  # Transparent plot background
        legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend background
        legend.key = element_rect(fill = "transparent", color = NA) # Transparent legend keys
        )

ggsave("results/figures/cora_at_disco.png", bg = "transparent", width = 6, height = 4, units = "in", dpi = 300)


?ggsave

png(filename = "results/figures/cora_at_disco.png", height = 3, width = 5, res = 300, units = "in")
disco.cora.p
dev.off()


## What parks are ravens associated with?
nsurv <- dat %>% mutate(park.year = paste(park, year, sep = "-")) %>%
  group_by(park.year) %>%
  summarise(nsurv = n_distinct(survey_id))

cora <- dat %>% filter(year %in% c(2005:2019, 2021:2023),
                       park %in% c("Carkeek Park", "Seward Park", "Discovery Park", "Washington Park Arboretum",
                                   "Lincoln Park")) %>%
  group_by(year, park, bird.code) %>%
  summarise(count = sum(seen, heard, fly), .groups = "drop") %>%
  pivot_wider(names_from = bird.code, values_from = count, values_fill = list(count = 0)) %>%
  pivot_longer(-c(1, 2), names_to = "bird.code", values_to = "count") %>%
  filter(bird.code == "CORA") %>%
  mutate(park.year = paste(park, year, sep = "-"),
         yr = as.numeric(scale(year))) %>%
  left_join(., nsurv)

cora$park <- as.factor(cora$park)
levels(cora$park)

mod <- gam(count ~ s(yr) + park + offset(log(nsurv)), 
           family = nb(), # Negative Binomial (for count data)
           data = cora)

mod2 <- gam(count ~ s(yr) + s(park, bs = "re") + offset(log(nsurv)), 
               family = nb(), # Negative Binomial (for count data)
               data = cora)




# Create a new data frame for prediction
predict_data <- expand.grid(
  yr = seq(min(cora$yr), max(cora$yr), length.out = 100),
  #station.code = levels(as.factor(cora$station.code)),
  park = levels(as.factor(cora$park)),
  nsurv = mean(cora$nsurv))

predicted <- predict(mod, newdata = predict_data, type = "link", se.fit = TRUE)

# Convert predictions to a data frame and transform back to response scale
predicted_df <- cbind(predict_data, mu = predicted$fit, se = predicted$se.fit) %>%
  mutate(
    lower = exp(mu - 1.96 * se),
    upper = exp(mu + 1.96 * se),
    mu = exp(mu)
  )

sd_year <- sd(cora$year)
mean_year <- mean(cora$year)

predicted_df$year <- predicted_df$yr * sd_year + mean_year

mean_predicted <- predicted_df %>%
  group_by(year) %>%
  mutate(mean_mu = mean(mu),
         mean_upper = mean(upper), 
         mean_lower = mean(lower))

mean_observed <- cora %>% group_by(park, year) %>%
  mutate(mean_count = mean(count))

levels(mean_observed$park)

ggplot(data = mean_observed, aes(x = year, y = mean_count)) +
  geom_line(data = mean_predicted, aes(x = year, y = mean_mu, color = "Modeled trend"),
            linetype = "dashed") +
  geom_ribbon(data = mean_predicted, aes(x = year, ymin = mean_lower, ymax = mean_upper, fill = "95% Confidence interval"), 
              alpha = 0.2, inherit.aes = FALSE) +
  geom_point(aes(color = park), alpha = 0.5) +
  labs(title = "AVERAGE ANNUAL COUNT OF COMMON RAVEN", subtitle = "Per NBP Park", 
       x = "Year", y = "Average annual count") + 
  
  scale_color_manual(name = "Key",
                     values = c("Carkeek Park" = "#36BA3A",
                                #"Cheasty Greenspace" = "green",
                                "Discovery Park" = "#0A3C23",
                                #"Golden Gardens Park" = "steelblue",
                                #"Genesee Park" = "yellow",
                                #"Lake Forest Park" = "orange",
                                "Lincoln Park" = "#FF8E54",
                                #"Magnuson Park" = "blue",
                                "Seward Park" = "#88D498",
                                "Washington Park Arboretum" = "#FFD700",
                                "Modeled trend" = "#FFB98C"),
                     limits = c("Carkeek Park",
                                #"Cheasty Greenspace",
                                "Discovery Park",
                                #"Golden Gardens Park",
                                #"Genesee Park",
                                #"Lake Forest Park",
                                "Lincoln Park",
                                #"Magnuson Park",
                                "Seward Park",
                                "Washington Park Arboretum",
                                "Modeled trend")) + 
  scale_fill_manual(name = element_blank(),
                    values = c("95% Confidence interval" = "#FFB98C"),
                    limits = c("95% Confidence interval")) +
  
  guides(color = guide_legend(order = 1),
         fill = guide_legend(order = 2)) +
  theme_bcs() + 
  theme(legend.spacing = unit(-1.17, "lines"),
        legend.key.size = unit(0.25, "cm"), 
        #plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "in"),
        panel.background = element_rect(fill = "transparent", color = NA), # Transparent panel
        plot.background = element_rect(fill = "transparent", color = NA),  # Transparent plot background
        legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend background
        legend.key = element_rect(fill = "transparent", color = NA) # Transparent legend keys
  )


ggsave("results/figures/cora_at_all_parks.png", bg = "transparent", width = 6, height = 4, units = "in", dpi = 300)

summary(mod)

cora %>% group_by(park) %>%
  summarise(total = sum(count))
