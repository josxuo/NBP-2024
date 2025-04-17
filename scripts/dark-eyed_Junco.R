#### Common Raven at Discovery Park

library(tidyverse)
library(readxl)

nbp <- read_xlsx("data/processed/nbp_tidy_jan_24.xlsx")
covs <- read_csv("data/processed/circ_no_overlap_covariates.csv")

nsurv <- nbp %>% mutate(station.year = paste(station.code, year, sep = "-")) %>%
  #filter(month %in% c(5, 6)) %>% 
  group_by(station.year) %>%
  summarise(nsurv = n_distinct(survey_id))

focal <- nbp %>% filter(#park %in% c("Magnuson Park", "Discovery Park"), 
                        #month %in% c(5, 6),
                        year %in% c(2005:2019, 2021:2023)) %>%
  group_by(year, park, station.code, bird.code) %>%
  summarise(count = sum(seen, heard, fly), .groups = "drop") %>%
  pivot_wider(names_from = bird.code, values_from = count, values_fill = list(count = 0)) %>%
  pivot_longer(-c(1:3), names_to = "bird.code", values_to = "count") %>%
  filter(bird.code == "VASW") %>%
  mutate(station.year = paste(station.code, year, sep = "-"),
         yr = as.numeric(scale(year)),
         park = as.factor(park)) %>%
  left_join(., nsurv)

head(focal)

plot(focal$year, focal$count)

focal[is.na(focal$station.code), ]$station.code <- "UNK"


library(glmmTMB)

mod.null <- glmmTMB(count ~ 1, offset = log(nsurv), 
                    zi = ~ 1, family = nbinom1, data = focal)
mod1 <- glmmTMB(count ~ yr, offset = log(nsurv), 
               zi = ~ 1, family = nbinom1(), data = focal)

mod2 <- glmmTMB(count ~ yr, offset = log(nsurv),
               zi = ~ 1, family = nbinom2(), data = focal)
mod <- glmmTMB(count ~ yr + (1 | station.code), offset = log(nsurv),
                zi = ~ yr + (1 |station.code), family = nbinom1(), data = focal)

AIC(mod.null, mod1, mod2, mod)


exp(-0.3365/sd_year)


summary(mod)

plot(simulateResiduals(mod))
testZeroInflation(mod)

focal$precipitation <- as.factor(focal$precipitation)

set.seed(123)  # For reproducibility

# Observed data
observed_counts <- c(`Light Rain` = 362, None = 3178, Rain = 289, Snow = 24)
total_count <- sum(observed_counts)

# Convert to proportions
proportions <- observed_counts / total_count

# Sample size (change this to your desired sample size)
sample_size <- 3100  

# Generate a new sample
sampled_levels <- sample(names(observed_counts), size = sample_size, replace = TRUE, prob = proportions)

# Create a data frame
df_sampled <- data.frame(Condition = factor(sampled_levels, levels = names(observed_counts)))

# Check proportions in the new sample
table(df_sampled$Condition) / sample_size  # Should be close to the original proportions

# View first few rows
head(df_sampled)

library(DHARMa)

AIC(mod.null, mod)

plot(simulateResiduals(mod))
testZeroInflation(mod)
summary(mod)

## diagnostics look OK

# Create a new data frame for prediction
predict_data <- expand.grid(
  yr = seq(min(focal$yr), max(focal$yr), length.out = 100),
  station.code = levels(as.factor(focal$station.code)),
  #park = levels(as.factor(focal$park)),
  nsurv = mean(focal$nsurv)
)

predict_data$precipitation <- df_sampled$Condition

# Generate predictions on the link scale (logit)
predicted <- predict(mod, newdata = predict_data, type = "link", se.fit = TRUE)


# Convert predictions to a data frame and transform back to response scale
predicted_df <- cbind(predict_data, mu = predicted$fit, se = predicted$se.fit) %>%
  mutate(
    lower = exp(mu - 1.96 * se),
    upper = exp(mu + 1.96 * se),
    mu = exp(mu)
  )

sd_year <- sd(focal$year)
mean_year <- mean(focal$year)

predicted_df$year <- predicted_df$yr * sd_year + mean_year

mean_predicted <- predicted_df %>%
  group_by(yr) %>%
  mutate(mean_mu = mean(mu),
         mean_upper = mean(upper), 
         mean_lower = mean(lower))

mean_observed <- focal %>% group_by(year) %>%
  summarise(mean_count = mean(count))


ggplot(data = mean_observed, aes(x = park, y = mean_count)) +
  geom_point() +
  #geom_smooth(method = "loess", color = bcs_colors["peach"], se = FALSE) +
  #geom_point(data = mean_predicted, aes(x = park, y = mean_mu)) +
  geom_errorbar(data = mean_predicted, aes(x = park, y = mean_mu, ymin = mean_lower, ymax = mean_upper), 
                inherit.aes = FALSE) +
  theme_bcs()


ggplot(data = mean_observed, aes(x = year, y = mean_count)) +
  #geom_line(data = mean_predicted, aes(x = year, y = mean_mu, color = "Modeled trend"),
           # linetype = "dashed") +
  #geom_ribbon(data = mean_predicted, aes(x = year, ymin = mean_lower, ymax = mean_upper, fill = "95% Confidence interval"), 
             # alpha = 0.2, inherit.aes = FALSE) +
  geom_point(color = bcs_colors["dark green"]) +
  scale_color_manual(name = element_blank())+
                     #values = unname(chart_colors))+
  #scale_fill_manual(name = element_blank(),
  #                  values = c("95% Confidence Interval" = bcs_colors["dark green"])) +
  labs(title = "AVERAGE ANNUAL COUNT OF VAUX'S SWIFT", subtitle = "Per count station at NBP Parks", 
       x = "Year", y = "Average annual count") + 
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
nsurv <- nbp %>% filter(month %in% c(5, 6)) %>%
  mutate(station.code = as.factor(ifelse(is.na(station.code), "Unk", station.code)),
                                  station.year = paste(station.code, year, sep = "-")) %>%
  group_by(station.year) %>%
  summarise(nsurv = n_distinct(survey_id))

focal<- nbp %>% filter(year %in% c(2005:2019, 2021:2023),
                       park %in% c("Discovery Park", "Magnuson Park")) %>%
  group_by(year, station.code, bird.code) %>%
  summarise(count = sum(seen, heard, fly), .groups = "drop") %>%
  pivot_wider(names_from = bird.code, values_from = count, values_fill = list(count = 0)) %>%
  pivot_longer(-c(1, 2), names_to = "bird.code", values_to = "count") %>%
  filter(bird.code == "BUOR") %>%
  mutate(station.code = as.factor(ifelse(is.na(station.code), "Unk", station.code)),
         station.year = paste(station.code, year, sep = "-"),
         yr = as.numeric(scale(year))) %>%
  left_join(., nsurv)

focal$station.code <- as.factor(focal$station.code)
levels(focal$station.code)
sum(is.na(focal$station.code))

library(mgcv)
mod <- gam(count ~ yr + s(station.code, bs = "re") + offset(log(nsurv)), 
           family = nb(), # Negative Binomial (for count data)
           data = focal)

mod2 <- gam(count ~ s(yr) + s(park, bs = "re") + offset(log(nsurv)), 
               family = nb(), # Negative Binomial (for count data)
               data = cora)


summary(mod)

# Create a new data frame for prediction
predict_data <- expand.grid(
  yr = seq(min(focal$yr), max(focal$yr), length.out = 100),
  #station.code = levels(as.factor(cora$station.code)),
  station.code = levels(as.factor(focal$station.code)),
  nsurv = mean(focal$nsurv))

predicted <- predict(mod, newdata = predict_data, type = "link", se.fit = TRUE)

# Convert predictions to a data frame and transform back to response scale
predicted_df <- cbind(predict_data, mu = predicted$fit, se = predicted$se.fit) %>%
  mutate(
    lower = exp(mu - 1.96 * se),
    upper = exp(mu + 1.96 * se),
    mu = exp(mu)
  )

sd_year <- sd(focal$year)
mean_year <- mean(focal$year)

predicted_df$year <- predicted_df$yr * sd_year + mean_year

mean_predicted <- predicted_df %>%
  group_by(year) %>%
  mutate(mean_mu = mean(mu),
         mean_upper = mean(upper), 
         mean_lower = mean(lower))

mean_observed <- focal %>% group_by(year) %>%
  mutate(mean_count = mean(count))

levels(mean_observed$park)

ggplot(data = mean_observed, aes(x = year, y = mean_count)) +
  geom_line(data = mean_predicted, aes(x = year, y = mean_mu, color = "Modeled trend"),
            linetype = "dashed") +
  geom_ribbon(data = mean_predicted, aes(x = year, ymin = mean_lower, ymax = mean_upper, fill = "95% Confidence interval"), 
              alpha = 0.2, inherit.aes = FALSE) +
  geom_point() +
  labs(title = "AVERAGE ANNUAL COUNT OF COMMON RAVEN", subtitle = "Per NBP Park", 
       x = "Year", y = "Average annual count") + 
  
  #scale_color_manual(name = "Key",
  #                   values = c("Carkeek Park" = "#36BA3A",
                                #"Cheasty Greenspace" = "green",
  #                              "Discovery Park" = "#0A3C23",
                                #"Golden Gardens Park" = "steelblue",
                                #"Genesee Park" = "yellow",
                                #"Lake Forest Park" = "orange",
   #                             "Lincoln Park" = "#FF8E54",
                                #"Magnuson Park" = "blue",
    #                            "Seward Park" = "#88D498",
     #                           "Washington Park Arboretum" = "#FFD700",
      #                          "Modeled trend" = "#FFB98C"),
       #              limits = c("Carkeek Park",
                                #"Cheasty Greenspace",
        #                        "Discovery Park",
                                #"Golden Gardens Park",
                                #"Genesee Park",
                                #"Lake Forest Park",
         #                       "Lincoln Park",
                                #"Magnuson Park",
          #                      "Seward Park",
           #                     "Washington Park Arboretum",
            #                    "Modeled trend")) + 
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


## Aerial insectivores

func <- read_xlsx("data/processed/nbp_species_traits.xlsx")

AI <- func %>% select(bird.code, AI)

nbp.ai <- nbp %>% left_join(., AI) %>% mutate(
  station.code = as.factor(ifelse(is.na(station.code), "Unk", station.code)))


nsurv <- nbp.ai %>% group_by(year, station.code) %>%
  summarise(nsurv = n_distinct(survey_id), .groups = "drop") %>%
  mutate(station.year = paste(station.code, year, sep = "-")) %>%
  select(station.year, nsurv)

srm_ai <- nbp.ai %>% filter(year %in% c(2005:2019, 2021:2023)) %>%
  group_by(year, park, station.code, AI) %>%
  summarise(count = sum(seen + heard + fly), .groups = "drop") %>%
  pivot_wider(names_from = AI, values_from = count, values_fill = list(count = 0)) %>%
  pivot_longer(-c(1:3), names_to = "AI", values_to = "count") %>%
  filter(AI == "AI") %>%
  mutate(yr = as.numeric(scale(year)),
         station.year = paste(station.code, year, sep = "-"),
         station.code2 = ifelse(is.na(station.code), "Unk", station.code)) %>%
  left_join(., nsurv, join_by("station.year" == "station.year"))

levels(srm_ai$station.code)
view(srm_ai)

str(srm_ai)

mod.null <- glmmTMB(count ~ 1, offset = log(nsurv), family = poisson, data = srm_ai)
mod1 <- glmmTMB(count ~ yr, offset = log(nsurv), family = poisson, data = srm_ai)
mod2 <- glmmTMB(count ~ yr + (1 | station.code), offset = log(nsurv), family = poisson, data = srm_ai)
mod3 <- glmmTMB(count ~ yr + park + (1 | station.code), offset = log(nsurv), family = poisson, data = srm_ai)
mod4 <- glmmTMB(count ~ yr + park + (1 | station.code), offset = log(nsurv), family = nbinom1(), data = srm_ai)

AIC(mod.null, mod1, mod2, mod3, mod4)

summary(mod3)


plot(simulateResiduals(mod4))
testZeroInflation(mod4)

summary

# Create a new data frame for prediction
predict_data <- expand.grid(
  yr = seq(min(srm_ai$yr), max(srm_ai$yr), length.out = 100),
  park = levels(as.factor((srm_ai$park))),
  station.code = levels(as.factor((srm_ai$station.code))),
  nsurv = mean(srm_ai$nsurv))

predict_data$park <- as.factor(predict_data$park)
predict_data$station.code <- as.factor(predict_data$station.code)

predict_data$nsurv
class(predict_data)

predicted <- predict(mod4, newdata = predict_data, type = "link", se.fit = TRUE)

# Convert predictions to a data frame and transform back to response scale
predicted_df <- cbind(predict_data, mu = predicted$fit, se = predica

sd_year <- sd(focal$year)
mean_year <- mean(focal$year)

predicted_df$year <- predicted_df$yr * sd_year + mean_year

mean_predicted <- predicted_df %>%
  group_by(year) %>%
  mutate(mean_mu = mean(mu),
         mean_upper = mean(upper), 
         mean_lower = mean(lower))

mean_observed <- focal %>% group_by(year) %>%
  mutate(mean_count = mean(count))

levels(mean_observed$park)

ggplot(data = mean_observed, aes(x = year, y = mean_count)) +
  geom_line(data = mean_predicted, aes(x = year, y = mean_mu, color = "Modeled trend"),
            linetype = "dashed") +
  geom_ribbon(data = mean_predicted, aes(x = year, ymin = mean_lower, ymax = mean_upper, fill = "95% Confidence interval"), 
              alpha = 0.2, inherit.aes = FALSE) +
  geom_point() +
  labs(title = "AVERAGE ANNUAL COUNT OF COMMON RAVEN", subtitle = "Per NBP Park", 
       x = "Year", y = "Average annual count") + 
  
  #scale_color_manual(name = "Key",
  #                   values = c("Carkeek Park" = "#36BA3A",
  #"Cheasty Greenspace" = "green",
  #                              "Discovery Park" = "#0A3C23",
  #"Golden Gardens Park" = "steelblue",
  #"Genesee Park" = "yellow",
  #"Lake Forest Park" = "orange",
  #                             "Lincoln Park" = "#FF8E54",
  #"Magnuson Park" = "blue",
  #                            "Seward Park" = "#88D498",
  #                           "Washington Park Arboretum" = "#FFD700",
  #                          "Modeled trend" = "#FFB98C"),
  #              limits = c("Carkeek Park",
  #"Cheasty Greenspace",
  #                        "Discovery Park",
  #"Golden Gardens Park",
  #"Genesee Park",
  #"Lake Forest Park",
  #                       "Lincoln Park",
  #"Magnuson Park",
  #                      "Seward Park",
  #                     "Washington Park Arboretum",
  #                    "Modeled trend")) + 
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
