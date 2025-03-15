library(rjags)
load.module("glm")
library(tidyverse)
library(readxl)

summarise<-function(x){
  if(!is.matrix(x)){ x<-as.matrix(x) }
  low<-apply(x,2,quantile,probs=0.025)
  upp<-apply(x,2,quantile,probs=0.975)
  mean<-apply(x,2,mean)
  out<-cbind(mean,low,upp)
  colnames(out)<-c("mean","lower","upper")
  signif(out,digits=4)
}

dat <- read_xlsx("data/b_intermediate_data/nbp_tidy_jan_24.xlsx")
circ.lu <- read_csv("data/c_analysis_data/nbp_circ_codes.csv") %>% na.exclude() %>% mutate(sort.col = row_number(.))
circ.filt <- read_csv("data/c_analysis_data/circ_no_overlap_covariates.csv")

dat <- left_join(dat, circ.lu)

nsurv <- dat %>%
  filter(year > 2004, 
         !year %in% c(2020, 2021, 2024)) %>%
  group_by(code, year, month) %>%
  reframe(nsurv = n_distinct(survey_id)) %>%
  pivot_wider(names_from = code, values_from = nsurv) %>%
  arrange(year, month) %>%
  replace(is.na(.), 0) %>%
  filter(month %in% c(4, 5, 6))

nsurv

## not really enough. Maybe look at just a few months?

years <- c(2004:2019, 2022, 2023)


d <- dat %>%
  filter(year > 2004,
         code %in% circ.filt$Station,
         !year %in% years,
         month == 12,
         !str_detect(species, pattern = " sp.| x "))

complete <- d %>%
  group_by(code, year) %>%
  reframe(nsurv = n_distinct(survey_id)) %>%
  pivot_wider(names_from = year, values_from = nsurv) %>%
  filter(!is.na(code)) %>%
  replace(is.na(.), 0) %>%
  mutate(row_sum = apply(.[,-1], 1, sum)) %>%
  filter(row_sum == 17) %>%
  pull(code)


complete
## Jan = 9, Feb = 24, March = 27, April = 39, May = 41, June = 35, July = 28
## Aug = 8, Sept = 21, Oct = 9, Nov = 14, Dec = 8

#### MAY IS THE WINNER! ####

d <- dat %>%
  filter(year > 2004,
         code %in% complete,
         !year %in% c(2020, 2021, 2024),
         month == 5,
         !str_detect(species, pattern = " sp.| x "))

S <- d %>%
  group_by(code, year) %>%
  reframe(count = sum(seen, heard, fly), S = n_distinct(species))


ggplot(S, aes(y = count, group = year)) +
  geom_boxplot()

ggplot(S, aes(y = S, group = year)) +
  geom_boxplot()

ggplot(S, aes(x = year, y = log(count))) + 
  geom_point() +
  geom_smooth(method = "lm")

S$yr <- S$year - min(S$year) + 1

count.mod <- lm(log(count) ~ yr, data = S)
summary(count.mod)

(exp(2.515 -0.0239) - exp(2.515 -0.0239*17)) / exp(2.515 -0.0239)  ## model shows may counts are down 31%

S2 <- d %>%
  group_by(year) %>%
  mutate(obs = seen + heard + fly) %>%
  reframe(meanCount = mean(obs), S = n_distinct(species))

S2$yr <- S2$year - min(S2$year) + 1

count.mod2 <- lm(log(meanCount) ~ yr, data = S2)

ggplot(S2, aes(x = year, y = log(meanCount))) +
  geom_point() + 
  geom_smooth(method = "lm")

summary(count.mod2)

# now need a sample of units with no change

dat3 <- d %>%
  filter(year %in% c(2019, 2023),
         month %in% c(4, 5, 6),
         code %in% hrcd$Station,
         !str_detect(species, pattern = " sp.| x "))

completeness3 <- dat3 %>%
  group_by(code, year) %>%
  reframe(nsurv = n_distinct(survey_id)) %>%
  pivot_wider(names_from = year, values_from = nsurv) %>%
  filter(`2019` == 3 & `2023` == 3,
         !code %in% completeness2$code)

no.change.codes <- completeness3[sample(dim(completeness3)[1], 12), ]$code

sus <- c(change.codes, no.change.codes)