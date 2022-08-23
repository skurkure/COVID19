library(tidyverse)
library(EpiEstim)
library(deSolve)


#read file
vaccination_rate_data <- read.csv("COVID19/covid19vaccinesbycounty.csv")
agg_vax <- aggregate(cumulative_fully_vaccinated~administered_date, data = vaccination_rate_data, sum)
agg_vax
flat_agg_vax <- agg_vax %>%
  filter(administered_date >= "2021-03-01" & administered_date <= "2021-04-30") %>%
  pull(cumulative_fully_vaccinated)
flat_agg_vax


agg_vax %>% ggplot(aes(x=as.Date(administered_date))) + geom_point(aes(y=cumulative_fully_vaccinated))
flat_agg_vax_only <- agg_vax %>%
  mutate(date = as.Date(administered_date)) %>%
  filter(administered_date >= "2021-03-01" & administered_date <= "2021-04-30")

flat_agg_vax_only %>% ggplot(aes(x=as.Date(administered_date))) + geom_point(aes(y=cumulative_fully_vaccinated))
flat_agg_vax
flat_agg_vax[15]-flat_agg_vax[14]
flat_agg_vax[30]-flat_agg_vax[29]
flat_lm <- lm(cumulative_fully_vaccinated~administered_date, data=flat_agg_vax_only)

flat_agg_vax_only <- flat_agg_vax_only  %>%
  mutate(date = as.Date(administered_date))
flat_agg_vax_only <- flat_agg_vax_only %>%
  mutate(time=1:61)
flat_agg_vax_only %>%
  ggplot(aes(x=time, y=cumulative_fully_vaccinated)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE)
cumulative_fully_vaccinated <- flat_agg_vax_only$cumulative_fully_vaccinated
time <- 1:61
summary(lm(cumulative_fully_vaccinated~time))$coefficients

#532550 per day
1/532550
#epsilon = 1.877758e-06
1.87778e-06 * 0.95
1.877758e-06 * 0.05


