library(tidyverse)
library(EpiEstim)
library(deSolve)


#read file
vaccination_rate_data <- read.csv("COVID19/covid19vaccinesbycounty.csv")
agg_vax <- aggregate(cumulative_fully_vaccinated~administered_date, data = vaccination_rate_data, sum)
agg_vax
omicron_agg_vax <- agg_vax %>%
  filter(administered_date >= "2022-01-01" & administered_date <= "2022-02-28") %>%
  pull(cumulative_fully_vaccinated)
omicron_agg_vax


agg_vax %>% ggplot(aes(x=as.Date(administered_date))) + geom_point(aes(y=cumulative_fully_vaccinated))
omicron_agg_vax_only <- agg_vax %>%
  mutate(date = as.Date(administered_date)) %>%
  filter(administered_date >= "2022-01-01" & administered_date <= "2022-02-28")

omicron_agg_vax_only %>% ggplot(aes(x=as.Date(administered_date))) + geom_point(aes(y=cumulative_fully_vaccinated))
omicron_agg_vax
omicron_agg_vax[15]-omicron_agg_vax[14]
omicron_agg_vax[30]-omicron_agg_vax[29]
omicron_lm <- lm(cumulative_fully_vaccinated~administered_date, data=omicron_agg_vax_only)

omicron_agg_vax_only <- omicron_agg_vax_only  %>%
  mutate(date = as.Date(administered_date))
omicron_agg_vax_only <- omicron_agg_vax_only %>%
  mutate(time=1:59)
omicron_agg_vax_only %>%
  ggplot(aes(x=time, y=cumulative_fully_vaccinated)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE)
cumulative_fully_vaccinated <- omicron_agg_vax_only$cumulative_fully_vaccinated
time <- 1:59
summary(lm(cumulative_fully_vaccinated~time))$coefficients
#57763.01 per day
1/57763.01
#1.731212e-05 is epsilon
1.731212e-05*0.95
#1.644651e-05
1.731212e-05*0.05
#8.65606e-07
