library(tidyverse)
library(EpiEstim)
library(deSolve)


#read file
covidvax <- read.csv("COVID19/covid19postvaxstatewidestats.csv")
head(covidvax, 1)


covidvax <- covidvax %>%
  mutate(Submission_Date = as.Date(date))

covidvax[is.na(covidvax)] <- 0


head(covidvax$vaccinated_deaths,5)


agg_vax_death <- aggregate(vaccinated_deaths~Submission_Date, data=covidvax, sum)
agg_vax_death


flat_agg_vax_death <- agg_vax_death %>%
  filter(Submission_Date >= "2021-03-01" & Submission_Date <= "2021-04-30") %>%
  pull(vaccinated_deaths)
flat_agg_vax_death


agg_vax_death %>% ggplot(aes(x=as.Date(Submission_Date))) + geom_point(aes(y=vaccinated_deaths))
flat_agg_vax_death_only <- agg_vax_death %>%
  mutate(date = as.Date(Submission_Date)) %>%
  filter(date >= "2021-03-01" & date <= "2021-04-30")

flat_agg_vax_death_only %>% ggplot(aes(x=as.Date(Submission_Date))) + geom_point(aes(y=vaccinated_deaths))
flat_agg_vax_death
flat_agg_vax_death[15]-flat_agg_vax_death[14]
flat_agg_vax_death[30]-flat_agg_vax_death[29]
flat_lm <- lm(vaccinated_deaths~date, data=flat_agg_vax_death_only)

flat_agg_vax_death_only <- flat_agg_vax_death_only  %>%
  mutate(date = as.Date(Submission_Date))
flat_agg_vax_death_only <- flat_agg_vax_death_only %>%
  mutate(time=1:61)
flat_agg_vax_death_only %>%
  ggplot(aes(x=time, y=vaccinated_deaths)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE)
vaccinated_deaths <- flat_agg_vax_death_only$vaccinated_deaths
time <- 1:61
summary(lm(vaccinated_deaths~time))$coefficients

