library(tidyverse)
library(EpiEstim)
library(deSolve)


#read file
covidvax <- read.csv("COVID19/covid19postvaxstatewidestats.csv")
head(covidvax, 1)


covidvax <- covidvax %>%
  mutate(Submission_Date = as.Date(date))

covidvax[is.na(covidvax)] <- 0


head(covidvax$vaccinated_cases,5)


agg_vax_cases <- aggregate(vaccinated_cases~Submission_Date, data=covidvax, sum)
agg_vax_cases


flat_agg_vax_cases <- agg_vax_cases %>%
  filter(Submission_Date >= "2021-03-01" & Submission_Date <= "2021-04-30") %>%
  pull(vaccinated_cases)
flat_agg_vax_cases


agg_vax_cases %>% ggplot(aes(x=as.Date(Submission_Date))) + geom_point(aes(y=vaccinated_cases))
flat_agg_vax_cases_only <- agg_vax_cases %>%
  mutate(date = as.Date(Submission_Date)) %>%
  filter(date >= "2021-03-01" & date <= "2021-04-30")

flat_agg_vax_cases_only %>% ggplot(aes(x=as.Date(Submission_Date))) + geom_point(aes(y=vaccinated_cases))
flat_agg_vax_cases
flat_agg_vax_cases[15]-flat_agg_vax_cases[14]
flat_agg_vax_cases[30]-flat_agg_vax_cases[29]
flat_lm <- lm(vaccinated_cases~date, data=flat_agg_vax_cases_only)

flat_agg_vax_cases_only <- flat_agg_vax_cases_only  %>%
  mutate(date = as.Date(Submission_Date))
flat_agg_vax_cases_only <- flat_agg_vax_cases_only %>%
  mutate(time=1:61)
flat_agg_vax_cases_only %>%
  ggplot(aes(x=time, y=vaccinated_cases)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE)
vaccinated_cases <- flat_agg_vax_cases_only$vaccinated_cases
time <- 1:61
summary(lm(vaccinated_cases~time))$coefficients

