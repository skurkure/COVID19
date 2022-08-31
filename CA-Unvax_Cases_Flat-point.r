library(tidyverse)
library(EpiEstim)
library(deSolve)


#read file
covidvax <- read.csv("COVID19/covid19postvaxstatewidestats.csv")
head(covidvax, 1)


covidvax <- covidvax %>%
  mutate(Submission_Date = as.Date(date))

covidvax[is.na(covidvax)] <- 0


head(covidvax$unvaccinated_cases,5)


agg_unvax_cases <- aggregate(unvaccinated_cases~Submission_Date, data=covidvax, sum)
agg_unvax_cases


flat_agg_unvax_cases <- agg_unvax_cases %>%
  filter(Submission_Date >= "2021-03-01" & Submission_Date <= "2021-04-30") %>%
  pull(unvaccinated_cases)
flat_agg_unvax_cases


agg_unvax_cases %>% ggplot(aes(x=as.Date(Submission_Date))) + geom_point(aes(y=unvaccinated_cases))
flat_agg_unvax_cases_only <- agg_unvax_cases %>%
  mutate(date = as.Date(Submission_Date)) %>%
  filter(date >= "2021-03-01" & date <= "2021-04-30")

flat_agg_unvax_cases_only %>% ggplot(aes(x=as.Date(Submission_Date))) + geom_point(aes(y=unvaccinated_cases))
flat_agg_unvax_cases
flat_agg_unvax_cases[15]-flat_agg_unvax_cases[14]
flat_agg_unvax_cases[30]-flat_agg_unvax_cases[29]
flat_lm <- lm(unvaccinated_cases~date, data=flat_agg_unvax_cases_only)

flat_agg_unvax_cases_only <- flat_agg_unvax_cases_only  %>%
  mutate(date = as.Date(Submission_Date))
flat_agg_unvax_cases_only <- flat_agg_unvax_cases_only %>%
  mutate(time=1:61)
flat_agg_unvax_cases_only %>%
  ggplot(aes(x=time, y=unvaccinated_cases)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE)
unvaccinated_cases <- flat_agg_unvax_cases_only$unvaccinated_cases
time <- 1:61
summary(lm(unvaccinated_cases~time))$coefficients

