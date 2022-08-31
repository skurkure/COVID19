library(tidyverse)
library(EpiEstim)
library(deSolve)


#read file
covidvax <- read.csv("COVID19/covid19postvaxstatewidestats.csv")
head(covidvax, 1)


covidvax <- covidvax %>%
  mutate(Submission_Date = as.Date(date))

covidvax[is.na(covidvax)] <- 0


head(covidvax$vaccinated_hosp,5)


agg_vax_hosp <- aggregate(vaccinated_hosp~Submission_Date, data=covidvax, sum)
agg_vax_hosp


flat_agg_vax_hosp <- agg_vax_hosp %>%
  filter(Submission_Date >= "2021-03-01" & Submission_Date <= "2021-04-30") %>%
  pull(vaccinated_hosp)
flat_agg_vax_hosp


agg_vax_hosp %>% ggplot(aes(x=as.Date(Submission_Date))) + geom_point(aes(y=vaccinated_hosp))
flat_agg_vax_hosp_only <- agg_vax_hosp %>%
  mutate(date = as.Date(Submission_Date)) %>%
  filter(date >= "2021-03-01" & date <= "2021-04-30")

flat_agg_vax_hosp_only %>% ggplot(aes(x=as.Date(Submission_Date))) + geom_point(aes(y=vaccinated_hosp))
flat_agg_vax_hosp
flat_agg_vax_hosp[15]-flat_agg_vax_hosp[14]
flat_agg_vax_hosp[30]-flat_agg_vax_hosp[29]
flat_lm <- lm(vaccinated_hosp~date, data=flat_agg_vax_hosp_only)

flat_agg_vax_hosp_only <- flat_agg_vax_hosp_only  %>%
  mutate(date = as.Date(Submission_Date))
flat_agg_vax_hosp_only <- flat_agg_vax_hosp_only %>%
  mutate(time=1:61)
flat_agg_vax_hosp_only %>%
  ggplot(aes(x=time, y=vaccinated_hosp)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE)
vaccinated_hosp <- flat_agg_vax_hosp_only$vaccinated_hosp
time <- 1:61
summary(lm(vaccinated_hosp~time))$coefficients

