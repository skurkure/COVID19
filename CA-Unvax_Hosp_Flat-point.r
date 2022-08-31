library(tidyverse)
library(EpiEstim)
library(deSolve)


#read file
covidvax <- read.csv("COVID19/covid19postvaxstatewidestats.csv")
head(covidvax, 1)


covidvax <- covidvax %>%
  mutate(Submission_Date = as.Date(date))

covidvax[is.na(covidvax)] <- 0


head(covidvax$unvaccinated_hosp,5)


agg_unvax_hosp <- aggregate(unvaccinated_hosp~Submission_Date, data=covidvax, sum)
agg_unvax_hosp


flat_agg_unvax_hosp <- agg_unvax_hosp %>%
  filter(Submission_Date >= "2021-03-01" & Submission_Date <= "2021-04-30") %>%
  pull(unvaccinated_hosp)
flat_agg_unvax_hosp


agg_unvax_hosp %>% ggplot(aes(x=as.Date(Submission_Date))) + geom_point(aes(y=unvaccinated_hosp))
flat_agg_unvax_hosp_only <- agg_unvax_hosp %>%
  mutate(date = as.Date(Submission_Date)) %>%
  filter(date >= "2021-03-01" & date <= "2021-04-30")

flat_agg_unvax_hosp_only %>% ggplot(aes(x=as.Date(Submission_Date))) + geom_point(aes(y=unvaccinated_hosp))
flat_agg_unvax_hosp
flat_agg_unvax_hosp[15]-flat_agg_unvax_hosp[14]
flat_agg_unvax_hosp[30]-flat_agg_unvax_hosp[29]
flat_lm <- lm(unvaccinated_hosp~date, data=flat_agg_unvax_hosp_only)

flat_agg_unvax_hosp_only <- flat_agg_unvax_hosp_only  %>%
  mutate(date = as.Date(Submission_Date))
flat_agg_unvax_hosp_only <- flat_agg_unvax_hosp_only %>%
  mutate(time=1:61)
flat_agg_unvax_hosp_only %>%
  ggplot(aes(x=time, y=unvaccinated_hosp)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE)
unvaccinated_hosp <- flat_agg_unvax_hosp_only$unvaccinated_hosp
time <- 1:61
summary(lm(unvaccinated_hosp~time))$coefficients

