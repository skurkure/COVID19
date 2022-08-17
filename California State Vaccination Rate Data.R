library(tidyverse)
library(ggplot2)
library(EpiEstim)
library(deSolve)
install.packages("Dict")
library(Dict)
county_covid_data <- read.csv("covid19cases_test.csv")
county_covid_data <- county_covid_data %>%
  mutate(State = "California")
head(county_covid_data)
unique(county_covid_data$area)
#the three extras are California, Out of State, and Unknown
sums <- county_covid_data %>%
  group_by(date) %>% 
  summarize(state_cum_cases = sum(cumulative_cases))
head(sums, 50)


cases <- dict()
#then add in all the dates in the key column of the dictionary
case_date_aggregator <-for i in county_covid_data$date:
  if i.date == cases.date():
    i.cumulative_cases += cases.cumulative_cases()


