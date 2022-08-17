library(tidyverse)
library(magrittr)
library(lubridate)
library(tibble)
library(ggplot2)
library(ggthemes)
library(hrbrthemes)
library(rvest)
library(gt)
library(deSolve)
library(EpiEstim)
library(incidence)
library(distcrete)
library(epitrix)
library(projections)
library(dse)
library(rmutil)
#set directory
#note on the data: Note: On November 29, 2021, the denominator for calculating vaccine coverage has been changed from age 16+ to age 12+ to reflect new vaccine eligibility criteria.
#The previous dataset based on age 16+ denominators has been uploaded as an archived table.
#The California Department of Public Health (CDPH) is identifying vaccination status of cases and deaths by analyzing the state immunization registry and registry of confirmed COVID-19 cases. 
#Post-vaccination cases, also referred to as vaccine breakthrough cases, are individuals who have a positive SARS-Cov-2 molecular test at least 14 days after they have completed their full one-dose or two-dose vaccination series. 
#Tracking cases of COVID-19 that occur after vaccination is important for monitoring the impact of immunization campaigns.
#While COVID-19 vaccines are safe and effective, some cases are still expected in persons who have been vaccinated, as no vaccine is 100% effective.
#Post-vaccination infection data is updated weekly. For comparison, data on cases and deaths among unvaccinated individuals is also included. 
#Partially vaccinated individuals are excluded.
#Note: To account for reporting and processing delays, there is a 7 day lag in provided data (for example, for data through 10/9/2021, only data through 10/2/2021 will be made available). 
#For hospitalizations and deaths, there is an even greater lag in reporting, so more recent data should be used with caution. 
#For display on the public dashboard, there is an additional 7-day lag for hospitalizations (14 days total) and an additional 14-day lag for death data (21 days total).
#Note that this lag is separate from the difference in dates between data processing and updates to the website (in the above example, data through 10/2/2021 would be updated on the website on 10/16/2021).
#open file
covidvax <- read.csv("/Users/sidkurkure/Downloads/COVID-19/covid19postvaxstatewidestats.csv")

covidvax[is.na(covidvax)] = 0


write.csv(covidvax,"/Users/sidkurkure/Downloads/COVID-19\\covid19postvaxstatewidestats.csv", row.names = TRUE)

