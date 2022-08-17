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
library(data.table)
#expanded model-add time element to unvaccinated and vaccinated transitions to Recovered or Dead compartments
pmodel <- function(time, init, parameters) {
  with(as.list(c(init, parameters)), {
    dSv <- -beta1 * Iv * Sv / N + alpha * R + epsilon * sigma * Su
    dSu <- -beta2 * Iu * Su / N - epsilon * (1-sigma) * Su - epsilon * sigma * Su
    dIv <-  beta1 * Iv * Sv / N - (1-kappa1) * lambda1 * Iv - kappa1 * rho1 * Iv
    dIu <-  beta2 * Iu * Su / N - (1-kappa2) * lambda2 * Iu - kappa2 * rho2 * Iu
    dR <- (1-kappa1) * lambda1 * Iv + (1-kappa2) * lambda2 * Iu - alpha * R + epsilon * (1-sigma) * Su
    dD <- kappa1 * rho1 * Iv + kappa2 * rho2 * Iu
    
    return(list(c(dSv, dSu, dIv, dIu, dR, dD)))
  })
}


covidvax <- read.csv("~/Downloads/COVID19/covid19postvaxstatewidestats.csv")
head(covidvax$date, 50)
#add total cases column
covidvax <- covidvax %>%
  mutate(total_cases = unvaccinated_cases + vaccinated_cases)
covidvax <- covidvax %>%
  mutate(Submission_Date = as.Date(date))
#population 
N <- covidvax$population[335] #total population on Jan 1, 2022
N_vax <- covidvax$population_vaccinated
N_unvax <- covidvax$population_unvaccinated
#exploratory data analysis
daily_cases_plot <- covidvax %>%
  ggplot(aes(x=Submission_Date,y=total_cases)) + geom_point()  +  
  labs(y="Daily Cases", x = "Submission Date",
       title="COVID-19 Observed Daily Incidence, State of California",
       subtitle = "Beginning February 1, 2021")
daily_cases_plot

#cumulative data and plot
covidvax <- covidvax %>%
  mutate("cum_cases" = cumsum(total_cases))
cum_tot_cases_plot <- covidvax %>%
  ggplot(aes(x=Submission_Date,y=cum_cases)) + geom_point()  +  
  labs(y="Cumulative Cases", x = "Submission Date",
       title="COVID-19 Observed Cumulative Incidence, State of California",
       subtitle = "Beginning February 1, 2021")
covidvax <- covidvax %>%
  mutate(total_deaths = unvaccinated_deaths + vaccinated_deaths)
cum_tot_cases_plot

covidvax <- covidvax %>%
  mutate(population = population_unvaccinated + population_vaccinated)
covidvax <- covidvax %>%
  mutate(cum_unvax_deaths = cumsum(unvaccinated_deaths))
covidvax <- covidvax %>%
  mutate(cum_vax_deaths = cumsum(vaccinated_deaths))
covidvax <- covidvax %>%
  mutate(cum_unvax_deaths = cumsum(unvaccinated_deaths))

#vax data
covidvax <- covidvax %>%
  mutate(cum_vax = cumsum(vaccinated_cases))
cum_vax_cases_plot <- covidvax %>%
  ggplot(aes(x=Submission_Date,y=cum_vax)) + geom_point()  +  
  labs(y="Cumulative Cases", x = "Submission Date",
       title="COVID-19 Vaccinated Observed Cumulative Incidence, State of California",
       subtitle = "Beginning February 1, 2021")
cum_vax_cases_plot

#unvax data
head(covidvax$unvaccinated_cases,100)
covidvax <- covidvax %>%
  mutate(cum_unvax = cumsum(unvaccinated_cases))
cum_unvax_cases_plot <- covidvax %>%
  ggplot(aes(x=Submission_Date,y=cum_unvax)) + geom_point()  +  
  labs(y="Cumulative Cases", x = "Submission Date",
       title="COVID-19 Unvaccinated Cumulative Incidence, State of California",
       subtitle = "Beginning February 1, 2021")
cum_unvax_cases_plot

#I compartments
which(colnames(covidvax) == "cum_vax")
#covidvax$I_vax_active <- c(rep(NA, length.out = 18),covidvax[19:nrow(covidvax), 25] - covidvax[1:(nrow(covidvax)-1), 25])

setDT(covidvax)[,I_vax_active:=cum_vax-shift(cum_vax,18,type="lag")]
covidvax$I_vax_active
setDT(covidvax)[,I_unvax_active:=cum_unvax-shift(cum_vax,18,type="lag")]
covidvax$I_unvax_active

#R compartment
covidvax <- covidvax %>%
  mutate(R_vax = cum_vax - cum_vax_deaths) %>%
  mutate(R_vax = lag(R_vax, 18)) %>%
  na.omit()
covidvax$R_vax
covidvax <- covidvax %>%
  mutate(R_unvax = cum_unvax - cum_unvax_deaths) %>%
  mutate(R_unvax = lag(R_unvax, 18)) %>%
  na.omit()
covidvax$R_unvax
covidvax <- covidvax %>%
  mutate(R_total = R_vax + R_unvax)
covidvax$R_total

omicron <- covidvax %>%
  filter(Submission_Date >="2022-01-01" & Submission_Date <= "2022-2-28")
omicron_vax_death <- omicron %>%
  pull(cum_vax_deaths)
omicron <- omicron %>%
  mutate(time=1:59)

omicron$date
Day <- 1:(length(omicron))
which(covidvax$Submission_Date == "2022-01-01") #335

#figuring out initial conditions
omicron$I_vax_active[1] #285276 for Iv
omicron$I_unvax_active[1] #928838 for Iu
omicron$R_total[1] #1329269 for R
omicron$cum_vax_deaths[1] + omicron$cum_unvax_deaths[1] #17556 for D
omicron$population_vaccinated[1] - omicron$cum_vax[1] #16356633 for Sv
omicron$population_unvaccinated[1] - omicron$cum_unvax[1] #4412601 for Su


#variable definitions
#beta1: infection rate from susceptible-vaccinated (Sv) to infected-vaccinated (Iv)
#beta2: infection rate from susceptible-unvaccinated (Su) to infected-unvaccinated (Iu)
#kappa1: case fatality rate for infected-vaccinated
#kappa2: case fatality rate for infected-unvaccinated
#rho1: inverse of time to death-vaccinated
#rho2: inverse of time to death-unvaccinated
#epsilon: vaccination rate 
#sigma: vaccine inefficacy rate
#set up parameters
init <- c(Sv = 16356633, Su = 4412601, Iv = 285276, Iu = 928838, R = 1329269, D = 17556)
parameters <- c(N = 34218000, beta1 = 0.355, beta2 = 0.0380, kappa1 = 0.0183, kappa2 = 0.0118, epsilon = 9.904923e-06, sigma = 0.05, rho1 = 0.0333, rho2 = 0.0144, alpha = 1/183, lambda1=0.134, lambda2=0.00380)
#beta1 = 0.355, beta2 = 0.130
times      <- seq(0, 500, by = 1)
#solve the system of ODE's
out <- ode(y = init, times = times, func = pmodel, parms = parameters)
#convert to dataframe
out <- as.data.frame(out)
#check values
head(out, 100)
#visualize
pmodel_graph_july <- out %>%
  ggplot(aes(x = times)) +
  labs(y="Individuals",
       x="Time",
       title="Modified SIRD Vaccination Model, Simulated") +
  geom_line(aes(y= Sv, color = "Susceptible-vaccinated")) +
  geom_line(aes(y= Su, color = "Susceptible-unvaccinated")) +
  geom_line(aes(y= Iv, color = "Infectious-vaccinated")) +
  geom_line(aes(y= Iu, color = "Infectious-unvaccinated")) +
  geom_line(aes(y= R, color = "Recovered")) +
  geom_line(aes(y= D, color = "Dead")) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(name = "",
                     breaks = c("Susceptible-vaccinated", "Susceptible-unvaccinated", "Infectious-vaccinated", "Infectious-unvaccinated", "Recovered", "Dead"),
                     values = c("red", "orange", "black", "brown", "green", "blue")) + theme_bw()
pmodel_graph_july

#see if July data matches

times_31 <- seq(1, 31, by = 1)
#solve the system of ODE's
out_2 <- ode(y = init, times = times_31, func = pmodel, parms = parameters)
#convert to dataframe
out_2 <- as.data.frame(out_2)
pmodel_graph_july_with_data <- ggplot(out_2, aes(x=times_31)) +
  labs(y="Individuals",
       x="Time",
       title="Modified SIRD Vaccination Model, Simulated") +
  geom_line(aes(y= Sv, color = "Susceptible-vaccinated")) +
  geom_line(aes(y= Su, color = "Susceptible-unvaccinated")) +
  geom_line(aes(y= Iv, color = "Infectious-vaccinated")) +
  geom_line(aes(y= Iu, color = "Infectious-unvaccinated")) +
  geom_line(aes(y= R, color = "Recovered")) +
  geom_line(aes(y= D, color = "Dead")) +
  geom_point(data = july, aes(y = population_vaccinated - cum_vax, color = "Susceptible-vaccinated")) +
  geom_point(data = july, aes(y = population_unvaccinated - cum_unvax, color = "Susceptible-unvaccinated")) +
  geom_point(data = july, aes(y = I_vax_active, color = "Infectious-vaccinated")) +
  geom_point(data = july, aes(y = I_unvax_active, color = "Infectious-unvaccinated")) +
  geom_point(data = july, aes(y = R_total, color = "Recovered")) +
  geom_point(data = july, aes(y = cum_vax_deaths + cum_unvax_deaths, color = "Dead")) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(name = "",
                     breaks = c("Susceptible-vaccinated", "Susceptible-unvaccinated", "Infectious-vaccinated", "Infectious-unvaccinated", "Recovered", "Dead"),
                     values = c("red", "orange", "black", "brown", "green", "blue")) + theme_bw()
pmodel_graph_july_with_data
#lines are simulated, points are actual values

pmodel_graph_july_with_data_no_S <- ggplot(out_2, aes(x=times_31)) +
  labs(y="Individuals",
        x="Time",
        title="Modified SIRD Vaccination Model, Simulated") +
  geom_line(aes(y= Iv, color = "Infectious-vaccinated")) +
  geom_line(aes(y= Iu, color = "Infectious-unvaccinated")) +
  geom_line(aes(y= R, color = "Recovered")) +
  geom_line(aes(y= D, color = "Dead")) +
  geom_point(data = july, aes(y = I_vax_active, color = "Infectious-vaccinated")) +
  geom_point(data = july, aes(y = I_unvax_active, color = "Infectious-unvaccinated")) +
  geom_point(data = july, aes(y = R_total, color = "Recovered")) +
  geom_point(data = july, aes(y = cum_vax_deaths + cum_unvax_deaths, color = "Dead")) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(name = "",
                      breaks = c("Susceptible-vaccinated", "Susceptible-unvaccinated", "Infectious-vaccinated", "Infectious-unvaccinated", "Recovered", "Dead"),
                      values = c("red", "orange", "black", "brown", "green", "blue")) + theme_bw()
pmodel_graph_july_with_data_no_S


#so now the values fit pretty well, but nothing really matches Iu very well. I'm going to see how the values compare if I use data values from July to August
july_aug <- covidvax %>%
  filter(Submission_Date >="2021-07-01" & Submission_Date <= "2021-8-31")
july_aug <- july_aug %>%
  mutate(time=1:62)
#variable definitions
#beta1: infection rate from susceptible-vaccinated (Sv) to infected-vaccinated (Iv)
#beta2: infection rate from susceptible-unvaccinated (Su) to infected-unvaccinated (Iu)
#kappa1: case fatality rate for infected-vaccinated
#kappa2: case fatality rate for infected-unvaccinated
#rho1: inverse of time to death-vaccinated
#rho2: inverse of time to death-unvaccinated
#epsilon: vaccination rate 
#sigma: vaccine inefficacy rate
#set up parameters
#N = 22632532 on Jan 1, 2022 ???
init <- c(Sv = 20291966, Su = 11012058, Iv = 3476, Iu = 315837, R = 303812, D = 9523)
parameters_ja <- c(N = 34218000, beta1 = 0.355, beta2 = 0.0450, kappa1 = 0.0183, kappa2 = 0.0118, epsilon = 9.904923e-06, sigma = 0.05, rho1 = 0.0333, rho2 = 0.0144, alpha = 1/183, lambda1=0.134, lambda2=0.00380)
times      <- seq(0, 500, by = 1)
#solve the system of ODE's
out_ja <- ode(y = init, times = times, func = pmodel, parms = parameters_ja)
#convert to dataframe
out_ja <- as.data.frame(out_ja)
#check values
head(out_ja, 100)
#visualize
pmodel_graph_july_aug <- out_ja %>%
  ggplot(aes(x = times)) +
  labs(y="Individuals",
       x="Time",
       title="Modified SIRD Vaccination Model, Simulated") +
  geom_line(aes(y= Sv, color = "Susceptible-vaccinated")) +
  geom_line(aes(y= Su, color = "Susceptible-unvaccinated")) +
  geom_line(aes(y= Iv, color = "Infectious-vaccinated")) +
  geom_line(aes(y= Iu, color = "Infectious-unvaccinated")) +
  geom_line(aes(y= R, color = "Recovered")) +
  geom_line(aes(y= D, color = "Dead")) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(name = "",
                     breaks = c("Susceptible-vaccinated", "Susceptible-unvaccinated", "Infectious-vaccinated", "Infectious-unvaccinated", "Recovered", "Dead"),
                     values = c("red", "orange", "black", "brown", "green", "blue")) + theme_bw()
pmodel_graph_july_aug

#see if July and August data matches

times_62 <- seq(1, 62, by = 1)
#solve the system of ODE's
out_ja_2 <- ode(y = init, times = times_62, func = pmodel, parms = parameters_ja)
#convert to dataframe
out_ja_2 <- as.data.frame(out_ja_2)
pmodel_graph_july_aug_with_data <- ggplot(out_ja_2, aes(x=times_62)) +
  labs(y="Individuals",
       x="Time",
       title="Modified SIRD Vaccination Model, Simulated") +
  geom_line(aes(y= Sv, color = "Susceptible-vaccinated")) +
  geom_line(aes(y= Su, color = "Susceptible-unvaccinated")) +
  geom_line(aes(y= Iv, color = "Infectious-vaccinated")) +
  geom_line(aes(y= Iu, color = "Infectious-unvaccinated")) +
  geom_line(aes(y= R, color = "Recovered")) +
  geom_line(aes(y= D, color = "Dead")) +
  geom_point(data = july_aug, aes(y = population_vaccinated - cum_vax, color = "Susceptible-vaccinated")) +
  geom_point(data = july_aug, aes(y = population_unvaccinated - cum_unvax, color = "Susceptible-unvaccinated")) +
  geom_point(data = july_aug, aes(y = I_vax_active, color = "Infectious-vaccinated")) +
  geom_point(data = july_aug, aes(y = I_unvax_active, color = "Infectious-unvaccinated")) +
  geom_point(data = july_aug, aes(y = R_total, color = "Recovered")) +
  geom_point(data = july_aug, aes(y = cum_vax_deaths + cum_unvax_deaths, color = "Dead")) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(name = "",
                     breaks = c("Susceptible-vaccinated", "Susceptible-unvaccinated", "Infectious-vaccinated", "Infectious-unvaccinated", "Recovered", "Dead"),
                     values = c("red", "orange", "black", "brown", "green", "blue")) + theme_bw()
pmodel_graph_july_aug_with_data
#lines are simulated, points are actual values

pmodel_graph_july_aug_with_data_no_S <- ggplot(out_ja_2, aes(x=times_62)) +
  labs(y="Individuals",
       x="Time",
       title="Modified SIRD Vaccination Model, Simulated") +
  geom_line(aes(y= Iv, color = "Infectious-vaccinated")) +
  geom_line(aes(y= Iu, color = "Infectious-unvaccinated")) +
  geom_line(aes(y= R, color = "Recovered")) +
  geom_line(aes(y= D, color = "Dead")) +
  geom_point(data = july_aug, aes(y = I_vax_active, color = "Infectious-vaccinated")) +
  geom_point(data = july_aug, aes(y = I_unvax_active, color = "Infectious-unvaccinated")) +
  geom_point(data = july_aug, aes(y = R_total, color = "Recovered")) +
  geom_point(data = july_aug, aes(y = cum_vax_deaths + cum_unvax_deaths, color = "Dead")) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(name = "",
                     breaks = c("Susceptible-vaccinated", "Susceptible-unvaccinated", "Infectious-vaccinated", "Infectious-unvaccinated", "Recovered", "Dead"),
                     values = c("red", "orange", "black", "brown", "green", "blue")) + theme_bw()
pmodel_graph_july_aug_with_data_no_S
#slightly different rates if using the July-Aug model. Some rates model the data better than others.
#For July data:
#parameters <- c(N = 34218000, beta1 = 0.355, beta2 = 0.130, kappa1 = 0.0183, kappa2 = 0.0118, epsilon = 9.904923e-06, sigma = 0.05, rho1 = 0.0333, rho2 = 0.0144, alpha = 1/183, lambda1=0.134, lambda2=0.00380)
#For July-August data:
#parameters_ja <- c(N = 34218000, beta1 = 0.355, beta2 = 0.0450, kappa1 = 0.0183, kappa2 = 0.0118, epsilon = 9.904923e-06, sigma = 0.05, rho1 = 0.0333, rho2 = 0.0144, alpha = 1/183, lambda1=0.134, lambda2=0.00380)
#the only difference is in beta2
#now to play with the vaccination rate(epsilon) and see what happens
#I think I will just use July data rates.


# #Epsilon = 5e-05
# #set up parameters
# init <- c(Sv = 20291966, Su = 11012058, Iv = 3476, Iu = 315837, R = 303812, D = 9523)
# parameters3 <- c(N = 34218000, beta1 = 0.355, beta2 = 0.0380, kappa1 = 0.0183, kappa2 = 0.0118, epsilon = 5e-05, sigma = 0.05, rho1 = 0.0333, rho2 = 0.0144, alpha = 1/183, lambda1=0.134, lambda2=0.00380)
# #beta1 = 0.355, beta2 = 0.130
# times <- seq(0, 500, by = 1)
# #solve the system of ODE's
# out_3 <- ode(y = init, times = times, func = pmodel, parms = parameters3)
# #convert to dataframe
# out_3 <- as.data.frame(out_3)
# #check values
# head(out_3, 100)
# #visualize
# pmodel_graph_july_epsilon_3 <- out_3 %>%
#   ggplot(aes(x = times)) +
#   labs(y="Individuals",
#        x="Time",
#        title="Modified SIRD Vaccination Model, Simulated") +
#   geom_line(aes(y= Sv, color = "Susceptible-vaccinated")) +
#   geom_line(aes(y= Su, color = "Susceptible-unvaccinated")) +
#   geom_line(aes(y= Iv, color = "Infectious-vaccinated")) +
#   geom_line(aes(y= Iu, color = "Infectious-unvaccinated")) +
#   geom_line(aes(y= R, color = "Recovered")) +
#   geom_line(aes(y= D, color = "Dead")) +
#   scale_y_continuous(labels = scales::comma) +
#   scale_color_manual(name = "",
#                      breaks = c("Susceptible-vaccinated", "Susceptible-unvaccinated", "Infectious-vaccinated", "Infectious-unvaccinated", "Recovered", "Dead"),
#                      values = c("red", "orange", "black", "brown", "green", "blue")) + theme_bw()
# pmodel_graph_july_epsilon_3

#Rate 2: Epsilon = 1e-04
#set up parameters
init <- c(Sv = 20291966, Su = 11012058, Iv = 3476, Iu = 315837, R = 303812, D = 9523)
parameters4 <- c(N = 34218000, beta1 = 0.355, beta2 = 0.0380, kappa1 = 0.0183, kappa2 = 0.0118, epsilon = 1e-04, sigma = 0.05, rho1 = 0.0333, rho2 = 0.0144, alpha = 1/183, lambda1=0.134, lambda2=0.00380)
#beta1 = 0.355, beta2 = 0.130
times      <- seq(0, 500, by = 1)
#solve the system of ODE's
out_4 <- ode(y = init, times = times, func = pmodel, parms = parameters4)
#convert to dataframe
out_4 <- as.data.frame(out_4)
#check values
head(out_4, 100)
#visualize
pmodel_graph_july_epsilon_4 <- out_4 %>%
  ggplot(aes(x = times)) +
  labs(y="Individuals",
       x="Time",
       title="Modified SIRD Vaccination Model, Simulated") +
  geom_line(aes(y= Sv, color = "Susceptible-vaccinated")) +
  geom_line(aes(y= Su, color = "Susceptible-unvaccinated")) +
  geom_line(aes(y= Iv, color = "Infectious-vaccinated")) +
  geom_line(aes(y= Iu, color = "Infectious-unvaccinated")) +
  geom_line(aes(y= R, color = "Recovered")) +
  geom_line(aes(y= D, color = "Dead")) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(name = "",
                     breaks = c("Susceptible-vaccinated", "Susceptible-unvaccinated", "Infectious-vaccinated", "Infectious-unvaccinated", "Recovered", "Dead"),
                     values = c("red", "orange", "black", "brown", "green", "blue")) + theme_bw()
pmodel_graph_july_epsilon_4

#Epsilon = 1e-03
#set up parameters
init <- c(Sv = 20291966, Su = 11012058, Iv = 3476, Iu = 315837, R = 303812, D = 9523)
parameters4 <- c(N = 34218000, beta1 = 0.355, beta2 = 0.0380, kappa1 = 0.0183, kappa2 = 0.0118, epsilon = 1e-03, sigma = 0.05, rho1 = 0.0333, rho2 = 0.0144, alpha = 1/183, lambda1=0.134, lambda2=0.00380)
#beta1 = 0.355, beta2 = 0.130
times      <- seq(0, 500, by = 1)
#solve the system of ODE's
out_4 <- ode(y = init, times = times, func = pmodel, parms = parameters4)
#convert to dataframe
out_4 <- as.data.frame(out_4)
#check values
head(out_4, 100)
#visualize
pmodel_graph_july_epsilon_10x <- out_4 %>%
  ggplot(aes(x = times)) +
  labs(y="Individuals",
       x="Time",
       title="Modified SIRD Vaccination Model, Simulated") +
  geom_line(aes(y= Sv, color = "Susceptible-vaccinated")) +
  geom_line(aes(y= Su, color = "Susceptible-unvaccinated")) +
  geom_line(aes(y= Iv, color = "Infectious-vaccinated")) +
  geom_line(aes(y= Iu, color = "Infectious-unvaccinated")) +
  geom_line(aes(y= R, color = "Recovered")) +
  geom_line(aes(y= D, color = "Dead")) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(name = "",
                     breaks = c("Susceptible-vaccinated", "Susceptible-unvaccinated", "Infectious-vaccinated", "Infectious-unvaccinated", "Recovered", "Dead"),
                     values = c("red", "orange", "black", "brown", "green", "blue")) + theme_bw()
pmodel_graph_july_epsilon_10x

#Rate: epsilon = 0.01
#set up parameters
init <- c(Sv = 20291966, Su = 11012058, Iv = 3476, Iu = 315837, R = 303812, D = 9523)
parameters5 <- c(N = 34218000, beta1 = 0.355, beta2 = 0.0380, kappa1 = 0.0183, kappa2 = 0.0118, epsilon = 0.01, sigma = 0.05, rho1 = 0.0333, rho2 = 0.0144, alpha = 1/183, lambda1=0.134, lambda2=0.00380)
#beta1 = 0.355, beta2 = 0.130
times      <- seq(0, 500, by = 1)
#solve the system of ODE's
out_5 <- ode(y = init, times = times, func = pmodel, parms = parameters5)
#convert to dataframe
out_5 <- as.data.frame(out_5)
#check values
head(out_5, 100)
#visualize
pmodel_graph_july_epsilon_5 <- out_5 %>%
  ggplot(aes(x = times)) +
  labs(y="Individuals",
       x="Time",
       title="Modified SIRD Vaccination Model, Simulated") +
  geom_line(aes(y= Sv, color = "Susceptible-vaccinated")) +
  geom_line(aes(y= Su, color = "Susceptible-unvaccinated")) +
  geom_line(aes(y= Iv, color = "Infectious-vaccinated")) +
  geom_line(aes(y= Iu, color = "Infectious-unvaccinated")) +
  geom_line(aes(y= R, color = "Recovered")) +
  geom_line(aes(y= D, color = "Dead")) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(name = "",
                     breaks = c("Susceptible-vaccinated", "Susceptible-unvaccinated", "Infectious-vaccinated", "Infectious-unvaccinated", "Recovered", "Dead"),
                     values = c("red", "orange", "black", "brown", "green", "blue")) + theme_bw()
pmodel_graph_july_epsilon_5

#Rate 6: Epsilon = 0.1
#set up parameters
init <- c(Sv = 20291966, Su = 11012058, Iv = 3476, Iu = 315837, R = 303812, D = 9523)
parameters6 <- c(N = 34218000, beta1 = 0.355, beta2 = 0.0380, kappa1 = 0.0183, kappa2 = 0.0118, epsilon = 0.1, sigma = 0.05, rho1 = 0.0333, rho2 = 0.0144, alpha = 1/183, lambda1=0.134, lambda2=0.00380)
#beta1 = 0.355, beta2 = 0.130
times      <- seq(0, 500, by = 1)
#solve the system of ODE's
out_6 <- ode(y = init, times = times, func = pmodel, parms = parameters6)
#convert to dataframe
out_6 <- as.data.frame(out_6)
#check values
head(out_6, 100)
#visualize
pmodel_graph_july_epsilon_6 <- out_6 %>%
  ggplot(aes(x = times)) +
  labs(y="Individuals",
       x="Time",
       title="Modified SIRD Vaccination Model, Simulated") +
  geom_line(aes(y= Sv, color = "Susceptible-vaccinated")) +
  geom_line(aes(y= Su, color = "Susceptible-unvaccinated")) +
  geom_line(aes(y= Iv, color = "Infectious-vaccinated")) +
  geom_line(aes(y= Iu, color = "Infectious-unvaccinated")) +
  geom_line(aes(y= R, color = "Recovered")) +
  geom_line(aes(y= D, color = "Dead")) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(name = "",
                     breaks = c("Susceptible-vaccinated", "Susceptible-unvaccinated", "Infectious-vaccinated", "Infectious-unvaccinated", "Recovered", "Dead"),
                     values = c("red", "orange", "black", "brown", "green", "blue")) + theme_bw()
pmodel_graph_july_epsilon_6


#Make Su, Sv, Iu, Iv, R, D graphs now
#Set dates for each of the dataframes so far
out <- out %>%
  mutate(date = seq(as.Date("2021/07/01"), by = "day", length.out = 501))
out_3 <- out_3 %>%
  mutate(date = seq(as.Date("2021/07/01"), by = "day", length.out = 501))
out_4 <- out_4 %>%
  mutate(date = seq(as.Date("2021/07/01"), by = "day", length.out = 501))
out_5 <- out_5 %>%
  mutate(date = seq(as.Date("2021/07/01"), by = "day", length.out = 501))
out_6 <- out_6 %>%
  mutate(date = seq(as.Date("2021/07/01"), by = "day", length.out = 501))


#Su graph
Su_graph <- 
  ggplot() +
  labs(y = "Individuals", x = "Date", title = "Susceptible and Unvaccinated Compartment") +
  geom_line( data = out, aes(x = date, y = Su), color = "black") +
  geom_line( data = out_4, aes(x = date, y = Su), color = "blue") +
  geom_line( data = out_5, aes(x = date, y = Su), color = "green") +
  geom_line( data = out_6, aes(x = date, y = Su), color = "red") +
  theme_bw()
Su_graph

#Sv graph
Sv_graph <- 
  ggplot() +
  labs(y = "Individuals", x = "Date", title = "Susceptible and Vaccinated Compartment") +
  geom_line( data = out, aes(x = date, y = Sv), color = "black") +
  geom_line( data = out_4, aes(x = date, y = Sv), color = "blue") +
  geom_line( data = out_5, aes(x = date, y = Sv), color = "green") +
  geom_line( data = out_6, aes(x = date, y = Sv), color = "red") +
  theme_bw()
Sv_graph

#Iu graph
Iu_graph <- 
  ggplot() +
  labs(y = "Individuals", x = "Date", title = "Infectious and Unvaccinated Compartment") +
  geom_line( data = out, aes(x = date, y = Iu), color = "black") +
  geom_line( data = out_4, aes(x = date, y = Iu), color = "blue") +
  geom_line( data = out_5, aes(x = date, y = Iu), color = "green") +
  geom_line( data = out_6, aes(x = date, y = Iu), color = "red") +
  theme_bw()
Iu_graph

#Iv graph
Iv_graph <- 
  ggplot() +
  labs(y = "Individuals", x = "Date", title = "Susceptible and Unvaccinated Compartments") +
  geom_line( data = out, aes(x = date, y = Iv), color = "black") +
  geom_line( data = out_4, aes(x = date, y = Iv), color = "blue") +
  geom_line( data = out_5, aes(x = date, y = Iv), color = "green") +
  geom_line( data = out_6, aes(x = date, y = Iv), color = "red") +
  theme_bw()
Iv_graph

#R graph
R_graph <- 
  ggplot() +
  labs(y = "Individuals", x = "Date", title = "Susceptible and Unvaccinated Compartments") +
  geom_line( data = out, aes(x = date, y = R), color = "black") +
  geom_line( data = out_4, aes(x = date, y = R), color = "blue") +
  geom_line( data = out_5, aes(x = date, y = R), color = "green") +
  geom_line( data = out_6, aes(x = date, y = R), color = "red") +
  theme_bw()
R_graph

#D graph
D_graph <- 
  ggplot() +
  labs(y = "Individuals", x = "Date", title = "Susceptible and Unvaccinated Compartments") +
  geom_line( data = out, aes(x = date, y = D), color = "black") +
  geom_line( data = out_4, aes(x = date, y = D), color = "blue") +
  geom_line( data = out_5, aes(x = date, y = D), color = "green") +
  geom_line( data = out_6, aes(x = date, y = D), color = "red") +
  theme_bw()
D_graph