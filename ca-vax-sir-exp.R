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

#Puglisi model
pmodel <- function(time, init, parameters) {
  with(as.list(c(init, parameters)), {
    dSv <- -beta1 * Iv * Sv / N + zeta * Su / N
    dSu <- -beta2 * Iu * Su / N - zeta * Su / N
    dIv <-  beta1 * Iv * Sv / N - gamma1 * Iv - kappa1 * Iv
    dIu <-  beta2 * Iu * Su / N - gamma2 * Iu - kappa2 * Iu
    dR <- gamma1 * Iv + gamma2 * Iu 
    dD <- kappa1 * Iv + kappa2 * Iu
    
    return(list(c(dSv, dSu, dIv, dIu, dR, dD)))
  })
}
#variable definitions
#beta1: infection rate from susceptible-vaccinated (Sv) to infected-vaccinated (Iv)
#beta2: infection rate from susceptible-unvaccinated (Su) to infected-unvaccinated (Iu)
#kappa1: case fatality rate for infected-vaccinated
#kappa2: case fatality rate for infected-unvaccinated
#gamma1: death rate for infected-vaccinated
#gamma2: death rate for infected-unvaccinated
#zeta: vaccination rate
#set up parameters
init <- c(Sv = 33922997, Su = 4877123, Iv = 3, Iu = 523, R = 295000, D = 0)
parameters <- c(N = 34218000, beta1 = 0.858, beta2 = 0.9, kappa1 = 0.014, kappa2 = 0.03, gamma1 = 0.082, gamma2 = 0.182, zeta = 0.0003, alpha = 0.8419, lambda1=1/10.3,lambda2=1/16.7)
times      <- seq(0, 500, by = 1)
#solve the system of ODE's
out <- ode(y = init, times = times, func = pmodel, parms = parameters)
#convert to dataframe
out <- as.data.frame(out)
#check values
head(out, 100)
#visualize
pmodel_graph <- out %>%
  ggplot(aes(x = times)) +
  labs(y="Individuals",
       x="Time",
       title="Modified SIRD vaccination model, Simulated") +
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
pmodel_graph
#what it looks like without S line
pmodel_graph_noS <- out %>%
  ggplot(aes(x = times)) +
  labs(y="Individuals",
       x="Time",
       title="Modified SIRD vaccination model, Simulated") +
  geom_line(aes(y= Iv, color = "Infectious-vaccinated")) +
  geom_line(aes(y= Iu, color = "Infectious-unvaccinated")) +
  geom_line(aes(y= R, color = "Recovered")) +
  geom_line(aes(y= D, color = "Dead")) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(name = "",
                     breaks = c("Infectious-vaccinated", "Infectious-unvaccinated", "Recovered", "Dead"),
                     values = c("black", "brown", "green", "blue")) + theme_bw()
pmodel_graph_noS
#expanded model-add time element to unvaccinated and vaccinated transitions to Recovered or Dead compartments
pmodel <- function(time, init, parameters) {
  with(as.list(c(init, parameters)), {
    dSv <- -beta1 * Iv * Sv / N + alpha * R + epsilon * sigma * Su / N
    dSu <- -beta2 * Iu * Su / N - epsilon * (1-sigma) * Su / N - epsilon * sigma * Su / N
    dIv <-  beta1 * Iv * Sv / N - (1-kappa1) * lambda1 * Iv - kappa1 * rho1 * Iv
    dIu <-  beta2 * Iu * Su / N - (1-kappa2) * rho1 * Iu - kappa2 * rho2 * Iu
    dR <- (1-kappa1) * lambda1 * Iv + (1-kappa2) * lambda2 * Iu - alpha * R + epsilon * (1-sigma) * Su / N
    dD <- kappa1 * rho1 * Iv + kappa2 * rho2 * Iu
    
    return(list(c(dSv, dSu, dIv, dIu, dR, dD)))
  })
}
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
#hard coded alpha and lambda1,lambda2 taken from https://www.nejm.org/doi/full/10.1056/NEJMoa2107058
init <- c(Sv = 33922997, Su = 4877123, Iv = 3, Iu = 523, R = 295000, D = 0)
parameters <- c(N = 34218000, beta1 = 0.858, beta2 = 0.9, kappa1 = 0.014, kappa2 = 0.03, epsilon = 0.0003, sigma = 0.05, rho1 = 0.05, rho2 = 0.2, alpha = 0.8419, lambda1=1/10.3,lambda2=1/16.7)
times      <- seq(0, 500, by = 1)
#solve the system of ODE's
out <- ode(y = init, times = times, func = pmodel, parms = parameters)
#convert to dataframe
out <- as.data.frame(out)
#check values
head(out, 100)
#visualize
pmodel_graph <- out %>%
  ggplot(aes(x = times)) +
  labs(y="Individuals",
       x="Time",
       title="Modified SIRD vaccination model, Simulated") +
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
pmodel_graph
#what it looks like without S line
pmodel_graph_noS <- out %>%
  ggplot(aes(x = times)) +
  labs(y="Individuals",
       x="Time",
       title="Modified SIRD vaccination model, Simulated") +
  geom_line(aes(y= Iv, color = "Infectious-vaccinated")) +
  geom_line(aes(y= Iu, color = "Infectious-unvaccinated")) +
  geom_line(aes(y= R, color = "Recovered")) +
  geom_line(aes(y= D, color = "Dead")) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(name = "",
                     breaks = c("Infectious-vaccinated", "Infectious-unvaccinated", "Recovered", "Dead"),
                     values = c("black", "brown", "green", "blue")) + theme_bw()
pmodel_graph_noS
