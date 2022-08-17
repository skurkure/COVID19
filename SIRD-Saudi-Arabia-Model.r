library(deSolve)
library(tidyverse)

#begin with simple SIR function
#sir function
sir <- function(time, init, parameters) {
  with(as.list(c(init, parameters)), {
    dS <- -beta * S * I / N
    dI <-  beta * S * I / N - gamma * I
    dR <-  gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}

#set the parameters
init       <- c(S = 34217997, I = 3, R = 0)
parameters <- c(N = 34218000, beta = 0.182, gamma = 0.1)
times      <- seq(0, 500, by = 1)

#solve the system of ODE's
out <- ode(y = init, times = times, func = sir, parms = parameters)

#convert to dataframe
out <- as.data.frame(out)
#check values
head(out, 10)
#visualize
SIR_graph <- out %>%
  ggplot(aes(x = times)) +
  labs(y="Individuals",
       x="Time",
       title="COVID-19 Cases in Saudi Arabia, Simulated") +
  geom_line(aes(y= S, color = "Susceptible")) +
  geom_line(aes(y= I, color = "Infectious")) +
  geom_line(aes(y= R, color = "Removed")) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(name = "",
                     breaks = c("Susceptible", "Infectious", "Removed"),
                     values = c("red", "black", "green")) + theme_bw()
SIR_graph
## Alternate Plot Method from Source Code
# matplot(x = times, y = out, type = "l",
#         xlab = "Time", ylab = "Susceptible and Recovered", main = "SIR Model",
#         lwd = 1, lty = 1, bty = "l", col = 2:4)
# 
# ## Add legend
# legend(40, 0.7, c("Susceptible", "Infected", "Recovered"), pch = 1, col = 2:4, bty = "n")



#SIRD model with vax-Algeria model using Saudi Arabia data
sird <- function(time, init, parameters) {
  with(as.list(c(init, parameters)), {
    dS <- -beta * S * I / N + alpha * R - epsilon * S / N
    dI <-  beta * S * I / N  - kappa * rho * I - (1 - kappa) * lambda * I
    dR <- (1-kappa) * lambda * I - alpha * R + epsilon * S / N 
    dD <-  kappa * rho * I
    
    return(list(c(dS, dI, dR, dD)))
  })
}
#variable definitions
#beta: infection rate
#kappa: case fatality rate
#alpha: rate at which vaccine is ineffective and vaccinated individuals lose immunity
#epsilon: vaccination rate
#lambda: recovery time
#rho: time to death
#set up parameters
init_vax <- c(S = 33922997, I = 3, R = 295000, D = 0)
parameters_vax <- c(N = 34218000, alpha = 0.05, beta = 0.858, gamma = 0.182, kappa = 0.014, epsilon = 0.0003, lambda = 0.1, rho = .067)
times_vax      <- seq(0, 500, by = 1)

#solve the system of ODE's
out1 <- ode(y = init_vax, times = times_vax, func = sird, parms = parameters_vax)

#convert to dataframe
out1 <- as.data.frame(out1)
#check values
head(out1, 100)
#visualize
SIRD_graph <- out1 %>%
  ggplot(aes(x = times)) +
  labs(y="Individuals",
       x="Time",
       title="COVID-19 Cases in Saudi Arabia, Simulated") +
  geom_line(aes(y= S, color = "Susceptible")) +
  geom_line(aes(y= I, color = "Infectious")) +
  geom_line(aes(y= R, color = "Removed")) +
  geom_line(aes(y =D, color = "Dead")) +
scale_y_continuous(labels = scales::comma) +
  scale_color_manual(name = "",
                     breaks = c("Susceptible", "Infectious", "Removed", "Dead"),
                     values = c("red", "black", "green", "blue")) + theme_bw()
SIRD_graph

#what it looks like without S line
SIRD_graph_noS <- out1 %>%
  ggplot(aes(x = times)) +
  labs(y="Individuals",
       x="Time",
       title="COVID-19 Cases in Saudi Arabia, Simulated") +
  geom_line(aes(y= I, color = "Infectious")) +
  geom_line(aes(y= R, color = "Removed")) +
  geom_line(aes(y =D, color = "Dead")) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(name = "",
                     breaks = c("Infectious", "Removed", "Dead"),
                     values = c("black", "green", "blue")) + theme_bw()
SIRD_graph_noS

#expanded SIRD model
sird_expanded <- function(time, init, parameters) {
  with(as.list(c(init, parameters)), {
    dS_vax <- -beta_vax * S_vax * I_vax / N + epsilon * sigma * S_unvax / N + alpha * R
    dI_vax <- beta_vax * S_vax * I_vax / N - kappa_vax * rho_vax * I_vax - (1 - kappa_vax) * lambda_vax * R
    dS_unvax <- -beta_unvax * S_unvax * I_unvax / N - epsilon * (1 - sigma) * S_unvax - epsilon * sigma * S_unvax
    dI_unvax <- beta_unvax * S_unvax * I_unvax / N - kappa_unvax * rho_unvax * I_unvax - (1 - kappa_unvax) * lambda_unvax * R
    dR <- (1 - kappa_vax) * lambda_vax * I_vax + (1 - kappa_unvax) * lambda_unvax * I_unvax + epsilon * (1 - sigma) * S_unvax / N - alpha * R
    dD <- kappa_vax * rho_vax * I_vax + kappa_unvax * rho_unvax * I_unvax
    return(list(c(dS_vax, dS_unvax, dI_vax, dI_unvax, dR, dD)))
  })
}
#variable definitions
#beta_vax: infection rate of vaccinated individuals
#beta_unvax: infection rate of unvaccinated individuals
#kappa_vax: case fatality rate of vaccinated individuals
#kappa_unvax: case fatality rate of unvaccinated individuals
#alpha: rate at which vaccine is ineffective and removed individuals lose immunity
#epsilon: vaccination rate
#lambda_vax: recovery time of vaccinated individuals
#lambda_unvax: recovery time of unvaccinated individuals
#rho_vax: time to death for vaccinated individuals
#rho_unvax: time to death for unvaccinated individuals
#sigma: vaccine inefficacy rate
#set up parameters
init_expanded <- c(S_vax = 0, S_unvax = 33922997, I_vax = 0, I_unvax = 3, R = 295000, D = 0)
parameters_expanded <- c(N = 34218000, alpha = 0.05, beta_vax = 0.343, beta_unvax = 0.858, kappa_vax = 0.011, kappa_unvax = 0.017, epsilon = 0.0003, sigma = 0.05, lambda_vax = 0.14, rho_vax = 0.05, lambda_unvax = 0.071, rho_unvax = 0.2)
times <- seq(0, 500, by = 1)

#solve the system of ODE's
out2 <- ode(y = init_expanded, times = times, func = sird_expanded, parms = parameters_expanded)

#convert to dataframe
out2 <- as.data.frame(out2)
#check values
head(out2, 100)
#visualize
SIRD_expanded_graph <- out2 %>%
  ggplot(aes(x = times)) +
  labs(y="Individuals",
       x="Time",
       title="COVID-19 Cases in Saudi Arabia, Simulated") +
  geom_line(aes(y= S_vax, color = "Vaccinated and Susceptible")) +
  geom_line(aes(y= S_unvax, color = "Unvaccinated and Susceptible")) +
  geom_line(aes(y= I_vax, color = "Vaccinated and Infectious")) +
  geom_line(aes(y= I_unvax, color = "Unvaccinated and Infectious")) +
  geom_line(aes(y= R, color = "Recovered")) +
  geom_line(aes(y =D, color = "Dead")) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(name = "",
                     breaks = c("Vaccinated and Susceptible", "Unvaccinated and Susceptible", "Vaccinated and Infectious", "Unvaccinated and Infectious", "Removed", "Dead"),
                     values = c("red","maroon", "black", "grey", "green", "blue")) + theme_bw()
SIRD_expanded_graph
