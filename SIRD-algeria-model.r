
#packages.install("deSolve")
#packages.install("tidyverse")
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
       title="COVID-19 Cases in Algeria, Simulated") +
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



#SIRD model with vax-Algeria model using Algeria data
sird <- function(time, init, parameters) {
  with(as.list(c(init, parameters)), {
    dS <- -beta * S * I / N + alpha * R - epsilon * S / N
    dI <-  beta * S * I / N  - kappa * rho * I - (1 - kappa) * lambda * I
    dR <- (1-kappa) * lambda * I - alpha * R + epsilon * S / N 
    #dD <-  kappa * rho * I
    
    #return(list(c(dS, dI, dR, dD)))
    return(list(c(dS, dI, dR)))
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
init <- c(S = 43851044, I = 0, R = 43851044)
params <- c(N = 43851044, alpha = 0.05, beta = 0.0561215, gamma = 0.0455331, kappa = 0.014, epsilon = 0.0003, lambda = 0.1, rho = .067)
times      <- seq(0, 500, by = 1)

#solve the system of ODE's
out1 <- ode(y = init, times = times, func = sird, parms = params)

#convert to dataframe
out1 <- as.data.frame(out1)
#check values
head(out1, 100)
#visualize
SIRD_graph <- out1 %>%
  ggplot(aes(x = times)) +
  labs(y="Individuals",
       x="Time",
       title="COVID-19 Cases in Algeria, Simulated") +
  geom_line(aes(y= S, color = "Susceptible")) +
  geom_line(aes(y= I, color = "Infectious")) +
  geom_line(aes(y= R, color = "Removed")) +
  #geom_line(aes(y =D, color = "Dead")) +
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
       title="COVID-19 Cases in Algeria, Simulated") +
  geom_line(aes(y= I, color = "Infectious")) +
  geom_line(aes(y= R, color = "Removed")) +
  #geom_line(aes(y =D, color = "Dead")) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(name = "",
                     breaks = c("Infectious", "Removed", "Dead"),
                     values = c("black", "green", "blue")) + theme_bw()
SIRD_graph_noS
