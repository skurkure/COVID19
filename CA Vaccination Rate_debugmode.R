library(tidyverse)
library(ggplot2)
library(EpiEstim)
library(deSolve)

vaccination_rate_data <- read.csv("covid19vaccinesbycounty.csv")
agg_vax <- aggregate(cumulative_fully_vaccinated~administered_date, data = vaccination_rate_data, sum)
agg_vax
july_agg_vax <- agg_vax %>%
  filter(administered_date >= "2021-07-01" & administered_date <= "2021-07-31") %>%
  pull(cumulative_fully_vaccinated)
july_agg_vax
#find rate of vaccination for July
#S unvaxed to R 

SR_unvax <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- - esm * S/N
    dR <- esm * S/N
    list(c(dS, dR))
  })
}
N <- 39350000
days <- 1:length(july_agg_vax)
days
#residual sum of squares function
RSS7 <- function(parameters) {
  names(parameters) <- c("esm")
  out <- ode(y = init7, times = days, func = SR_unvax, parms = parameters)
  fit <- out[ , 3]
  sum(july_agg_vax - fit)^2
}
covidvax <- read.csv("covid19postvaxstatewidestats.csv")
covidvax <- covidvax %>%
  mutate(Submission_Date = as.Date(date))
covidvax <- covidvax %>%
  mutate(cum_unvax_deaths = cumsum(unvaccinated_deaths))
covidvax <- covidvax %>%
  mutate(cum_unvax = cumsum(unvaccinated_cases))
covidvax <- covidvax %>%
  mutate(cum_unvax_recovered = cum_unvax - cum_unvax_deaths)
july <- covidvax %>%
  filter(Submission_Date >="2021-07-01" & Submission_Date <= "2021-7-31")
init7 <- c(S = july$population_unvaccinated[1]-july$cum_unvax[1], R = july$cum_unvax_recovered[1])
july$population_unvaccinated[1]
july$cum_unvax[1]
july$cum_unvax_recovered[1]
#Nelder-Mead”, “BFGS”, “CG”, “L-BFGS-B”, “SANN”, “Brent"
#optimization
Opt7 <- optim(c(0.5), RSS7, method = "Nelder-Mead")
#check for convergence
Opt7$message

Opt_par7 <- setNames(Opt7$par, c("esm"))


#that can't be right-fix this
Opt_par7
#SANN
# esm 
# 33.90357 
#SANN
# esm 
# 496.7213
#CG
# esm 
# 1516.631 
#Nelder-Mead
# esm 
# 38.85 
#ggplot(aes(x=date)) + geom_point(aes(y=cumulative_fully_vaccinated))
date = agg_vax$administered_date
cumulative_fully_vaccinated = agg_vax$cumulative_fully_vaccinated
agg_vax %>% ggplot(aes(x=date)) + geom_point(aes(y=cumulative_fully_vaccinated))


#with including N, got 0.55 for Nelder Mead, 501.8158 for SANN, 20620707841 for CG
#
#