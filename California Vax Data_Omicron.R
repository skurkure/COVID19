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

#set directory
#source: https://data.chhs.ca.gov/dataset/covid-19-time-series-metrics-by-county-and-state
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
covidvax <- read.csv("~/Downloads/covid19postvaxstatewidestats_updated.csv")
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
covidvax <- covidvax %>%
  mutate(population = population_unvaccinated + population_vaccinated)

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


#figuring out I and R compartments
#I compartments
which(colnames(covidvax) == "cum_vax") #column 25
#covidvax$I_vax_active <- c(rep(NA, length.out = 18),covidvax[19:nrow(covidvax), 25] - covidvax[1:(nrow(covidvax)-1), 25])

setDT(covidvax)[,I_vax_active:=cum_vax-shift(cum_vax,18,type="lag")]
covidvax$I_vax_active
setDT(covidvax)[,I_unvax_active:=cum_unvax-shift(cum_vax,18,type="lag")]
covidvax$I_unvax_active

#according to https://pubmed.ncbi.nlm.nih.gov/34154563/#:~:text=The%20estimated%20mean%20time%20between,%3B%20heterogeneity%20P%20%3D%200.85).,
#the average time to death is 18 days. i'll take cases from 18 days prior and subtract
#deaths from the day 18 days later.

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


#these are for the cumulative R cases
covidvax <- covidvax %>%
  mutate(cum_vax_recovered = cum_vax - cum_vax_deaths)
covidvax <- covidvax %>%
  mutate(cum_unvax_recovered = cum_unvax - cum_unvax_deaths)
omicron <- covidvax %>%
  filter(Submission_Date >="2022-01-01" & Submission_Date <= "2022-2-28")

#for recovered, i have to subtract the deaths from the cumulative cases (as well as the active cases). i don't think they
#occur on the same day, so I will subtract deaths from a few days after from when cases are reported

#i will plot out cases and deaths to get a general sense of the two trends

omicron_vax_cases_deaths_plot <- omicron %>%
  ggplot(aes(x=date)) + geom_point(aes(y=cum_vax), color = "green") + geom_point(aes(y=cum_vax_deaths), color= "black")
omicron_vax_cases_deaths_plot
omicron %>%
  ggplot(aes(x=date)) + geom_point(aes(y=cum_vax_deaths), color= "black")
#fit the data

#Full Model
sird_expanded <- function(time, init, parameters) {
  with(as.list(c(init, parameters)), {
    dS_vax <- -beta_vax * S_vax * I_vax / N + epsilon * sigma * S_unvax + alpha * R 
    dI_vax <- beta_vax * S_vax * I_vax / N - kappa_vax * rho_vax * I_vax - (1 - kappa_vax) * lambda_vax * I_vax
    dS_unvax <- -beta_unvax * S_unvax * I_unvax / N - epsilon * (1 - sigma) * S_unvax - epsilon * sigma * S_unvax
    dI_unvax <- beta_unvax * S_unvax * I_unvax / N - kappa_unvax * rho_unvax * I_unvax - (1 - kappa_unvax) * lambda_unvax * I_unvax
    dR <- (1 - kappa_vax) * lambda_vax * I_vax + (1 - kappa_unvax) * lambda_unvax * I_unvax + epsilon * (1 - sigma) * S_unvax / N  - alpha * R 
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
#rho_vax: inverse of time to death for vaccinated individuals
#rho_unvax: inverse of time to death for unvaccinated individuals
#sigma: vaccine inefficacy rate


#will probably need to find individual rates one by one, so we won't use this for now


#set arbitrary initial conditions-try July when cases begin spiking
#this is the Delta wave
N <- 22632532
omicron <- covidvax %>%
  filter(Submission_Date >="2022-01-01" & Submission_Date <= "2022-2-28")
omicron_vax_death <- omicron %>%
  pull(cum_vax_deaths)
omicron$date
Day <- 1:(length(omicron))
which(covidvax$Submission_Date == "2022-01-01") #335
omicron_data <- covidvax %>%
  filter(Submission_Date >= "2022-01-01" & Submission_Date < "2022-02-28") %>%
  pull(cum_cases)
omicron_data

#vax pop only
#the infectious person infects one person on average every 2 days and is infectious for 5 days.
#for example, β=1/2 day⁻¹=0.5 day⁻¹ and γ=1/5 day⁻¹=0.2 day⁻¹
#time in days for predictions
sird_starting_date <- "2020-07-01"
tail(covidvax, 1) #Dec 4, 2021, row 307
sird_ending_date <- "2021-12-04"
t <- 151:307
#vax data
omicron_vax_data <- covidvax %>%
  filter(Submission_Date >= "2022-01-01" & Submission_Date <= "2022-02-28") %>%
  pull(cum_vax)
#nonvax data now
omicron_unvax_data <- covidvax %>%
  filter(Submission_Date >= "2022-01-01" & Submission_Date <= "2022-02-28") %>%
  pull(cum_unvax)

#S to I vax
SI_vax <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- - beta * S * I / N
    dI <- beta * S * I / N
    list(c(dS, dI))
  })
}
days <- 1:length(omicron_vax_death)
#residual sum of squares function
RSS_new <- function(parameters) {
  names(parameters) <- c("beta")
  out <- ode(y = init_new, times = days, func = SI_vax, parms = parameters)
  fit <- out[ , 3]
  sum((omicron_vax_data - fit)^2)
}

init_new <- c(S = omicron$population_vaccinated[1]-omicron$cum_vax[1], I = omicron$I_vax_active[1])
#Nelder-Mead”, “BFGS”, “CG”, “L-BFGS-B”, “SANN”, “Brent"
#optimization
Opt_new <- optim(c(0.5), RSS_new, method = "Brent")
#check for convergence
Opt_new$message

Opt_par_new <- setNames(Opt_new$par, c("beta"))
Opt_par_new
#SANN
# beta 
# 0.05151367
#Nelder-Mead
# beta 
# 0.05151367
#Brent
# beta 
# 0.1855103 


#S to I unvax
SI_unvax <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- - beta * S * I / N
    dI <- beta * S * I / N
    list(c(dS, dI))
  })
}
days <- 1:length(omicron_vax_death)
#residual sum of squares function
RSS_new2 <- function(parameters) {
  names(parameters) <- c("beta")
  out <- ode(y = init_new2, times = days, func = SI_unvax, parms = parameters)
  fit <- out[ , 3]
  sum((omicron_unvax_data - fit)^2)
}

init_new2 <- c(S = omicron$population_unvaccinated[1]-omicron$cum_unvax[1], I = omicron$I_unvax_active[1])
#Nelder-Mead”, “BFGS”, “CG”, “L-BFGS-B”, “SANN”, “Brent"
#optimization
Opt_new2 <- optim(c(0.5), RSS_new2, method = "Brent", upper = c(1), lower = c(0))
#check for convergence
Opt_new2$message

Opt_par_new2 <- setNames(Opt_new2$par, c("beta"))
Opt_par_new2
#SANN
# beta 
# 0.0380466 
#Nelder-Mead
# beta 
# 0.03789062 
#Brent
# beta 
# 0.03803438 

#I to D vax
omicron_vax_death <- omicron %>%
  pull(cum_vax_deaths)
ID_vax <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dI <- -kr * I
    dD <- kr * I
    list(c(dI, dD))
  })
}

days <- 1:length(omicron_vax_death)
days
#residual sum of squares function
RSS3 <- function(parameters) {
  names(parameters) <- c("kr")
  out <- ode(y = init3, times = days, func = ID_vax, parms = parameters)
  fit <- out[ , 3]
  sum((omicron_vax_death - fit)^2)
}

init3 <- c(I = omicron$I_vax_active[1], D = omicron$cum_vax_deaths[1])
#Nelder-Mead”, “BFGS”, “CG”, “L-BFGS-B”, “SANN”, “Brent"
#optimization
Opt3 <- optim(c(0.5), RSS3, method = "Brent", upper = c(1), lower = c(0))
#check for convergence
Opt3$message

Opt_par3 <- setNames(Opt3$par, c("kr"))
Opt_par3

#Nelder-Mead
# kr 
# 0.0006103516
#SANN
# kr 
# 0.0006388306 
#Brent
# kr 
# 0.0006077762 


omicron %>%
  ggplot((aes(x=date))) + geom_point(aes(y=cum_vax), color = "green") + geom_point(aes(y=cum_vax_deaths), color = "black")
omicron %>%
  ggplot((aes(x=date))) + geom_point(aes(y=cum_vax_deaths), color = "black")


#I to D unvax
omicron_unvax_death <- omicron %>%
  pull(cum_unvax_deaths)
ID_unvax <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dI <- -kr * I
    dD <- kr * I
    list(c(dI, dD))
  })
}

days <- 1:length(omicron_vax_death)
days
#residual sum of squares function
RSS4 <- function(parameters) {
  names(parameters) <- c("kr")
  out <- ode(y = init4, times = days, func = ID_unvax, parms = parameters)
  fit <- out[ , 3]
  sum((omicron_unvax_death - fit)^2)
}

init4 <- c(I = omicron$I_unvax_active[1], D = omicron$cum_unvax_deaths[1])
#Nelder-Mead”, “BFGS”, “CG”, “L-BFGS-B”, “SANN”, “Brent"
#optimization
Opt4 <- optim(c(0.5), RSS4, method = "Brent", upper = c(1), lower = c(0))
#check for convergence
Opt4$message

Opt_par4 <- setNames(Opt3$par, c("kr"))
Opt_par4
#Nelder-Mead
# kr 
# 0.0006077762 
#SANN
# kr 
# 0.0006077762 
#Brent
# kr 
# 0.0006077762 



#I to R vax
#(1-kappa) *lambda = gamma *lambda = gl
omicron_vax_recovered <- omicron %>%
  pull(R_vax)

IR_vax <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dI <- -gl * I
    dR <- gl * I
    list(c(dI, dR))
  })
}


#up to here the code should be good

#residual sum of squares function
RSS5 <- function(parameters) {
  names(parameters) <- c("gl")
  out <- ode(y = init5, times = days, func = IR_vax, parms = parameters)
  fit <- out[ , 3]
  sum((omicron_vax_recovered - fit)^2)
}

init5 <- c(I = omicron$I_vax_active[1], R = omicron$R_vax[1])
#Nelder-Mead”, “BFGS”, “CG”, “L-BFGS-B”, “SANN”, “Brent"
#optimization
Opt5 <- optim(c(0.5), RSS5, method = "Brent", upper = c(1), lower = c(0))
#check for convergence
Opt5$message

Opt_par5 <- setNames(Opt5$par, c("gl"))
Opt_par5

#SANN
# gl 
# 0.1320521 
#Nelder-Mead
# gl 
# 0.1320312
#Brent
# gl 
# 0.1320479 


#I to R unvax
omicron_unvax_recovered <- omicron %>%
  pull(R_unvax)
IR_unvax <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dI <- - gl * I
    dR <- gl * I
    list(c(dI, dR))
  })
}

#residual sum of squares function
RSS6 <- function(parameters) {
  names(parameters) <- c("gl")
  out <- ode(y = init6, times = days, func = IR_unvax, parms = parameters)
  fit <- out[ , 3]
  sum((omicron_unvax_recovered - fit)^2)
}

init6 <- c(I = omicron$I_unvax_active[1], R = omicron$R_unvax[1])
#Nelder-Mead”, “BFGS”, “CG”, “L-BFGS-B”, “SANN”, “Brent"
#optimization
Opt6 <- optim(c(0.5), RSS6, method = "Brent", upper = c(1), lower = c(0))
#check for convergence
Opt6$message

Opt_par6 <- setNames(Opt6$par, c("gl"))
Opt_par6

#SANN
# gl 
# 0.003799402 
#Nelder-Mead
# gl 
# 0.00380249
#Brent
# gl 
# 0.003800306 

#when were the other waves?

covidvax %>%
  ggplot(aes(x=date)) + geom_point(aes(y=total_cases), color = "black") + geom_point(aes(y=vaccinated_cases), color = "green") + geom_point(aes(y=unvaccinated_cases), color = "red")

#only catching end of alpha wave because that's when the vaccine came out
#I think Omicron wave is around the end of this dataset?
#see if the county dataset has more dates so I can see the trends there

#end goal: plot out ~3 different time points and fiddle with the vaccination rate on each to see its effect



