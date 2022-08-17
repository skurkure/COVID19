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
covidvax <- read.csv("covid19postvaxstatewidestats.csv")
head(covidvax$date, 50)
#add total cases column
covidvax <- covidvax %>%
  mutate(total_cases = unvaccinated_cases + vaccinated_cases)
covidvax <- covidvax %>%
  mutate(Submission_Date = as.Date(date))
#population 
N <- covidvax$population[151] #total population on July 1
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

#fit the data


#Full Model
sird_expanded <- function(time, init, parameters) {
  with(as.list(c(init, parameters)), {
    dS_vax <- -beta_vax * S_vax * I_vax / N + epsilon * sigma * S_unvax / N + alpha * R 
    dI_vax <- beta_vax * S_vax * I_vax / N - kappa_vax * rho_vax * I_vax - (1 - kappa_vax) * lambda_vax * R 
    dS_unvax <- -beta_unvax * S_unvax * I_unvax / N - epsilon * (1 - sigma) * S_unvax - epsilon * sigma * S_unvax
    dI_unvax <- beta_unvax * S_unvax * I_unvax / N - kappa_unvax * rho_unvax * I_unvax - (1 - kappa_unvax) * lambda_unvax * R 
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

covidvax <- covidvax %>%
  mutate(population = population_unvaccinated + population_vaccinated)
covidvax <- covidvax %>%
  mutate(cum_unvax_deaths = cumsum(unvaccinated_deaths))
covidvax <- covidvax %>%
  mutate(cum_vax_deaths = cumsum(vaccinated_deaths))
covidvax <- covidvax %>%
  mutate(cum_unvax_deaths = cumsum(unvaccinated_deaths))

SIR <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- -beta * I * S / N
    dI <- beta * I * S / N - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}
#set arbitrary initial conditions-try July when cases beginning spiking
N <- 39350000
july <- covidvax %>%
  filter(Submission_Date >="2021-07-01" & Submission_Date <= "2021-7-31")
july$date
Day <- 1:(length(july))
which(covidvax$Submission_Date == "2021-07-01") #151
july_data <- covidvax %>%
  filter(Submission_Date >= "2021-07-01" & Submission_Date < "2021-08-01") %>%
  pull(cum_cases)
july_data
init <- c(S = covidvax$population_unvaccinated[151], I = covidvax$cum_cases[151], R = covidvax$population_vaccinated[151])
#residual sum of squares function
RSS <- function(parameters) {
  names(parameters) <- c("beta", "gamma")
  out <- ode(y = init, times = Day, func = SIR, parms = parameters)
  fit <- out[ , 3]
  sum((july_data - fit)^2)
}

#Nelder-Mead”, “BFGS”, “CG”, “L-BFGS-B”, “SANN”, “Brent"

#optimization
Opt <- optim(c(0.5, 0.5), RSS, method = "L-BFGS-B", lower = c(0,0), upper = c(1,1))

#check for convergence
Opt$message

Opt_par <- setNames(Opt$par, c("beta", "gamma"))
Opt_par
# beta     gamma 
# 0.5410568 0.1375219 
#so this is for total cases

#vax pop only
SIRD_vax <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- -beta * I * S / N
    dI <- beta * I * S / N - gamma * I - kappa * I
    dR <- gamma * I
    dD <- kappa * I
    list(c(dS, dI, dR, dD))
  })
}

#this time get the vax data only


july_vax_data <- covidvax %>%
  filter(Submission_Date >= "2021-07-01" & Submission_Date <= "2021-07-31") %>%
  pull(cum_vax)
init1 <- c(S = covidvax$population_vaccinated[151]-covidvax$cum_vax[134], I = covidvax$cum_vax[151], R = covidvax$cum_vax[134], D = covidvax$cum_vax_deaths[151])
#residual sum of squares function
RSS1 <- function(parameters) {
  names(parameters) <- c("beta", "gamma", "kappa")
  out <- ode(y = init1, times = Day, func = SIRD_vax, parms = parameters)
  fit <- out[ , 3]
  sum((july_vax_data - fit)^2)
}
#Nelder-Mead”, “BFGS”, “CG”, “L-BFGS-B”, “SANN”, “Brent
#optimization
Opt1 <- optim(c(0.5, 0.5, 0.5), RSS1, method = "CG")
#check for convergence
Opt1$message

Opt_par1 <- setNames(Opt1$par, c("beta", "gamma", "kappa"))
Opt_par1
#L-BFGS-B
# beta     gamma     kappa
# 0.7773811 0.2226248 0.2226248
#BFGS
# beta     gamma     kappa
# 0.6948938 0.1959882 0.1959880
#how to interpret:
#the infectious person infects one person on average every 2 days and is infectious for 5 days.
#for example, β=1/2 day⁻¹=0.5 day⁻¹ and γ=1/5 day⁻¹=0.2 day⁻¹
#time in days for predictions
sird_starting_date <- "2020-07-01"
tail(covidvax, 1) #Dec 4, 2021, row 307
sird_ending_date <- "2021-12-04"
t <- 151:307
#get fitted values from SIR model
fitted_cumulative_incidence <- data.frame(ode(y = init1, times = t, func = SIRD, parms = Opt_par))
fitted_cumulative_incidence
#not really sure how R goes down. also the gamma and kappa are not equal, but the computer wouldn't know that without data
#to here it doesn't work

#try just SIR
SIR_vax_a <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- -beta * I * S / N
    dI <- beta * I * S / N - gamma * I
    dR <- gamma * I
    list(c(dS, dI, dR))
  })
}
#residual sum of squares function
RSSa <- function(parameters) {
  names(parameters) <- c("beta", "gamma")
  out <- ode(y = inita, times = Day, func = SIR_vax_a, parms = parameters)
  fit <- out[ , 3]
  sum((july_vax_data - fit)^2)
}
inita <- c(S = covidvax$population_vaccinated[151]-covidvax$cum_vax[134], I = covidvax$cum_vax[151], R = covidvax$cum_vax[134])
#Nelder-Mead”, “BFGS”, “CG”, “L-BFGS-B”, “SANN”, “Brent"
#optimization
Opta <- optim(c(0.5, 0.5), RSSa, method = "CG")
#check for convergence
Opta$message

Opt_para <- setNames(Opta$par, c("beta", "gamma"))
Opt_para

#CG
# beta     gamma 
# 0.6185610 0.3427913 
#these are probably the correct rates

#nonvax data now


july_unvax_data <- covidvax %>%
  filter(Submission_Date >= "2021-07-01" & Submission_Date <= "2021-07-31") %>%
  pull(cum_unvax)
init2 <- c(S = covidvax$population_unvaccinated[151]-covidvax$cum_unvax[134], I = covidvax$cum_unvax[151], R = covidvax$cum_unvax[134], D = covidvax$cum_unvax_deaths[151])
#residual sum of squares function
RSS2 <- function(parameters) {
  names(parameters) <- c("beta", "gamma", "kappa")
  out <- ode(y = init2, times = Day, func = SIRD_vax, parms = parameters)
  fit <- out[ , 3]
  sum((july_unvax_data - fit)^2)
}
#Nelder-Mead”, “BFGS”, “CG”, “L-BFGS-B”, “SANN”, “Brent
#optimization
Opt2 <- optim(c(0.5, 0.5, 0.5), RSS2, method = "L-BFGS-B", lower = c(0,0), upper = c(1,1))
#check for convergence
Opt2$message

Opt_par2 <- setNames(Opt2$par, c("beta", "gamma", "kappa"))
Opt_par2

# L-BFGS-B
# beta      gamma      kappa 
# 0.30941737 0.04769642 0.04769734 
# beta      gamma      kappa 
# 0.10769907 0.01376049 0.01376040 
#BFGS
# beta         gamma         kappa 
# 1.426830e-04  2.392663e+04 -2.392664e+04 
#CG
# beta      gamma      kappa 
# 0.41064130 0.06317886 0.06317856 

#just SIR model-unvax
#residual sum of squares function
initb <- c(S = covidvax$population_unvaccinated[151]-covidvax$cum_unvax[134], I = covidvax$cum_unvax[151], R = covidvax$cum_unvax[134])
RSSb <- function(parameters) {
  names(parameters) <- c("beta", "gamma")
  out <- ode(y = initb, times = Day, func = SIR_vax_a, parms = parameters)
  fit <- out[ , 3]
  sum((july_unvax_data - fit)^2)
}
#Nelder-Mead”, “BFGS”, “CG”, “L-BFGS-B”, “SANN”, “Brent
#optimization
Optb <- optim(c(0.5, 0.5), RSSb, method = "CG")
#check for convergence
Optb$message

Opt_parb <- setNames(Optb$par, c("beta", "gamma"))
Opt_parb
#CG
# beta     gamma 
# 0.5896847 0.1848488 


SIRD <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dS <- -beta * I * S / N + sigma * R
    dI <- beta * I * S / N - gamma * I - kappa * I
    dR <- gamma * I - sigma * R
    dD <- kappa * I
    list(c(dS, dI, dR, dD))
  })
}



covidvax <- covidvax %>%
  mutate(population = population_unvaccinated + population_vaccinated)





#try getting R and D from here
#start with D_vax
july_vax_death <- july %>%
  pull(cum_vax_deaths)
ID_vax <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dI <- -kr * I
    dD <- kr * I
    list(c(dI, dD))
  })
}

days <- 1:length(july_vax_death)
days
#residual sum of squares function
RSS3 <- function(parameters) {
  names(parameters) <- c("kr")
  out <- ode(y = init3, times = days, func = ID_vax, parms = parameters)
  fit <- out[ , 3]
  sum((july_vax_death - fit)^2)
}

init3 <- c(I = july$cum_vax[1], D = july$cum_vax_deaths[1])
#Nelder-Mead”, “BFGS”, “CG”, “L-BFGS-B”, “SANN”, “Brent"
#optimization
Opt3 <- optim(c(0.5), RSS3, method = "Nelder-Mead")
#check for convergence
Opt3$message

Opt_par3 <- setNames(Opt3$par, c("kr"))
Opt_par3

#SANN
# kr 
# 0.000186189 
#Brent
# kr 
# 0.0001941406 
#Nelder-Mead
# kr 
# 0.0001953125 
july %>%
  ggplot((aes(x=date))) + geom_point(aes(y=cum_vax), color = "green") + geom_point(aes(y=cum_vax_deaths), color = "black")
july %>%
  ggplot((aes(x=date))) + geom_point(aes(y=cum_vax_deaths), color = "black")


#[2,3]
# > Opt_par3
# kr 
# 0.00265921 
# > Opt_par3
# kr 
# 0.001330355 


#I to D unvax
july_unvax_death <- july %>%
  pull(cum_unvax_deaths)
ID_unvax <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dI <- -kr * I
    dD <- kr * I
    list(c(dI, dD))
  })
}

days <- 1:length(july_vax_death)
days
#residual sum of squares function
RSS4 <- function(parameters) {
  names(parameters) <- c("kr")
  out <- ode(y = init4, times = days, func = ID_unvax, parms = parameters)
  fit <- out[ , 3]
  sum((july_unvax_death - fit)^2)
}

init4 <- c(I = july$cum_unvax[1], D = july$cum_unvax_deaths[1])
#Nelder-Mead”, “BFGS”, “CG”, “L-BFGS-B”, “SANN”, “Brent"
#optimization
Opt4 <- optim(c(0.5), RSS4, method = "Nelder-Mead")
#check for convergence
Opt4$message

Opt_par4 <- setNames(Opt3$par, c("kr"))
Opt_par4
#SANN
# kr 
# 0.0001953125 
#Brent
# kr 
# 0.0001953125 
#Nelder-Mead
# kr 
# 0.0001953125 

#for recovered, i have to subtract the deaths from the cases. i don't think they
#occur on the same day, so I will subtract deaths from a few days after from when
#cases are reported

#i will plot out cases and deaths to get a general sense of the two trends

july_vax_cases_deaths_plot <- july %>%
  ggplot(aes(x=date)) + geom_point(aes(y=cum_vax), color = "green") + geom_point(aes(y=cum_vax_deaths), color= "black")
july_vax_cases_deaths_plot
july %>%
  ggplot(aes(x=date)) + geom_point(aes(y=cum_vax_deaths), color= "black")
#can't really tell anything

#according to https://pubmed.ncbi.nlm.nih.gov/34154563/#:~:text=The%20estimated%20mean%20time%20between,%3B%20heterogeneity%20P%20%3D%200.85).,
#the average time to death is 18 days. i'll take cases from 18 days prior and subtract
#deaths from the day 18 days later.

#nvm i'll just subtract one column from the other, it's simpler that way

covidvax <- covidvax %>%
  mutate(cum_vax_recovered = cum_vax - cum_vax_deaths)
covidvax <- covidvax %>%
  mutate(cum_unvax_recovered = cum_unvax - cum_unvax_deaths)
july <- covidvax %>%
  filter(Submission_Date >="2021-07-01" & Submission_Date <= "2021-7-31")
july_vax_recovered <- july$cum_vax_recovered
#I to R vax
#(1-kappa) *lambda = gamma *lambda = gl
IR_vax <- function(time, state, parameters) {
  par <- as.list(c(state, parameters))
  with(par, {
    dI <- - gl * I
    dR <- gl * I
    list(c(dI, dR))
  })
}


#residual sum of squares function
RSS5 <- function(parameters) {
  names(parameters) <- c("gl")
  out <- ode(y = init5, times = days, func = IR_vax, parms = parameters)
  fit <- out[ , 3]
  sum((july_vax_recovered - fit)^2)
}

init5 <- c(I = july$cum_vax[1], R = july$cum_vax_recovered[1])
#Nelder-Mead”, “BFGS”, “CG”, “L-BFGS-B”, “SANN”, “Brent"
#optimization
Opt5 <- optim(c(0.5), RSS5, method = "Nelder-Mead")
#check for convergence
Opt5$message

Opt_par5 <- setNames(Opt5$par, c("gl"))
Opt_par5

#Nelder-Mead
# gl 
# 0.1765625 
#BFGS
# gl 
# 0.1765625 
#SANN
# gl 
#0.1765884 

#I to R unvax
july_unvax_recovered <- july$cum_unvax_recovered
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
  sum((july_unvax_recovered - fit)^2)
}

init6 <- c(I = july$cum_unvax[1], R = july$cum_unvax_recovered[1])
#Nelder-Mead”, “BFGS”, “CG”, “L-BFGS-B”, “SANN”, “Brent"
#optimization
Opt6 <- optim(c(0.5), RSS6, method = "BFGS")
#check for convergence
Opt6$message

Opt_par6 <- setNames(Opt5$par, c("gl"))
Opt_par6

#SANN
# gl 
# 0.1765884 
#Nelder-Mead
# gl 
# 0.1765884 
#BFGS
# gl 
# 0.1765884 

#when were the other waves?

covidvax %>%
  ggplot(aes(x=date)) + geom_point(aes(y=total_cases), color = "black") + geom_point(aes(y=vaccinated_cases), color = "green") + geom_point(aes(y=unvaccinated_cases), color = "red")

#only catching end of alpha wave because that's when the vaccine came out
#I think Omicron wave is around the end of this dataset?
#see if the county dataset has more dates so I can see the trends there

#end goal: plot out ~3 different time points and fiddle with the vaccination rate on each to see its effect