## Profiling Customer Television Watching Behaviour Using A Novel Bayesian
## Hierarchical Joint Model for Time-to-Event and Count Data
## by Moral et al.

## The code in this script file wrangles the data and fits the Bayesian model
## WARNING: the model takes approximately 50 hours to run

## loading libraries
library(tidyverse)
library(readxl)
library(coda)
library(bayesplot)
library(MCMCvis)
library(runjags)
library(gridExtra)
library(extraDistr)

## loading the customer data
active <- read_xlsx("active.xlsx")
cancelled <- read_xlsx("cancelled.xlsx")

active$RUID_KEY <- as.factor(active$RUID_KEY)
levels(active$RUID_KEY) <- paste0("A:", 1:20)

cancelled$RUID_KEY <- as.factor(cancelled$RUID_KEY)
levels(cancelled$RUID_KEY) <- paste0("C:", 1:20)

## merging both into a single dataset
full_data <- rbind(active, cancelled)
n_customers <- length(unique(full_data$RUID_KEY))

get_times <- function(x) {
  all_splits <- str_split(x, ",") %>%
    unlist %>%
    as.numeric
  n_states <- sum(!is.na(all_splits))
  times <- all_splits %>%
    na.omit %>%
    as.numeric
  return(list("N" = n_states,
              "Time" = times/60))
}

full_data_times <- lapply(full_data$SEQUENCE_CONTENT, get_times)
N_full_data <- sapply(full_data_times, function(x) x$N)
N_full_data_split <- split(N_full_data, full_data$RUID_KEY)
max_N <- max(lapply(N_full_data_split, length) %>% unlist)
N_matrix <- matrix(NA, ncol = n_customers, nrow = max_N)
for(i in 1:n_customers) {
  N_matrix[,i] <- c(N_full_data_split[[i]], rep(NA, max_N - length(N_full_data_split[[i]])))
}
colnames(N_matrix) <- names(N_full_data_split)

Time_full_data <- sapply(full_data_times, function(x) x$Time)
n_sessions <- apply(N_matrix, 2, function(x) sum(!is.na(x)))
Time_array <- array(NA, dim = c(max(N_full_data), max(n_sessions), n_customers))

Drama <- Movies <- Sport <- Factual <- Kids <- 
  Reality <- Crime <- News <- unknown <- array(NA, dim = c(max(N_full_data), max(n_sessions), n_customers))

for(c in 1:n_customers) {
  for(s in 1:n_sessions[c]) {
    for(i in 1:N_matrix[s,c]) {
      Drama[i,s,c] <- ifelse(full_data$GENRE[s*c] == "Drama", 1, 0)
      Movies[i,s,c] <- ifelse(full_data$GENRE[s*c] == "Movies", 1, 0)
      Sport[i,s,c] <- ifelse(full_data$GENRE[s*c] == "Sport", 1, 0)
      Factual[i,s,c] <- ifelse(full_data$GENRE[s*c] == "Factual", 1, 0)
      Kids[i,s,c] <- ifelse(full_data$GENRE[s*c] == "Kids", 1, 0)
      Reality[i,s,c] <- ifelse(full_data$GENRE[s*c] == "Reality", 1, 0)
      Crime[i,s,c] <- ifelse(full_data$GENRE[s*c] == "Crime", 1, 0)
      News[i,s,c] <- ifelse(full_data$GENRE[s*c] == "News", 1, 0)
      unknown[i,s,c] <- ifelse(full_data$GENRE[s*c] == "unknown", 1, 0)
    }
  }
}

Time_full_data_split <- split(Time_full_data, full_data$RUID_KEY)

for(c in 1:n_customers) {
  for(s in 1:n_sessions[c]) {
    Time_array[,s,c] <- c(Time_full_data_split[[c]][[s]],
                          rep(NA, dim(Time_array)[1] - length(Time_full_data_split[[c]][[s]])))
  }
}

## JAGS model
model <- "model.txt"
jagsscript <- cat("
  model {
  
    ## Initialising
    for(c in 1:n_customers) {
      for(s in 1:n_sessions[c]) {
        mu[1,s,c] = Time[1,s,c]
        Time_pred[1,s,c] = Time[1,s,c]
      }
    }

    ## Likelihood
    for(c in 1:n_customers) {
      for(s in 1:n_sessions[c]) {
        ## Likelihood for number of events - truncated Poisson
        N[s,c] ~ dpois(lambda[s,c])T(1,)
        ## linear predictor
        log(lambda[s,c]) = d[c]
        
        for(i in 2:N[s,c]) {
          ## Likelihood for time until next event - gamma
          Time[i,s,c] ~ dgamma(psi[c], psi[c]/mu[i,s,c])
          ## linear predictor
          log(mu[i,s,c]) = eta[i,s,c] + omega[i,s,c]
          ## regressors - random intercepts per genre
          eta[i,s,c] = g[1,c] * Drama[i,s,c] +
                       g[2,c] * Movies[i,s,c] +
                       g[3,c] * Sport[i,s,c] +
                       g[4,c] * Factual[i,s,c] +
                       g[5,c] * Kids[i,s,c] +
                       g[6,c] * Reality[i,s,c] +
                       g[7,c] * Crime[i,s,c] +
                       g[8,c] * News[i,s,c]
          ## AR component
          omega[i,s,c] = p[s,c] * log(Time[i-1,s,c])
        }
      }
      
      ## Priors for random effects
      d[c] ~ dnorm(delta, tau_delta)
      for(genre in 1:n_genre) {
        g[genre,c] ~ dnorm(gamma[genre], tau_gamma[genre])
      }
      for(s in 1:n_sessions[c]) {
        p[s,c] ~ dnorm(phi[c], tau_phi)
      }
    }
    
    ## Priors for model parameters
    delta ~ dnorm(0, 0.001)
    tau_delta ~ dt(0, 0.01, 1)T(0,)
    sd_delta = pow(tau_delta, -1/2)
    for(genre in 1:n_genre) {
      gamma[genre] ~ dnorm(0, 0.001)
      tau_gamma[genre] ~ dt(0, 0.01, 1)T(0,)
      sd_gamma[genre] = pow(tau_gamma[genre], -1/2)
    }
    for(c in 1:n_customers) {
      phi[c] ~ dnorm(0, 0.001)
      psi[c] ~ dgamma(0.001, 0.001)
    }
    tau_phi ~ dt(0, 0.01, 1)T(0,)
    sd_phi = pow(tau_phi, -1/2)
  }
", file = model)

# initialization
initfunction <- function(chain) {
  return(switch(chain,
                "1" = list(
                  .RNG.name="base::Super-Duper",
                  .RNG.seed=1),
                "2" = list(
                  .RNG.name="base::Super-Duper",
                  .RNG.seed=2),
                "3" = list(
                  .RNG.name="base::Super-Duper",
                  .RNG.seed=3),
  ))
}

# JAGS data
## reducing dimension of the data to max. 300 journeys per customer and 300 events per journey
N_matrix2 <- N_matrix[1:300,]
for(i in 1:300) {
  for(j in 1:40) {
    N_matrix2[i,j] <- min(N_matrix2[i,j], 300)
  }
}
jags_data <- list("n_customers" = n_customers, # scalar
                  "n_sessions" = sapply(n_sessions, function(x) min(x, 300)), # vector, customer
                  "n_genre" = 8, # scalar
                  "Time" = Time_array[1:300,1:300,] + .00001, # 3-d array, event x session x customer
                  "Drama" = Drama[1:300,1:300,], # 3-d dummy array, event x session x customer
                  "Movies" = Movies[1:300,1:300,], # 3-d dummy array, event x session x customer
                  "Sport" = Sport[1:300,1:300,], # 3-d dummy array, event x session x customer
                  "Factual" = Factual[1:300,1:300,], # 3-d dummy array, event x session x customer
                  "Kids" = Kids[1:300,1:300,], # 3-d dummy array, event x session x customer
                  "Reality" = Reality[1:300,1:300,], # 3-d dummy array, event x session x customer
                  "Crime" = Crime[1:300,1:300,], # 3-d dummy array, event x session x customer
                  "News" = News[1:300,1:300,], # 3-d dummy array, event x session x customer
                  "N" = N_matrix2) # matrix, session x customer

# parameters to be monitored
jags_params <- c("delta","phi","sd_delta","sd_phi","psi",
                 "d","p","g","gamma","sd_gamma")

# run JAGS in parallel (4 cores)
nChains <- 3
nAdaptSteps <- 1000
nBurninSteps <- 2000
nThinSteps <- 10
nUseSteps <- 6000

## WARNING: the following piece of code takes approximately 50 hours to run
runJagsOut <- run.jags(method = "parallel",
                       model = model,
                       monitor = jags_params,
                       data = jags_data,
                       n.chains = nChains,
                       adapt = nAdaptSteps,
                       burnin = nBurninSteps,
                       sample = ceiling(nUseSteps/nChains),
                       thin = nThinSteps,
                       summarise = FALSE,
                       plots = FALSE,
                       inits = initfunction)

# coda samples - MCMC
coda_samples <- as.mcmc.list(runJagsOut)

# parameter estimates
est_delta <- MCMCsummary(coda_samples, params = c("delta","sd_delta"), round = 4)
est_gamma <- MCMCsummary(coda_samples, params = c("gamma","sd_gamma"), round = 4)
est_phi <- MCMCsummary(coda_samples, params = c("phi","sd_phi"), round = 4)
est_psi <- MCMCsummary(coda_samples, params = c("psi"), round = 4)
est_g <- MCMCsummary(coda_samples, params = c("g"), round = 4)

save.image("moral_et_al_all_model_outputs.RData")