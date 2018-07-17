
# setwd("~/Desktop/R-Desktop/Banquo2")
rm(list=ls())

# libraries and options for them
library(rstan)
rstan_options(auto_write = TRUE)

# load the rethinking package which has some useful tools for analyzing stan models
library(rethinking)

# note that for more than one chain it can be parallelized across chains with:
# options(mc.cores = parallel::detectCores())

# define a few convenience variables for our fake data
N <- 50
S <- 15
T <- 1
Nsites <- 50
sites <- as.numeric(gl(Nsites,N))

# fake observed relative abundances
observed <- matrix(runif(S * Nsites), Nsites, S)
observed <- sweep(observed, 2, 1:S, '*')
observed <- sweep(observed, 1, rowSums(observed), '/')

# fake relative abundances from traitspace
traitspace <- matrix(1.0/S, Nsites * N, S)

# fake traits
traits <- matrix(runif(S * T), S, T)

# shove everything important into a named list
stan.data <- list(
  N = nrow(traitspace),
  S = ncol(traitspace),
  Nsites = nrow(observed),
  T = ncol(traits),
  sites = sites,
  traits = traits,
  observed = observed,
  traitspace = traitspace
)

# DEBUG: we should think about 'clever' starting values to help convergence of chains

# fit banquo model in stan
banquo.fit <- stan(
  file = "banquo.stan",
  data = stan.data,
  iter = 1000,
  chains = 1,
  # verbose = TRUE,
  boost_lib = "/usr/include/boost",
  eigen_lib = "/usr/include/eigen3"
)
