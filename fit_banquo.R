
setwd("~/Desktop/R-Desktop/Banquo2")
rm(list=ls())

# libraries and options for them
library(rstan)
rstan_options(auto_write = TRUE)

# to install the package below
# install.packages(c("coda","mvtnorm","devtools","loo"))
library(devtools)
 devtools::install_github("rmcelreath/rethinking")
library(rethinking)

# define a few convenience variables
N <- 10
S <- 15
T <- 1
Nsites <- 2
sites <- as.numeric(gl(2,5))
# fake observed relative abundances
observed <- matrix(runif(S * Nsites), Nsites, S)
observed <- sweep(observed, 2, 1:S, '*')
observed <- sweep(observed, 1, rowSums(observed), '/')

# fake relative abundances from traitspace
traitspace <- matrix(1.0/S, N, S)

# fake traits
# note not really used yet
traits <- matrix(runif(S * T), S, T)

# shove everything important into a named list
stan.data <- list(
  N = nrow(traitspace),
  S = ncol(traitspace),
Nsites = nrow(observed),
 sites = sites,
  # T = ncol(traits),
  # traits = traits,
  observed = observed,
  traitspace = traitspace
)

# fit Daniel's stan model
banquo.fit <- stan(
  file = "lib/stan/banquo.stan",
  data = stan.data,
  iter = 500,
  chains = 1
)
# fit Daniel's stan model
banquo_simp.fit <- stan(
  file = "lib/stan/banquo.stan",
  data = stan.data,
  iter = 500,
  chains = 1
)
# take a lot at some of the inferred parameters
precis(banquo.fit, depth=2, pars=c("logmean_alphaii", "sigma_alphaii", "alpha"))
