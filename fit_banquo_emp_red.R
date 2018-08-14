
# clear the workspace
rm(list=ls())

# load required libraries
library(rstan); rstan_options(auto_write = TRUE)
library(rethinking)

# read in the data
load("stan_data.Rdata")

# compile the model but do nothing else
banquo.model <- stan_model(
  file="banquo.stan",
  verbose=FALSE
  # boost_lib = "/usr/include/boost",
  # eigen_lib = "/usr/include/eigen3"
)

# fit banquo model in stan
banquo.fit <- sampling(
  banquo.model,
  data = stan.data,
  iter = 20,
  chains = 1,
  init = list(list(
    logmean_alphaii = 1,
    sigma_alphaii = 0.001,
    alphaij_intercept = 0.01,
    alphaij_center = 0,
    alphaij_width = 0.01
  ))
)

#id <- floor(100000*runif(1))
#save(list = ls(), file = paste0("results/results_", id, ".Rdata"))
