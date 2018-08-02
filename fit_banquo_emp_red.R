
rm(list=ls())
# Librairies
library(rstan)
rstan_options(auto_write = TRUE)
library(rethinking)


load("stan_data.Rdata")

# fit banquo model in stan
banquo.fit <- stan(
  file = "banquo_tmp.stan",
  data = stan.data,
  iter = 20,
  chains = 1,
  init = list(list( logmean_alphaii =1,
                    sigma_alphaii =0.001,
                    alphaij_intercept =0.01,
                    alphaij_center =0,
                    alphaij_width =0.01)) ,
  verbose = TRUE
  #   boost_lib = "/usr/include/boost",
  #  eigen_lib = "/usr/include/eigen3"
)
#id <- floor(100000*runif(1))
#save(list = ls(), file = paste0("results/results_", id, ".Rdata"))
