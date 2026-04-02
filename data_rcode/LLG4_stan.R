library(rstan)
library(tidyverse)


iowa_rec18_dum <- read.csv("data_rcode/iowa18.csv")
# x is for the survival time variable
# z is for the cure rate variable
iowa18_stan <- list(N = nrow(iowa_rec18_dum), t = iowa_rec18_dum$time_to_recurrence,
                    is_censored = iowa_rec18_dum$is_censored,
                    cured = iowa_rec18_dum$cured, x = cbind(1, iowa_rec18_dum[, -c(3,4,5,6,10, 13, 21, 22, 31, 32:38, 41:43)]),
                    z = iowa_rec18_dum[, -c(3,4,5,6,10, 13, 21, 22, 31, 32:43)], K = 22, M = 25)

llg_stan <- stan_model("stan_final/LLG4.stan")
llg_stan_fit <- sampling(llg_stan, data = iowa18_stan, seed = 12345, chains = 4, iter = 2000, init = 0,
                          pars = c("betaU", "betaC", "sigma", "alpha_fg", "lambda_fg", "lp__", "log_lik"))

saveRDS(llg_stan_fit, "results/data_result/llg4.rds")
