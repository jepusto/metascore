source("R/simulation-functions.R")

dat <- r_SMD(
  studies = 80, 
  mean_effect = 0.2, 
  sd_effect = 0.1,
  n_sim = n_beta(20,80, 1, 1.5), 
  p_RR = .1
)

fit_meta(dat, method = "FE")
rma(yi = g, vi = Va, data = dat, method = "FE")
fit_meta(dat, method = "ML")
rma(yi = g, vi = Va, data = dat, method = "ML")
fit_meta(dat, method = "REML")
rma(yi = g, vi = Va, data = dat, method = "REML")
fit_meta(dat, method = "REML", weights = 1 / Va)
rma(yi = g, vi = Va, weights = 1 / Va, data = dat, method = "REML")

mod <- fit_meta(dat, method = "REML")
alpha <- .025
type <- c("TES-norm","TES-binom","parametric","robust")
