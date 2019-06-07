devtools::load_all()
source("paper-3PSM-and-excess-significant-test/simulation-functions.R")

runSim(
  100, 
  studies = 80, 
  mean_effect = 0.2, 
  sd_effect = 0.1,
  n_sim = n_beta(20, 80, 1, 1.5), 
  p_RR = 1,
  test_alpha = .025, 
  methods = c("FE","ML","REML","WLS"),
  score_types = c("TES-norm","TES-binom","parametric","robust"),
  LRT = TRUE, k_min = 3L, 
)

dat <- r_SMD(
  studies = 80, 
  mean_effect = 0.2, 
  sd_effect = 0.1,
  n_sim = n_beta(20,80, 1, 1.5), 
  p_RR = .1
)

fit_meta(dat, method = "FE")
rma(yi = g, vi = Va, data = dat, method = "FE")
fit_ML <- fit_meta(dat, method = "ML")
fit_ML
rma(yi = g, vi = Va, data = dat, method = "ML")
fit_meta(dat, method = "REML")
rma(yi = g, vi = Va, data = dat, method = "REML")
fit_meta(dat, method = "REML", weights = 1 / Va)
rma(yi = g, vi = Va, weights = 1 / Va, data = dat, method = "REML")


test_for_selection(dat)


alpha <- .025
methods <- c("FE","ML","REML","WLS")
score_types <- c("TES-norm","TES-binom","parametric","robust")
LRT <- TRUE
k_min <- 3
tol <- 10^-3
LRT_method <- "L-BFGS-B"

LRT_3PSM(fit_ML)
LRT_VHSM(fit_ML, two_sided = FALSE)

