rm(list=ls())
source("paper-3PSM-and-excess-significance-test/simulation-functions.R")

res <- runSim(
  10000,
  studies = 50,
  mean_effect = 0.4,
  sd_effect = 0.1,
  n_sim = n_beta(20, 100, 1, 1.5),
  p_RR = .5,
  test_alpha = .025,
  methods = c("FE","ML","REML","WLS"),
  score_types = c("TES-norm","TES-binom","parametric","robust"),
  LRT = TRUE, k_min = 2L,
)
res
