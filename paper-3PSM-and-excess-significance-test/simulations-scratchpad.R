rm(list=ls())
source("paper-3PSM-and-excess-significance-test/simulation-functions.R")

res <- runSim(
  20,
  studies = 20,
  mean_effect = 0.1,
  sd_effect = 0,
  n_sim = n_beta(20, 100, 1, 1.5),
  p_thresholds = .025,
  p_RR = 1,
  seed = 408055689
)

res
