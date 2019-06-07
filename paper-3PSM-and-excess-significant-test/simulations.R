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
