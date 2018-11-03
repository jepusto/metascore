rm(list=ls())

library(dplyr)
library(purrr)

#--------------------------------------------------------
# Simulation conditions: no selection 
#--------------------------------------------------------
source("R/VHSM-likelihood.R")
source("R/fit-VHSM.R")
source("R/score-tests.R")
source("R/simulation-functions.R")
source("R/wild-bootstrap.R")
source_obj <- ls()

set.seed(20181004)

design_factors <- list(
  studies = c(20, 40, 80, 120),
  mean_effect = seq(-0.5, 1.5, 0.1), 
  sd_effect = c(0.0, 0.1, 0.2, 0.4),
  p_thresholds = .025, 
  p_RR = seq(0, 1, 0.1),
  replicate = 1:4
)

lengths(design_factors)
prod(lengths(design_factors))

params <-
  cross_df(design_factors) %>%
  filter(p_RR == 0 | mean_effect %in% c(0, 0.4, 0.8)) %>%
  mutate(
    reps = 50,
    seed = round(runif(1) * 2^30) + 1:n()
  ) %>%
  sample_frac() 

nrow(params)

score_test_types <- list(
  two_sided = FALSE,
  type = c("parametric","robust"), 
  info = "expected",
  prior_mass = c(0, 0.5)
) %>%
  cross_df() %>%
  filter(type == "robust" | prior_mass == 0)

LRT_types <- list(
  two_sided = FALSE,
  k_min = c(0L, 2L)
) %>%
  cross_df()

n_sim <- n_beta(20, 120, 1, 3)

#--------------------------------------------------------
# run simulations in parallel
#--------------------------------------------------------
library(Pusto)

# cluster <- start_parallel(source_obj = source_obj)
# 
# results <- 
#    evaluate_by_row(
#     params[1:20,], 
#     runSim, 
#     n_sim = n_sim, 
#     score_test_types = score_test_types,
#     LRT_types = LRT_types,
#     boot_n_sig = TRUE,
#     boot_qscore = FALSE
#   )

cluster <- start_parallel(source_obj = source_obj, setup = "register")

system.time(
  results <- plyr::mdply(
    params[1:544,], 
    runSim,
    n_sim = n_sim,
    score_test_types = score_test_types,
    LRT_types = LRT_types,
    boot_n_sig = TRUE,
    boot_qscore = FALSE,
    .parallel = TRUE)
)

stop_parallel(cluster)

#--------------------------------------------------------
# Save results and details
#--------------------------------------------------------

session_info <- sessionInfo()
run_date <- date()

save(params, results, session_info, run_date, 
     file = "auxiliary/score-test-simulation-results.Rdata")

