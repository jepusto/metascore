rm(list=ls())

library(dplyr)
library(purrr)
devtools::load_all()

#--------------------------------------------------------
# Simulation conditions: no selection 
#--------------------------------------------------------
source("R/VHSM-likelihood.R")
source("R/score-tests.R")
source("R/simulation-functions.R")
source("R/wild-bootstrap.R")
source_obj <- ls()

set.seed(20181004)

design_factors <- list(
  studies = c(20, 40, 80, 120, 200),
  n_factor = 2L,
  mean_effect = seq(-0.5, 1.5, 0.1), 
  sd_effect = c(0.0, 0.1, 0.2, 0.4),
  p_thresholds = .025, 
  p_RR = 1L,
  replicate = 1:4
)

lengths(design_factors)
prod(lengths(design_factors))

params <-
  cross_df(design_factors) %>%
  mutate(
    reps = 1000,
    seed = round(runif(1) * 2^30) + 1:n()
  ) %>%
  sample_frac() 

score_test_types <- list(
  type = c("parametric","subscore","robust"), 
  info = c("expected")
) %>%
  cross_df()


n_sim <- n_beta(20, 120, 1, 3)

#--------------------------------------------------------
# run simulations in parallel
#--------------------------------------------------------

library(Pusto)

cluster <- start_parallel(source_obj = source_obj, register = TRUE)

tm <- system.time(
  results <- plyr::mdply(params, .f = runSim, 
                         n_sim = n_sim, 
                         boot_n_sig = TRUE,
                         .parallel = TRUE)
)

tm 
parallel::stopCluster(cluster)



#--------------------------------------------------------
# Save results and details
#--------------------------------------------------------

session_info <- sessionInfo()
run_date <- date()

save(params, results, session_info, run_date, 
     file = "auxiliary/score-test-simulation-results.Rdata")
