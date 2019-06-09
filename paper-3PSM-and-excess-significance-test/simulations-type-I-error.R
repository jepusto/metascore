rm(list=ls())
source("paper-3PSM-and-excess-significance-test/simulation-functions.R")

runSim(
  1000,
  studies = 50,
  mean_effect = 1,
  sd_effect = 0.3,
  n_sim = n_beta(20, 80, 1, 1.5),
  p_RR = 1,
  test_alpha = .025,
  methods = c("FE","ML","REML","WLS"),
  score_types = c("TES-norm","TES-binom","parametric","robust"),
  LRT = TRUE, k_min = 3L,
)


source_obj <- ls()

#--------------------------------------------------------
# Simulation conditions
#--------------------------------------------------------
library(tidyverse)

set.seed(20190608)

design_factors <- list(
  studies = c(20, 50, 80),
  mean_effect = seq(0, 1.2, 0.1), 
  sd_effect = c(0, 0.1, 0.2, 0.3, 0.4),
  p_thresholds = .025, 
  p_RR = 1
)

lengths(design_factors)
prod(lengths(design_factors))

params <- 
  cross_df(design_factors) %>%
  mutate(
    reps = 20,
    seed = round(runif(1) * 2^30) + row_number()
  ) %>%
  sample_frac()

params

#--------------------------------------------------------
# run simulations in parallel - furrr
#--------------------------------------------------------

library(future)
plan(multiprocess)
library(furrr)

tm <- system.time(
  results <- 
    params %>%
    mutate(
      res = future_pmap(., .f = runSim, 
                        n_sim = n_beta(20,100,1,1.5),
                        methods = c("FE","ML","REML","WLS"),
                        score_types = c("TES-norm","TES-binom","parametric","robust"),
                        LRT = TRUE, k_min = 2L)
    ) %>%
    unnest()
)

tm
(8000 / unique(results$reps)) * (tm[[3]] / 60^2)

#--------------------------------------------------------
# run simulations in parallel - mdply
#--------------------------------------------------------

library(Pusto)

cluster <- start_parallel(source_obj = source_obj, register = TRUE)

tm <- system.time(
  results <- plyr::mdply(params, .f = runSim, 
                         n_sim = n_beta(20,100,1,1.5), 
                         methods = c("FE","ML","REML","WLS"),
                         score_types = c("TES-norm","TES-binom","parametric","robust"),
                         LRT = TRUE, k_min = 2L, 
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
     file = "paper-3PSM-and-excess-significance-test/Type-I-error-results.Rdata")
