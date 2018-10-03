rm(list=ls())

library(dplyr)
library(purrr)
devtools::load_all()


#------------------------------------------
# sample size distributions
#------------------------------------------

# sample sizes from Carter et al. (2015) via Inzlicht et al. (2017)
# Data available from https://osf.io/fcts8/

n0 <- c(129L, 86L, 90L, 54L, 71L, 69L, 65L, 56L, 54L, 74L, 45L, 43L, 41L, 40L, 
        41L, 41L, 56L, 39L, 37L, 35L, 57L, 38L, 26L, 38L, 38L, 37L, 37L, 42L, 
        25L, 35L, 30L, 23L, 33L, 33L, 33L, 51L, 32L, 33L, 31L, 33L, 31L, 29L, 
        32L, 40L, 30L, 30L, 29L, 30L, 29L, 28L, 24L, 27L, 17L, 26L, 24L, 25L, 
        26L, 27L, 24L, 24L, 24L, 24L, 25L, 25L, 20L, 24L, 22L, 22L, 20L, 21L, 
        21L, 20L, 20L, 20L, 20L, 20L, 19L, 19L, 19L, 19L, 18L, 18L, 18L, 18L, 
        17L, 17L, 15L, 16L, 16L, 15L, 15L, 15L, 15L, 15L, 15L, 16L, 15L, 14L, 
        16L, 14L, 14L, 14L, 14L, 14L, 14L, 13L, 14L, 14L, 13L, 13L, 12L, 11L, 
        10L, 11L, 11L, 11L, 10L, 10L)
n1 <- c(122L, 109L, 90L, 108L, 71L, 68L, 63L, 60L, 53L, 27L, 45L, 42L, 42L, 
        40L, 38L, 38L, 22L, 38L, 40L, 41L, 19L, 38L, 50L, 38L, 37L, 36L, 36L, 
        27L, 44L, 33L, 37L, 44L, 33L, 33L, 33L, 15L, 33L, 29L, 31L, 28L, 30L, 
        32L, 29L, 20L, 30L, 30L, 30L, 29L, 29L, 27L, 31L, 27L, 34L, 24L, 26L, 
        25L, 23L, 22L, 24L, 24L, 24L, 23L, 21L, 21L, 26L, 20L, 22L, 21L, 23L, 
        21L, 20L, 20L, 20L, 20L, 20L, 20L, 20L, 19L, 19L, 18L, 18L, 18L, 18L, 
        18L, 16L, 16L, 18L, 16L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 13L, 14L, 
        14L, 12L, 14L, 14L, 14L, 14L, 13L, 13L, 14L, 13L, 13L, 13L, 12L, 12L, 
        11L, 12L, 11L, 11L, 10L, 10L, 10L)
Carter_total_n <- n0 + n1

Carter_n <- round((n0 + n1) / 2)


# dat <- r_SMD(studies = 100, mean_effect = 0.1, sd_effect = 0.1,
#              n_sim = n_empirical(Carter_n),
#              p_thresholds = .025, p_RR = 0.5)
# studies <- nrow(dat)

test_types <- 
  list(type = c("parametric","robust"), info = c("expected")) %>%
  cross()
  
# runSim(reps = 1000, studies = 100, mean_effect = 0.2, sd_effect = 0.1,
#        n_sim = n_empirical(Carter_n), n_factor = 3L, 
#        p_thresholds = .025, p_RR = 1, test_types = test_types)

#--------------------------------------------------------
# Simulation conditions: no selection 
#--------------------------------------------------------
source("R/VHSM-likelihood.R")
source("R/score-tests.R")
source("R/simulation-functions.R")
source_obj <- ls()

set.seed(20181002)

design_factors <- list(
  studies = c(20, 40, 80, 120),
  n_factor = 2L,
  mean_effect = seq(-0.5, 1, 0.1), 
  sd_effect = c(0.0, 0.01, 0.1, 0.2, 0.4),
  p_thresholds = .025, 
  p_RR = 1L,
  replicate = 1:8
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

#--------------------------------------------------------
# run simulations in parallel
#--------------------------------------------------------

library(Pusto)

cluster <- start_parallel(source_obj = source_obj, register = TRUE)

tm <- system.time(
  results <- plyr::mdply(params, .f = runSim, 
                         n_sim = n_empirical(Carter_n), 
                         test_types = test_types,
                         .parallel = TRUE)
)

tm 
parallel::stopCluster(cluster)



#--------------------------------------------------------
# Save results and details
#--------------------------------------------------------

session_info <- sessionInfo()
run_date <- date()

save(params, results, Carter_total_n, session_info, run_date, 
     file = "auxiliary/score-test-simulation-results.Rdata")
