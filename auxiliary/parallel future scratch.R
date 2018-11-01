library(tidyverse)
library(furrr)
devtools::load_all(".")
rm(list=ls())

design_factors <- list(
  studies = c(20, 40, 80, 120),
  mean_effect = seq(-0.5, 1.5, 0.1), 
  sd_effect = c(0.0, 0.1, 0.2, 0.4),
  p_thresholds = .025, 
  p_RR = seq(0, 1,0.1),
  replicate = 1:4
)

lengths(design_factors)
prod(lengths(design_factors))

params <-
  cross_df(design_factors) %>%
  filter(p_RR == 0 | mean_effect %in% c(0, 0.4, 0.8)) %>%
  mutate(
    reps = 10,
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


evaluate_by_row <- function(params, sim_function, ..., .progress = FALSE) {
  
  results_list <- 
    params %>%
    furrr::future_pmap(., .f = sim_function, 
                       ..., .progress = .progress)
  
  params %>%
    mutate(..results = results_list) %>%
    unnest(..results)

}

results <- 
  params %>% 
  filter(row_number() <= 5) %>%
  evaluate_by_row(
    sim_function = runSim,
    n_sim = n_sim, 
    score_test_types = score_test_types,
    LRT_types = LRT_types,
    boot_n_sig = FALSE,
    boot_qscore = FALSE
  )
