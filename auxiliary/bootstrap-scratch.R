library(tidyverse)
devtools::load_all()
rm(list=ls())

studies = 80
mean_effect = 1.0
sd_effect = 0.1
n_sim = n_beta(n_min = 20, n_max = 120, na = 1, nb = 3)
p_thresholds = .025
p_RR = 1
test_steps <- .025
# prior_mass <- 2 / 5

dat <- r_SMD(studies, mean_effect, sd_effect, n_sim, p_thresholds = p_thresholds, p_RR = p_RR)
mean(dat$p > .05)
model <- rma(yi = g, vi = Va, data = dat, method = "ML")

VHSM_score_test(model, steps = .025, type = "robust")

VHSM_score_test(model, steps = .025, type = "robust", bootstrap = TRUE, reps = 100, seed = 4)

VHSM_score_test(model, steps = .025, type = "robust", bootstrap = TRUE, seed = 4)

wild_bootstrap(model, 
               stat = VHSM_score_test, 
               steps = .025, type = "robust", return = "Q", 
               plot = TRUE, seed = 4)

VHSM_score_test(model, steps = .025)

wild_bootstrap(model, 
               stat = VHSM_score_test, 
               steps = .025, return = "Q", 
               reps = 1999L,
               plot = TRUE)

boots <- wild_bootstrap(model, 
               stat = VHSM_score_test, 
               steps = .025, type = "robust", return = "Q", 
               plot = FALSE, 
               return_boot_dist = TRUE)

