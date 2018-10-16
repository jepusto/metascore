library(tidyverse)
library(metafor)
devtools::load_all()
rm(list=ls())

studies = 40
mean_effect = 0.2
sd_effect = 0.1
n_sim = n_beta(n_min = 20, n_max = 120, na = 1, nb = 3)
p_thresholds = .025
p_RR = 0.1
test_steps <- .025
# prior_mass <- 2 / 5

dat <- r_SMD(studies, mean_effect, sd_effect, n_sim, p_thresholds = p_thresholds, p_RR = p_RR)
mean(dat$p > .05)
model <- rma(yi = g, vi = Va, data = dat, method = "ML")

booties <- bootstrap_n_sig(model, step = test_steps, return_boot_dist = TRUE)
plot(density(booties), main = "Number of significant studies under null model")
(S <- n_sig(yi = model$yi, sei = sqrt(model$vi), step = test_steps))
mean(S <= booties)
abline(v = S, col = "red")



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


#------------------------------------------------------------
# Compare bootstrap distribution to actual distribution
#------------------------------------------------------------
reps <- 1000
meta_dat <- rerun(reps, r_SMD(studies, mean_effect, sd_effect, n_sim, p_thresholds = p_thresholds, p_RR = p_RR))
test_types <- data_frame(type = "robust", info = "expected")
mods <- map(meta_dat, fit_meta)




test_results <- map_dfr(meta_dat, estimate_effects, 
                        test_types = test_types, .id = "id")

VHSM_score_test(mods[[1]], steps = .025, type = "robust", bootstrap = TRUE, reps = 100)

boot_results <- map_dfr(mods[1:5], VHSM_score_test, 
                        steps = .025, type = "robust", 
                        bootstrap = TRUE, reps = 10,
                        return_boot_dist = TRUE, .id = "id")


ggplot(test_results, aes(Q_score)) + 
  geom_density(fill = "blue", alpha = 0.2) + 
  theme_minimal()
