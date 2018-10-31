library(tidyverse)
library(metafor)
devtools::load_all()
rm(list=ls())

studies = 80
mean_effect = 0.2
sd_effect = 0.1
n_sim = n_beta(n_min = 20, n_max = 120, na = 1, nb = 3)
p_thresholds = .025
p_RR = 0.2
test_steps <- .025
score_test_types <- data_frame(type = "robust", info = "expected", prior_mass = 0)

dat <- r_SMD(studies, mean_effect, sd_effect, n_sim, p_thresholds = p_thresholds, p_RR = p_RR)
mean(dat$p > .05)
model <- rma(yi = g, vi = Va, data = dat, weights = 1 / Va, method = "REML")

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

boots <- wild_bootstrap(model, 
               stat = VHSM_score_test, 
               steps = .025, type = "robust", return = "Q", 
               plot = FALSE, 
               return_boot_dist = TRUE)


#------------------------------------------------------------
# Compare bootstrap distribution to actual distribution
#------------------------------------------------------------

studies = 80
mean_effect = 1.2
sd_effect = 0.1
n_sim = n_beta(n_min = 20, n_max = 120, na = 1, nb = 3)
p_thresholds = .025
p_RR = 1
test_steps <- .025
score_test_types <- data_frame(type = "robust", info = "expected", prior_mass = 0.5)

reps <- 1000
meta_dat <- rerun(reps, r_SMD(studies, mean_effect, sd_effect, n_sim, p_thresholds = p_thresholds, p_RR = p_RR))

# Score wild bootstrap

test_results <- map_dfr(meta_dat, estimate_effects, 
                        score_test_types = score_test_types,
                        boot_n_sig = FALSE, 
                        boot_qscore = FALSE, 
                        .id = "id")

f <- function(dat, reps = 100) {
  mod <- fit_meta_ML(dat)
  wild_bootstrap(mod, stat = VHSM_score_test, 
                 steps = .025, type = "robust", prior_mass = 0.5, 
                 return = "Q", reps = reps,
                 return_boot_dist = TRUE)
}

# f(meta_dat[[1]])

boot_results <- 
  map(meta_dat[1:50], f, reps = 200) %>%
  map_dfr(~ data.frame(Stat = .), .id = "id")

boot_quants <-
  boot_results %>%
  group_by(id) %>%
  summarise(q95 = quantile(Stat, .95))

ggplot(test_results, aes(Stat)) + 
  geom_density(fill = "blue", alpha = 0.2) + 
  geom_vline(xintercept = quantile(test_results$Stat, .95), color = "blue") + 
  stat_function(fun = dchisq, args = list(df = 1), color = "red", size = 2) + 
  coord_cartesian(xlim = c(0, 10), ylim = c(0,1.5)) + 
  geom_density(data = boot_results, aes(Stat), alpha = 0.2, fill = "orange") +
  # geom_vline(xintercept = quantile(boot_results$Stat, .95), color = "orange") +
  # geom_density(data = boot_results, aes(Stat, group = factor(id)), alpha = 0.2, color = "darkgrey") +
  # geom_vline(data = boot_quants, aes(xintercept = q95), color = "darkgrey") +
  theme_minimal() + 
  theme(legend.position = "none")


# n-sig wild bootstrap

test_results <- map_dfr(meta_dat, estimate_effects, 
                        boot_n_sig = TRUE, .id = "id")


bootstrap_n_sig(mods[[1]], step = .025, reps = 100)

boot_results <- map(mods[1:10], bootstrap_n_sig, 
                        step = .025,
                        return_boot_dist = TRUE)

boot_df <- map_dfr(boot_results, ~ data.frame(Stat = .), .id = "id")


test_results %>%
  filter(info == "ML") %>%
ggplot(aes(Stat)) + 
  geom_density(fill = "blue", alpha = 0.2) + 
  geom_density(data = boot_df, aes(Stat), alpha = 0.2, fill = "purple") + 
  geom_density(data = boot_df, aes(Stat, color = factor(id))) + 
  theme_minimal()
