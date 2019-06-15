library(tidyverse)
load("paper-3PSM-and-excess-significance-test/Type-I-error-results.Rdata")

results_wide <-
  results %>%
  unnest() %>%
  select(-reps, -seed, -p_thresholds, -p_RR) %>%
  mutate(studies_fac = factor(studies))


#----------------------------
# non-convergence rates
#----------------------------

non_convergence_rates <-
  results_wide %>%
  filter(type=="LRT")

ggplot(non_convergence_rates, aes(mean_effect, pct_NA, color = studies_fac, linetype = studies_fac)) +
  geom_line() + 
  geom_point() + 
  facet_wrap(~ sd_effect) + 
  theme_light()

#----------------------------
# Uniform significance
#----------------------------

uniformity_rates <-
  results_wide %>%
  filter(type=="LRT") %>%
  select(studies_fac, mean_effect, sd_effect, pct_no_sig, pct_all_sig) %>%
  gather(cat, pct, starts_with("pct_")) %>%
  mutate(
    cat = recode(cat, pct_no_sig = "No sig. studies", pct_all_sig = "All sig. studies")
  )

ggplot(uniformity_rates, aes(mean_effect, pct, color = studies_fac, linetype = studies_fac)) +
  geom_line() + 
  geom_point() + 
  facet_grid(cat ~ sd_effect) + 
  theme_light()

#----------------------------
# TES
#----------------------------

TES <- 
  results_wide %>%
  filter(type %in% c("TES-norm","TES-binom"))

ggplot(TES, aes(mean_effect, reject_025, color = studies_fac, linetype = type, shape = type)) +
  geom_line() + 
  geom_point() +
  geom_hline(yintercept = .025) + 
  facet_grid(model ~ sd_effect, scales = "free_y") + 
  theme_light()

ggplot(TES, aes(mean_effect, reject_050, color = studies_fac, linetype = type, shape = type)) +
  geom_line() + 
  geom_point() +
  geom_hline(yintercept = .05) + 
  facet_grid(model ~ sd_effect, scales = "free_y") + 
  theme_light()

ggplot(TES, aes(mean_effect, reject_100, color = studies_fac, linetype = type, shape = type)) +
  geom_line() + 
  geom_point() +
  geom_hline(yintercept = .10) + 
  facet_grid(model ~ sd_effect, scales = "free_y") + 
  theme_light()

#----------------------------
# LRT
#----------------------------

LRT <- 
  results_wide %>%
  filter(type == "LRT") %>%
  select(studies_fac, mean_effect, sd_effect, starts_with("reject")) %>%
  gather("alpha","reject_rate", starts_with("reject")) %>%
  mutate(alpha = as.numeric(str_sub(alpha,-3,-1)) / 1000)

LRT_alphas <-
  LRT %>%
  select(sd_effect, alpha) %>%
  distinct()

ggplot(LRT, aes(mean_effect, reject_rate, color = studies_fac, linetype = studies_fac, shape = studies_fac)) +
  geom_line() + 
  geom_point() +
  expand_limits(y = 0) + 
  geom_hline(data = LRT_alphas, aes(yintercept = alpha)) + 
  facet_grid(alpha ~ sd_effect, scales = "free_y") + 
  theme_light()


#----------------------------
# parametric score test
#----------------------------

score_parametric <- 
  results_wide %>%
  filter(type == "parametric", model != "FE")

ggplot(score_parametric, aes(mean_effect, reject_025, color = studies_fac, linetype = studies_fac, shape = studies_fac)) +
  geom_line() + 
  geom_point() +
  geom_hline(yintercept = .025) + 
  expand_limits(y = 0) + 
  facet_grid(sd_effect ~ model, scales = "free_y") + 
  theme_light()

ggplot(score_parametric, aes(mean_effect, reject_050, color = studies_fac, linetype = studies_fac, shape = studies_fac)) +
  geom_line() + 
  geom_point() +
  geom_hline(yintercept = .050) + 
  expand_limits(y = 0) + 
  facet_grid(sd_effect ~ model, scales = "free_y") + 
  theme_light()

ggplot(score_parametric, aes(mean_effect, reject_100, color = studies_fac, linetype = studies_fac, shape = studies_fac)) +
  geom_line() + 
  geom_point() +
  geom_hline(yintercept = .100) + 
  expand_limits(y = 0) + 
  facet_grid(sd_effect ~ model, scales = "free_y") + 
  theme_light()

score_parametric <- 
  results_wide %>%
  filter(type == "parametric-tweaked", model != "FE")

ggplot(score_parametric, aes(mean_effect, reject_025, color = studies_fac, linetype = studies_fac, shape = studies_fac)) +
  geom_line() + 
  geom_point() +
  geom_hline(yintercept = .025) + 
  expand_limits(y = 0) + 
  facet_grid(sd_effect ~ model, scales = "free_y") + 
  theme_light()

ggplot(score_parametric, aes(mean_effect, reject_050, color = studies_fac, linetype = studies_fac, shape = studies_fac)) +
  geom_line() + 
  geom_point() +
  geom_hline(yintercept = .050) + 
  expand_limits(y = 0) + 
  facet_grid(sd_effect ~ model, scales = "free_y") + 
  theme_light()

ggplot(score_parametric, aes(mean_effect, reject_100, color = studies_fac, linetype = studies_fac, shape = studies_fac)) +
  geom_line() + 
  geom_point() +
  geom_hline(yintercept = .100) + 
  expand_limits(y = 0) + 
  facet_grid(sd_effect ~ model, scales = "free_y") + 
  theme_light()


#----------------------------
# robust score test
#----------------------------

score_robust <- 
  results_wide %>%
  filter(type == "robust", model != "FE")

ggplot(score_robust, aes(mean_effect, reject_025, color = studies_fac, linetype = studies_fac, shape = studies_fac)) +
  geom_line() + 
  geom_point() +
  geom_hline(yintercept = .025) + 
  expand_limits(y = 0) + 
  facet_grid(sd_effect ~ model, scales = "free_y") + 
  theme_light()

ggplot(score_robust, aes(mean_effect, reject_050, color = studies_fac, linetype = studies_fac, shape = studies_fac)) +
  geom_line() + 
  geom_point() +
  geom_hline(yintercept = .050) + 
  expand_limits(y = 0) + 
  facet_grid(sd_effect ~ model, scales = "free_y") + 
  theme_light()

ggplot(score_robust, aes(mean_effect, reject_100, color = studies_fac, linetype = studies_fac, shape = studies_fac)) +
  geom_line() + 
  geom_point() +
  geom_hline(yintercept = .100) + 
  expand_limits(y = 0) + 
  facet_grid(sd_effect ~ model, scales = "free_y") + 
  theme_light()

score_robust <- 
  results_wide %>%
  filter(type == "robust-tweaked", model != "FE")

ggplot(score_robust, aes(mean_effect, reject_025, color = studies_fac, linetype = studies_fac, shape = studies_fac)) +
  geom_line() + 
  geom_point() +
  geom_hline(yintercept = .025) + 
  expand_limits(y = 0) + 
  facet_grid(sd_effect ~ model, scales = "free_y") + 
  theme_light()

ggplot(score_robust, aes(mean_effect, reject_050, color = studies_fac, linetype = studies_fac, shape = studies_fac)) +
  geom_line() + 
  geom_point() +
  geom_hline(yintercept = .050) + 
  expand_limits(y = 0) + 
  facet_grid(sd_effect ~ model, scales = "free_y") + 
  theme_light()

ggplot(score_robust, aes(mean_effect, reject_100, color = studies_fac, linetype = studies_fac, shape = studies_fac)) +
  geom_line() + 
  geom_point() +
  geom_hline(yintercept = .100) + 
  expand_limits(y = 0) + 
  facet_grid(sd_effect ~ model, scales = "free_y") + 
  theme_light()

#----------------------------
# compare tests
#----------------------------

tests_selected <- 
  results_wide %>%
  filter(type=="LRT" | 
           (type=="parametric" & model == "ML") | 
           (type %in% c("robust","TES-binom","TES-norm") & model == "WLS"))

ggplot(tests_selected, aes(mean_effect, reject_025, color = studies_fac, shape = studies_fac)) +
  geom_line() + 
  geom_point() +
  expand_limits(y=0) + 
  geom_hline(yintercept = .025) + 
  facet_grid(type ~ sd_effect, scales = "free_y") + 
  theme_light()

ggplot(tests_selected, aes(mean_effect, reject_050, color = studies_fac, shape = studies_fac)) +
  geom_line() + 
  geom_point() +
  expand_limits(y=0) + 
  geom_hline(yintercept = .050) + 
  facet_grid(type ~ sd_effect, scales = "free_y") + 
  theme_light()

ggplot(tests_selected, aes(mean_effect, reject_100, color = studies_fac, shape = studies_fac)) +
  geom_line() + 
  geom_point() +
  expand_limits(y=0) + 
  geom_hline(yintercept = .100) + 
  facet_grid(type ~ sd_effect, scales = "free_y") + 
  theme_light()
