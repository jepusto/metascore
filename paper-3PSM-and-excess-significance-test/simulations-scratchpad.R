rm(list=ls())
source("paper-3PSM-and-excess-significance-test/simulation-functions.R")

res <- runSim(
  8000,
  studies = 20,
  mean_effect = 0.1,
  sd_effect = 0,
  n_sim = n_beta(20, 100, 1, 1.5),
  p_thresholds = .025,
  p_RR = 1,
  seed = 408055689
)

res

# debug Type I error results

load("paper-3PSM-and-excess-significance-test/Type-I-error-results.Rdata")

new_params <- 
  results %>% 
  filter(map_lgl(res, is.null)) %>%
  select(-res)

new_params

reps <- 8000
studies <- 20
mean_effect <- 0.1
sd_effect <- 0
n_sim <- n_beta(20, 100, 1, 1.5)
p_thresholds <- .025
p_RR <- 1
seed <- 408055689
test_alpha <- .025
methods <- c("FE","ML","REML","WLS")
score_types <- c("TES-norm","TES-binom","parametric","robust")


set.seed(seed)

dats <- 
  rerun(reps, {
    r_SMD(studies, mean_effect, sd_effect, n_sim, 
          p_thresholds = p_thresholds, p_RR = p_RR)
  })

test_res <- 
  map(dats, .f = possibly(test_for_selection, NULL), 
      alpha = test_alpha,
      methods = methods,
      score_types = score_types)  

test_res %>%
  bind_rows() %>%
  group_by(model, type) %>% 
  summarise(
    pct_NA = mean(is.na(p_val)),
    pct_no_sig = mean(sig==0),
    pct_all_sig = mean(sig==1),
    reject_025 = mean(p_val < .025, na.rm = TRUE),
    reject_050 = mean(p_val < .050, na.rm = TRUE),
    reject_100 = mean(p_val < .100, na.rm = TRUE)
  )

err <- which(map_lgl(test_res, is.null))
dat <- dats[[err]]
test_for_selection(dat)
test_for_selection(dats[[1]])
alpha <- .025
LRT <- TRUE
k_min <- 3L
tol <- 10^-3
LRT_method <- "L-BFGS-B"

if (length(fit_methods) > 0) {
  mods <- map(fit_methods, ~ fit_meta(dat = dat, method = .))
  names(mods) <- fit_methods
} else {
  mods <- list()
}

if ("WLS" %in% methods) {
  mods <- c(mods, list(WLS = fit_meta(dat, method = "REML", weights = 1 / Va)))
}

score_res <- map_dfr(mods, .f = possibly(simple_scores, otherwise = tibble(type = score_types, Z = NA, p_val = NA)), type = score_types, alpha = alpha, .id = "method")  

simple_scores(mods[[1]], type = score_types, alpha = alpha)
