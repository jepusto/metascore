
#------------------------------------------
# sample size distributions
#------------------------------------------

# simulate empirical sample size distribution

n_empirical <- function(n_vec) {
  function(studies, n_factor = 1L) sample(n_vec * n_factor, size = studies, replace = TRUE)
}

# simulate sample size distributions from shifted, scaled beta

n_beta <- function(n_min, n_max, na, nb) {
  n_diff <- n_max - n_min
  function(studies, n_factor = 1L) round(n_factor * (n_min + n_diff * rbeta(n = studies, shape1 = na, shape2 = nb)))
}


#--------------------------------------------
# simulate standardized mean differences
#--------------------------------------------

r_SMD <- function(studies, mean_effect, sd_effect, 
                  n_sim, n_factor = 1L, 
                  p_thresholds = .025, p_RR = .1) {
  
  # sample t-statistics until sufficient number of studies obtained
  dat <- data.frame()
  
  while (nrow(dat) < studies) {
    # simulate true effects
    delta_i <- rnorm(studies, mean = mean_effect, sd = sd_effect)
    
    # simulate sample sizes
    n_i <- n_sim(studies, n_factor)
    df <- 2 * (n_i - 1)
    
    # simulate t-statistics and p-values
    t_i <- rnorm(n = studies, mean = delta_i * sqrt(n_i / 2)) / sqrt(rchisq(n = studies, df = df) / df)
    p_onesided <- pt(t_i, df = df, lower.tail = FALSE)
    p_twosided <- 2 * pt(abs(t_i), df = df, lower.tail = FALSE)
    
    # effect censoring based on p-values
    p_observed <- c(1, p_RR)[cut(p_onesided, c(0, p_thresholds, 1), labels = FALSE, include.lowest = TRUE)]
    observed <- runif(studies) < p_observed
    
    # put it all together
    if (nrow(dat) + sum(observed) > studies) observed <- which(observed)[1:(studies - nrow(dat))]
    new_dat <- data.frame(n = n_i[observed], t = t_i[observed], p = p_twosided[observed])
    dat <- rbind(dat, new_dat)
  }
  
  # calculate standardized mean difference estimates (Hedges' g's) and transformed d's
  
  J_i <- with(dat, 1 - 3 / (8 * n - 9))
  
  dat <- within(dat, {
    
    d <- sqrt(2 / n) * t
    g <- J_i * d
    
    Vd <- 2 / n + g^2 / (4 * (n - 1))
    Vg <- J_i^2 * Vd
    sd <- sqrt(Vg)
    
    Va <- 2 / n
    sda <- sqrt(Va)
    
    Top <- (1:studies) %in% (order(n, decreasing = TRUE)[1:10])
  })
  
  return(dat)
}


#----------------------------------------
# Run all of the methods
#----------------------------------------

fit_meta_ML <- function(dat, max_iter = 100L, step_adj = 1L, tau2_min = NULL) {
  
  suppressPackageStartupMessages(
    require(metafor, quietly = TRUE, warn.conflicts = FALSE)
  )
  
  require(rlang, quietly = TRUE, warn.conflicts = FALSE)
  
  if (is.null(tau2_min)) tau2_min <- -min(dat$Va)
  
  dat <- enexpr(dat)
  
  rma_ML <- NULL
  fits <- 0L
  
  while (is.null(rma_ML) & fits <= 5) {
    
    control_list <- list(maxiter = max_iter, stepadj = step_adj, tau2.min = tau2_min)
    
    rma_call <- expr(rma(yi = g, vi = Va, data = !!dat, method = "ML", control = !!control_list))
    
    # rma_env <- child_env(caller_env(), control_list = control_list)
    
    rma_ML <- suppressWarnings(
      tryCatch(eval_bare(rma_call, caller_env()),
        error = function(e) NULL)
    )
    
    fits <- fits + 1L
    step_adj <- step_adj / 2L
  }
  
  rma_ML
}

fit_meta_FE_REML <- function(dat, max_iter = 100L, step_adj = 1L, tau2_min = -min(dat$Va)) {
  
  suppressPackageStartupMessages(
    require(metafor, quietly = TRUE, warn.conflicts = FALSE)
  )
  
  rma_ML <- NULL
  fits <- 0L
  
  while (is.null(rma_ML) & fits <= 5) {
    
    rma_ML <- suppressWarnings(
      tryCatch(
        rma(yi = g, vi = Va, weights = 1 / Va, data = dat, method = "REML", 
            control = list(maxiter = max_iter, stepadj = step_adj, tau2.min = tau2_min)),
        error = function(e) NULL)
    )
    
    fits <- fits + 1L
    step_adj <- step_adj / 2L
  }
  
  rma_ML
}

estimate_effects <- function(dat, 
                             test_steps = .025, 
                             score_test_types = NULL,
                             LRT_types = NULL,
                             boot_n_sig = FALSE,
                             boot_qscore = FALSE,
                             max_iter = 100L,
                             step_adj = 1L,
                             tau2_min = -min(dat$Va)) {
  
  suppressPackageStartupMessages(require(purrrlyr))
  
  rma_ML <- fit_meta_ML(dat, max_iter = max_iter, step_adj = step_adj, tau2_min = tau2_min)
  rma_FE_REML <- fit_meta_FE_REML(dat, max_iter = max_iter, step_adj = step_adj, tau2_min = tau2_min)
  
  mods <- list(ML = rma_ML, `FE-REML` = rma_FE_REML)
  
  n_nonsig <- nrow(dat) - n_sig(yi = dat$d, sei = dat$sda, step = test_steps)
  
  res <- data_frame()
  
  if (!is.null(score_test_types)) {
    res_score <- 
      score_test_types %>%
      invoke_rows(VHSM_score_test, .d = ., model = rma_ML, steps = test_steps, .to = "test_stats") %>%
      unnest(test_stats) %>%
      rename(Stat = Score_stat)
    
    res <- bind_rows(res, res_score)
  }
  
  if (!is.null(LRT_types)) {
    res_LRT <-
      LRT_types %>%
      invoke_rows(LRT_VHSM, .d = ., model = rma_ML, steps = test_steps, .to = "test_stats") %>%
      unnest(test_stats) %>%
      select(-edge_RE, -edge_VHSM, -RE, -VHSM) %>%
      rename(Stat = LRT) %>%
      mutate(
        type = "LRT"
      )
  
    res <- bind_rows(res, res_LRT)    
  }

  if (boot_n_sig) {
    res_n_sig <-
      map_dfr(mods, bootstrap_n_sig, step = test_steps, .id = "info") %>%
      mutate(two_sided = FALSE, type = "n-sig bootstrap", df = length(test_steps))
    
    res <- bind_rows(res, res_n_sig)
  }
  
  if (boot_qscore) {
    res_qscore <-
      map_dfr(mods, bootstrap_quick_score, steps = test_steps, .id = "info") %>%
      mutate(two_sided = TRUE, type = "quick score bootstrap", df = length(test_steps))
    
    res <- bind_rows(res, res_qscore)
  }
  
  res %>%
    mutate(non_sig = n_nonsig)
  
}


#------------------------------------------------------
# Simulation Driver
#------------------------------------------------------

runSim <- function(reps, 
                   studies, mean_effect, sd_effect, 
                   n_sim, n_factor = 1L, 
                   p_thresholds = .025, p_RR = 1,
                   test_steps = .025, 
                   score_test_types = NULL, 
                   LRT_types = NULL,
                   boot_n_sig = FALSE,
                   boot_qscore = FALSE,
                   seed = NULL, ...) {
  
  suppressPackageStartupMessages(require(purrr))
  suppressPackageStartupMessages(require(dplyr))
  suppressPackageStartupMessages(require(tidyr))
  
  if (!is.null(seed)) set.seed(seed)
  
  grouping_vars <- union(names(score_test_types), names(LRT_types))
  
  rerun(reps, {
    r_SMD(studies, mean_effect, sd_effect, n_sim, 
          p_thresholds = p_thresholds, p_RR = p_RR) %>%
      estimate_effects(test_steps = test_steps, 
                       score_test_types = score_test_types, 
                       LRT_types = LRT_types, 
                       boot_n_sig = boot_n_sig,
                       boot_qscore = boot_qscore)  
  }) %>%
    bind_rows() %>%
    group_by_at(.vars = grouping_vars) %>% 
    summarise(
      pct_all_sig = mean(non_sig == 0),
      pct_NA = mean(is.na(p_val)),
      reject_025 = mean(p_val < .025, na.rm = TRUE),
      reject_050 = mean(p_val < .050, na.rm = TRUE),
      reject_100 = mean(p_val < .100, na.rm = TRUE)
    )
}
