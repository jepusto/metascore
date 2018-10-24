
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

#-------------------------------------------
# 3-parameter selection model
#-------------------------------------------

fit_3PSM <- function(g, Vg, V_lab) {
  
  require(weightr, quietly = TRUE, warn.conflicts = FALSE)
  
  # count sig and non-sig p-values
  p_vals <- pnorm(g / sqrt(Vg), lower.tail = FALSE)
  p_counts <- mean(p_vals < .025)
  k <- length(g)
  
  # if too few significant p-values, return unadjusted estimates, set p(no missing) = 1
  
  if (p_counts < 2 / k) {
    rma_fit <- tryCatch(rma(yi = g, vi = Vg), error = function(e) NULL)
    if (is.null(rma_fit)) rma_fit <- rma(yi = g, vi = Vg, method = "FE")
    res <- data.frame(
      V = V_lab, estimator = "3PSM", se_method = "",
      est = rma_fit$b[[1]], 
      se = rma_fit$se, 
      p_val = rma_fit$pval,
      p_nomissing = 1
    )  
    return(res)
  }
  
  # if all significant p-values, set selection weight for non-sig p-values to 1 / (k + 1)
  selection_wts <- if (p_counts==1) c(1, 1 / (k + 1)) else NULL
  
  # fit selection model
  wf <- suppressWarnings(
    tryCatch(
      weightfunct(g, Vg, steps = c(.025, 1), weights = selection_wts), 
      error = function(e) NULL)
  )
  
  # return NA if errors
  if (is.null(wf)) {
    res <- data.frame(
      V = V_lab, estimator = "3PSM", se_method = "",
      est = NA, se = NA, p_val = NA,
      p_nomissing = NA
    )
    return(res)
  } 
  
  if (p_counts==1) {
    
    # if all significant effects
    chol_hess <- tryCatch(chol(wf[[2]]$hessian[1:2,1:2]), error = function(e) NULL)
    
    if (is.null(chol_hess)) {
      wf <- suppressWarnings(weightfunct(g, Vg, steps = c(.025, 1), weights = selection_wts, fe = TRUE))
      est <- wf[[2]]$par[1]
      SE <- sqrt(1 / wf[[2]]$hessian[1,1])[1]
    } else {
      est <- wf[[2]]$par[2]
      SE <- sqrt(diag(chol2inv(chol_hess)))[2]
    }
    
  } else {
    
    # if not all significant effects
    chol_hess <- tryCatch(chol(wf[[2]]$hessian), error = function(e) NULL)
    
    if (is.null(chol_hess)) {
      if (wf[[2]]$par[[1]] < 10^-4) {
        # if tau estimate is near zero, move to FE model
        wf <- suppressWarnings(weightfunct(g, Vg, steps = c(.025, 1), weights = selection_wts, fe = TRUE))
        est <- wf[[2]]$par[1]
        chol_hess <- tryCatch(chol(wf[[2]]$hessian), error = function(e) NULL)
        SE <- if (is.null(chol_hess)) NA else sqrt(diag(chol2inv(chol_hess)))[1]
      } else {
        # otherwise give up
        est <- wf[[2]]$par[2]
        SE <- Inf
      }
    } else {
      est <- wf[[2]]$par[2]
      SE <- sqrt(diag(chol2inv(chol_hess)))[2]
    }
  }
  
  if (p_counts==1) {
    LR_pval <- NA
  } else {
    lrchisq <- 2 * (abs(wf[[1]]$value - wf[[2]]$value))
    LR_pval <- 1 - pchisq(lrchisq, 1L)
  }
  
  res <- data.frame(
    V = V_lab,
    estimator = "3PSM",
    se_method = "",
    est = est,
    se = SE,
    p_val = (2 * pnorm(-abs(est / SE))),
    p_nomissing = LR_pval
  )
  
  return(res)
  
}


#----------------------------------------
# Run all of the methods
#----------------------------------------

fit_meta_ML <- function(dat, max_iter = 100L, step_adj = 1L, tau2_min = -min(dat$Va)) {
  
  suppressPackageStartupMessages(
    require(metafor, quietly = TRUE, warn.conflicts = FALSE)
  )
  
  rma_ML <- NULL
  fits <- 0L
  
  while (is.null(rma_ML) & fits <= 5) {
    
    rma_ML <- suppressWarnings(
      tryCatch(
        rma(yi = g, vi = Va, data = dat, method = "ML", 
            control = list(maxiter = max_iter, stepadj = step_adj, tau2.min = tau2_min)),
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
                             boot_n_sig = FALSE,
                             boot_qscore = FALSE,
                             max_iter = 100L,
                             step_adj = 1L,
                             tau2_min = -min(dat$Va)) {
  
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
      rename(Stat = Q_score)
    
    res <- bind_rows(res, res_score)
  }
  
  if (boot_n_sig) {
    res_n_sig <-
      map_dfr(mods, bootstrap_n_sig, step = test_steps, .id = "info") %>%
      mutate(type = "n-sig bootstrap", prior_mass = NA, df = length(test_steps))
    
    res <- bind_rows(res, res_n_sig)
  }
  
  if (boot_qscore) {
    res_qscore <-
      map_dfr(mods, bootstrap_quick_score, steps = test_steps, .id = "info") %>%
      mutate(type = "quick score bootstrap", prior_mass = NA, df = length(test_steps))
    
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
                   boot_n_sig = FALSE,
                   boot_qscore = FALSE,
                   seed = NULL, ...) {
  
  suppressPackageStartupMessages(require(purrr))
  suppressPackageStartupMessages(require(purrrlyr))
  suppressPackageStartupMessages(require(dplyr))
  suppressPackageStartupMessages(require(tidyr))
  
  if (!is.null(seed)) set.seed(seed)
  
  rerun(reps, {
    r_SMD(studies, mean_effect, sd_effect, n_sim, 
          p_thresholds = p_thresholds, p_RR = p_RR) %>%
      estimate_effects(test_steps = test_steps, 
                       score_test_types = score_test_types, 
                       boot_n_sig = boot_n_sig,
                       boot_qscore = boot_qscore)  
  }) %>%
    bind_rows() %>%
    group_by(type, info, prior_mass) %>% 
    summarise(
      pct_all_sig = mean(non_sig == 0),
      pct_NA = mean(is.na(p_val)),
      reject_025 = mean(p_val < .025, na.rm = TRUE),
      reject_050 = mean(p_val < .050, na.rm = TRUE),
      reject_100 = mean(p_val < .100, na.rm = TRUE)
    )
}
