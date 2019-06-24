
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
  
  # calculate standardized mean difference estimates (Hedges' g's)
  
  J_i <- with(dat, 1 - 3 / (8 * n - 9))
  
  dat <- within(dat, {
    
    d <- sqrt(2 / n) * t
    g <- J_i * d
    
    Vd <- 2 / n + g^2 / (4 * (n - 1))
    Vg <- J_i^2 * Vd
    sd <- sqrt(Vg)
    
    Va <- 2 / n
    sda <- sqrt(Va)
  })
  
  return(dat)
}

#----------------------------------------------------------------
# Calculate adjusted alpha
#----------------------------------------------------------------

tweak_alpha <- function(p_vals, alpha, k_min = 2L) {

  # adjust alpha
  sig <- p_vals < alpha
  n_sig <- sum(sig)
  
  if (n_sig < k_min) {
    p_ordered <- sort(p_vals)
    alpha <- mean(p_ordered[k_min + 0:1])
  } 
  
  if ((length(p_vals) - n_sig) < k_min) {
    p_ordered <- sort(p_vals, decreasing = TRUE)
    alpha <- mean(p_ordered[k_min + 0:1])
  }
  
  alpha
}

#----------------------------------------------------------------
# Likelihood ratio test for basic meta-analysis model with no covariates
#----------------------------------------------------------------

negloglik_3PSM <- function(theta, alpha, y, s, sig) {

  tau_sq <- theta[1]

  mu <- theta[2]

  if (is.null(alpha)) {
    alpha <- 0.5
    omega <- 1
  } else {
    omega <- theta[3]
  }

  weight_vec <- if_else(sig, 1, omega)
  
  eta_sq <- tau_sq + s^2
  
  c_i <- (s * qnorm(1 - alpha) - mu) / sqrt(eta_sq)
  A_i <- 1 - (1 - omega) * pnorm(c_i)
  
  log_lik <- sum(log(weight_vec)) - sum((y - mu)^2 / eta_sq) / 2 - sum(log(eta_sq)) / 2 - sum(log(A_i))
  
  -1 * log_lik
  
} 

fit_3PSM <- function(y, s, sig,
                     alpha = .025, 
                     tau_sq_start = mean(s^2) / 4, 
                     beta_start = mean(y),
                     tol = 10^-3, 
                     method = "L-BFGS-B", 
                     control = list()) {
  
  k <- length(y)

  if (is.null(alpha)) {
    par <- c(tau_sq_start, beta_start)
    lower <- c((tol - 1) * min(s)^2, -Inf)
    upper <- c(Inf, Inf)
  } else {
    par <- c(tau_sq_start, beta_start, 1)
    lower <- c((tol - 1) * min(s)^2, -Inf, tol)  
    upper <- c(Inf,Inf, 1 / tol)
  }
  
  optim_args <- list(par = par, fn = negloglik_3PSM, 
                     alpha = alpha, 
                     y = y, s = s, sig = sig,
                     method = method,
                     control = control)
  
  if (method == "L-BFGS-B") {
    optim_args$lower <- lower
    optim_args$upper <- upper
  } 
  
  do.call(optim, args = optim_args)
  
}

LRT_3PSM <- function(mod, 
                     alpha = .025, 
                     k_min = 2L, 
                     tol = 10^-3, 
                     method = "L-BFGS-B", 
                     use_gradient = TRUE, 
                     control = list()) {
  
  # count sig and non-sig p-values
  k <- mod$k
  y <- as.vector(mod$yi)
  s <- sqrt(mod$vi)
  p_vals <- pnorm(y / s, lower.tail = FALSE)
  n_sig <- sum(p_vals < alpha)
  
  # adjust alpha
  alpha <- tweak_alpha(p_vals, alpha = alpha, k_min = k_min) 
  sig <- p_vals < alpha
  
  # fit random effects model 
  
  null_opt <- fit_3PSM(y = y, s = s, sig = sig, 
                       alpha = NULL,
                       tau_sq_start = mod$tau2,
                       beta_start = as.vector(mod$b),
                       tol = tol, method = method, control = control)
  
  
  # fit selection model
  
  VHSM_opt <- fit_3PSM(y = y, s = s, sig = sig, 
                       alpha = alpha,
                       tau_sq_start = mod$tau2,
                       beta_start = as.vector(mod$b),
                       tol = tol, method = method, control = control)
  
  # Likelihood ratio test
  
  LRT <- 2 * (null_opt$value - VHSM_opt$value)

  Z_LRT <- sign(VHSM_opt$par[3] - 1) * sqrt(LRT)
  p_val <- pnorm(Z_LRT)
  
  tibble(
    model = "ML",
    type = "LRT",
    sig = n_sig / k,
    Z = Z_LRT,
    p_val = p_val
  )
}

#----------------------------------------------------------------
# Score tests for basic meta-analysis model with no covariates
#----------------------------------------------------------------

simple_scores <- function(mod, type = c("TES-norm","TES-binom","parametric","robust"), alpha = .025) {
  
  # pull required quantities from model
  k <- mod$k
  y_i <- as.numeric(mod$yi)
  e_i <- residuals(mod)
  v_i <- mod$vi
  s_i <- sqrt(v_i)
  mu <- as.numeric(mod$beta)
  tau_sq <- mod$tau2
  inv_var_i <- 1 / (v_i + tau_sq)
  w_i <- if (!is.null(mod$weights)) mod$weights else inv_var_i
  
  # calculate power per study  
  z_alpha <- qnorm(1 - alpha)
  c_i <- (s_i * z_alpha - mu) * sqrt(inv_var_i)
  Pwr_i <- pnorm(c_i, lower.tail = FALSE)
  
  # significance indicator
  O_i <- y_i / s_i > z_alpha
  
  Observed <- sum(O_i)
  Expected <- sum(Pwr_i)
  Score_pi <- Expected - Observed
  
  res <- tibble()
  
  # calculate TES
  
  if (any(c("TES-norm","TES-binom") %in% type)) {
    Z_TES <- Score_pi / sqrt(Expected * (k - Expected) / k)
    
    if ("TES-norm" %in% type) {
      res_TES <- tibble(
        type = "TES-norm",
        Z = Z_TES,
        p_val = pnorm(Z_TES)
      )
      res <- bind_rows(res, res_TES)    
    }
    
    if ("TES-binom" %in% type) {
      res_TES <- tibble(
        type = "TES-binom",
        Z = Z_TES,
        p_val = pbinom(k - Observed, size = k, prob = 1 - Expected / k)
      )
      res <- bind_rows(res, res_TES)
    }
  }
  
  if (any(c("parametric","robust") %in% type)) {
    
    FI_beta <- sum(w_i) # same as 1 / mod$se^2
    FI_tausq <- 1 / mod$se.tau2^2
    
    FI_pi <- sum(Pwr_i * (1 - Pwr_i))
    
    dnorm_c_i <- dnorm(c_i)
    FI_pi_beta <- -sum(dnorm_c_i * sqrt(inv_var_i))
    FI_pi_tausq <- -sum(c_i * dnorm_c_i * inv_var_i) / 2
    
    if ("parametric" %in% type) {
      V_parametric <- FI_pi - FI_pi_beta^2 / FI_beta - FI_pi_tausq^2 / FI_tausq
      Z_parametric = Score_pi / sqrt(V_parametric)
      
      res_parametric <- tibble(
        type = "parametric",
        Z = Z_parametric,
        p_val = pnorm(Z_parametric)
      )
      res <- bind_rows(res, res_parametric)
    }
    
    if ("robust" %in% type) {
      
      S_mu_i <- w_i * e_i
      S_tausq_i <- if (mod$method == "ML") {
        inv_var_i * (e_i^2 * inv_var_i - 1) / 2
      } else if (mod$method == "REML") { 
        inv_var_i * (e_i^2 * inv_var_i - 1 + inv_var_i / sum(inv_var_i)) / 2
      } else {
        0
      }
      F_i <- Pwr_i - O_i - S_mu_i * FI_pi_beta / FI_beta - S_tausq_i * FI_pi_tausq / FI_tausq
      
      V_robust <- sum(F_i^2)
      Z_robust = Score_pi / sqrt(V_robust)
      
      res_robust <- tibble(
        type = "robust",
        Z = Z_robust,
        p_val = pnorm(Z_robust)
      )
      
      res <- bind_rows(res, res_robust)
    }
  }
  
  res
} 


#----------------------------------------
# Run all of the methods
#----------------------------------------

fit_meta <- function(dat, method = "REML", weights = NULL, max_iter = 100L, step_adj = 1L, tau2_min = NULL) {
  
  suppressPackageStartupMessages(
    require(metafor, quietly = TRUE, warn.conflicts = FALSE)
  )
  
  require(rlang, quietly = TRUE, warn.conflicts = FALSE)
  
  if (is.null(tau2_min)) tau2_min <- -min(dat$Va)
  
  method <- enexpr(method)
  weights <- enexpr(weights)
  
  rma_fit <- NULL
  fits <- 0L
  
  while (is.null(rma_fit) & fits <= 5) {
    
    control_list <- list(maxiter = max_iter, stepadj = step_adj, tau2.min = tau2_min)
    
    rma_call <- 
      if (is.null(weights)) {
        expr(rma(yi = g, vi = Va, data = dat, method = !!method, control = !!control_list))  
      } else {
        expr(rma(yi = g, vi = Va, weights = !!weights, data = dat, method = !!method, control = !!control_list))
      }
    
    rma_fit <- suppressWarnings(
      tryCatch(eval_bare(rma_call, caller_env()), error = function(e) NULL)
    )
    
    fits <- fits + 1L
    step_adj <- step_adj / 2L
  }
  
  rma_fit
}

test_for_selection <- function(dat, alpha = .025, 
                               methods = c("FE","ML","REML","WLS"),
                               score_types = c("TES-norm","TES-binom","parametric","robust"),
                               tweak = c("no","yes"),
                               LRT = TRUE, k_min = 2L, tol = 10^-3, LRT_method = "L-BFGS-B"
                               ) {
  
  suppressPackageStartupMessages(
    require(dplyr, quietly = TRUE, warn.conflicts = FALSE)
  )
  suppressPackageStartupMessages(
    require(purrr, quietly = TRUE, warn.conflicts = FALSE)
  )

  fit_methods <- setdiff(methods, "WLS")
  
  if (length(fit_methods) > 0) {
    mods <- map(fit_methods, ~ fit_meta(dat = dat, method = .))
    names(mods) <- fit_methods
  } else {
    mods <- list()
  }
  
  if ("WLS" %in% methods) {
    mods <- c(mods, list(WLS = fit_meta(dat, method = "REML", weights = 1 / Va, tau2_min = 0)))
  }
  
  res <- tibble()
  
  if (!is.null(score_types) & "no" %in% tweak) {
    score_res <- map_dfr(mods, possibly(simple_scores, otherwise = tibble(type = score_types, Z = NA, p_val = NA)),
                         type = score_types, alpha = alpha, .id = "model")  
    res <- bind_rows(res, score_res)
  }
  
  if (!is.null(score_types) & "yes" %in% tweak) {
    alpha_tweak <- tweak_alpha(p_vals = pnorm(dat$g / dat$sda, lower.tail = FALSE), 
                               alpha = alpha, k_min = k_min) 
    
    tweak_res <- map_dfr(mods, possibly(simple_scores, otherwise = tibble(type = score_types, Z = NA, p_val = NA)),
                         type = score_types, alpha = alpha_tweak, .id = "model")  
    tweak_res$type <- paste(tweak_res$type, "tweaked", sep = "-")
    res <- bind_rows(res, tweak_res)
  }
  
  if (LRT & "ML" %in% methods) {
    possibly_LRT_3PSM <- possibly(LRT_3PSM, 
                                  otherwise = tibble(model = "ML", type = "LRT", 
                                                     sig = mean(pnorm(dat$g / dat$sda, lower.tail = FALSE) < alpha), 
                                                     Z = NA, p_val = NA))
    LRT_res <- possibly_LRT_3PSM(mod = mods$ML, alpha = alpha, k_min = k_min, tol = tol, method = LRT_method)
    res <- bind_rows(res, LRT_res)
  }
  
  res
  
} 


#------------------------------------------------------
# Simulation Driver
#------------------------------------------------------

runSim <- function(reps, 
                   studies, mean_effect, sd_effect, 
                   n_sim, n_factor = 1L, 
                   p_thresholds = .025, p_RR = 1,
                   test_alpha = .025, 
                   methods = c("FE","ML","REML","WLS"),
                   score_types = c("TES-norm","TES-binom","parametric","robust"),
                   tweak = c("no","yes"),
                   LRT = TRUE, k_min = 2L, tol = 10^-3, LRT_method = "L-BFGS-B",
                   seed = NULL, ...) {
  
  suppressPackageStartupMessages(require(purrr))
  suppressPackageStartupMessages(require(dplyr))
  suppressPackageStartupMessages(require(tidyr))
  
  if (!is.null(seed)) set.seed(seed)
 
  rerun(reps, {
    r_SMD(studies, mean_effect, sd_effect, n_sim, 
          p_thresholds = p_thresholds, p_RR = p_RR) %>%
      test_for_selection(alpha = test_alpha,
                         methods = methods,
                         score_types = score_types, tweak = tweak,
                         LRT = LRT, k_min = k_min, tol = tol, LRT_method = LRT_method)  
  }) %>%
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
}
