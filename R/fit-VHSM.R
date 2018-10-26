

#-------------------------------------------
# fit Vevea-Hedges selection model
#-------------------------------------------

fit_VHSM <- function(y, s, X, 
                     steps = .025, 
                     tau_sq_start = mean(s^2) / 4, 
                     beta_start = mean(y),
                     tol = 10^-3, 
                     method = "L-BFGS-B", 
                     use_gradient = TRUE, 
                     control = list()) {
  
  k <- length(y)
  p <- NCOL(X)
  q <- length(steps)
  
  if (is.null(steps)) {
    par <- c(tau_sq_start, beta_start)
    lower <- c((tol - 1) * min(s)^2, rep(-Inf, p))  
  } else {
    par <- c(tau_sq_start, beta_start, rep(1, q))
    lower <- c((tol - 1) * min(s)^2, rep(-Inf, p), rep(tol, q))  
  }
  
  optim_args <- list(par = par, fn = VHSM_negloglik_theta, 
                     steps = steps, 
                     y = y, s = s, X = X, 
                     method = method,
                     control = control)
  
  if (use_gradient) optim_args$gr <- VHSM_neg_score_theta
  
  if (method == "L-BFGS-B") optim_args$lower <- lower
  
  do.call(optim, args = optim_args)

}

#-------------------------------------------
# Likelihood ratio test
#-------------------------------------------

p_cat <- function(p_vals, steps) {
  cats <- cut(p_vals, breaks = c(0, steps, 1), include.lowest = TRUE)
  table(cats)
}

find_new_steps <- function(p_vals, steps, k_min) {
  
  if (length(p_vals) < k_min * (length(steps) + 1)) stop("Number of effect sizes is too small. Use fewer steps or smaller k_min.")
  
  n_count <- p_cat(p_vals, steps)
  new_steps <- steps
  p_ordered <- sort(p_vals, decreasing = TRUE)
  
  for (j in length(steps):1) {
    if (n_count[[1 + j]] < k_min) {
      new_steps[j] <- mean(p_ordered[p_ordered < steps[j]][0:1 + k_min - n_count[[1 + j]]])
      n_count[j] <- n_count[[j]] - k_min + n_count[[1 + j]]
      n_count[j + 1] <- k_min
    } 
  }
  
  new_steps
}

LRT_VHSM <- function(model, 
                     steps = .025, 
                     k_min = 3L, 
                     tol = 10^-3, 
                     method = "L-BFGS-B", 
                     use_gradient = TRUE, 
                     control = list()) {
  
  
  # count sig and non-sig p-values
  y <- as.vector(model$yi)
  s <- sqrt(model$vi)
  X <- model$X
  
  # fit random effects model 
  
  null_opt <- fit_VHSM(y = y, s = s, X = X, 
                       steps = NULL,
                       tau_sq_start = model$tau2,
                       beta_start = as.vector(model$b),
                       tol = tol, method = method, control = control)

  # adjust steps
  
  p_vals <- pnorm(y / s, lower.tail = FALSE)
  new_steps <- find_new_steps(p_vals = p_vals, steps = steps, k_min = k_min)
  
  # fit selection model
  
  VHSM_opt <- fit_VHSM(y = y, s = s, X = X, 
                       steps = new_steps,
                       tau_sq_start = model$tau2,
                       beta_start = as.vector(model$b),
                       tol = tol, method = method, control = control)
  
  # Likelihood ratio test
  
  LRT <- 2 * (null_opt$value - VHSM_opt$value)
  df <- length(steps)
  p_val <- pchisq(LRT, df = df, lower.tail = FALSE)
  
  tibble::data_frame(
    LRT = LRT,
    df = df,
    p_val = p_val,
    RE = list(null_opt),
    VHSM = list(VHSM_opt)
  )
  
}
