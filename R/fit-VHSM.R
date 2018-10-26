

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

LRT_VHSM <- function(model, 
                     steps = .025, 
                     k_min = 2, 
                     tol = 10^-3, 
                     method = "L-BFGS-B", 
                     use_gradient = TRUE, 
                     control = list()) {
  
  
  # count sig and non-sig p-values
  y <- as.vector(model$yi)
  s <- sqrt(model$vi)
  X <- model$X
  p_vals <- pnorm(y / s, lower.tail = FALSE)
  cats <- cut(p_vals, breaks = c(0, steps, 1), include.lowest = TRUE)
  ns_count <- table(cats)
  
  # if too many significant p-values, adjust step
  # new_step <- if (ns_count < k_min) mean(p_vals[rank(p_vals) %in% (k - k_min - 0:1)]) else step
  
  # fit random effects model 
  
  null_opt <- fit_VHSM(y = y, s = s, X = X, 
                       steps = NULL,
                       tau_sq_start = model$tau2,
                       beta_start = as.vector(model$b),
                       tol = tol, method = method, control = control)
  
  # fit selection model
  
  VHSM_opt <- fit_VHSM(y = y, s = s, X = X, 
                       steps = steps,
                       tau_sq_start = model$tau2,
                       beta_start = as.vector(model$b),
                       tol = tol, method = method, control = control)
  
  # Likelihood ratio test
  
  LRT <- 2 * (null_opt$value - VHSM_opt$value)
  df <- length(steps)
  p_val <- pchisq(LRT, df = df, lower.tail = FALSE)
  
  data.frame(
    LRT = LRT,
    df = df,
    p_val = p_val
  )
  
}
