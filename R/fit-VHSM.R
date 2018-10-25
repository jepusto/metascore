#-------------------------------------------
# fit Vevea-Hedges selection model
#-------------------------------------------

fit_VHSM <- function(model, steps = .025, tol = 10^-3, method = "L-BFGS-B", use_gradient = TRUE, control = list()) {
  
  y <- as.vector(model$yi)
  s <- sqrt(model$vi)
  X <- model$X
  k <- model$k
  p <- NCOL(X)
  q <- length(steps)
  
  if (is.null(steps)) {
    par <- c(model$tau2, as.vector(model$b))
    lower <- c((tol - 1) * min(s)^2, rep(-Inf, p))  
  } else {
    par <- c(model$tau2, as.vector(model$b), rep(1, q))
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

LRT_VHSM <- function(model, steps = .025, k_min = 2, method = "L-BFGS-B", control = list()) {
  
  
  # count sig and non-sig p-values
  y <- as.vector(model$yi)
  s <- sqrt(model$vi)
  p_vals <- pnorm(y / s, lower.tail = FALSE)
  cats <- cut(p_vals, breaks = c(0, steps, 1), include.lowest = TRUE)
  ns_count <- table(cats)
  
  # if too many significant p-values, adjust step
  new_step <- if (ns_count < k_min) mean(p_vals[rank(p_vals) %in% (k - k_min - 0:1)]) else step
  
  # fit random effects model 
  
  null_opt <- refit_RE(model = model, method = method, control = control)
  
  # fit selection model
  
  VHSM_opt <- fit_VHSM(model = model, steps = steps, method = method, control = control)
  
  # Likelihood ratio test
  
  LRT <- 2 * (null_opt$value - VHSM_opt$value)
  df <- length(steps)
  p_val <- pchisq(LRT, df = df, lower.tail = FALSE)
  
  list(
    RE = null_opt, 
    `3PSM` = VHSM_opt,
    LRT = LRT,
    df = df,
    p_val = p_val
  )
  
}
