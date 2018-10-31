

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
    upper <- rep(Inf, p + 1)
  } else {
    par <- c(tau_sq_start, beta_start, rep(1, q))
    lower <- c((tol - 1) * min(s)^2, rep(-Inf, p), rep(tol, q))  
    upper <- c(rep(Inf, p + 1), rep(1 / tol, q))
  }
  
  optim_args <- list(par = par, fn = VHSM_negloglik_theta, 
                     steps = steps, 
                     y = y, s = s, X = X, 
                     method = method,
                     control = control)
  
  if (use_gradient) optim_args$gr <- VHSM_neg_score_theta
  
  if (method == "L-BFGS-B") {
    optim_args$lower <- lower
    optim_args$upper <- upper
  } 
  
  res <- do.call(optim, args = optim_args)

  res$edge <- any(res$par == lower | res$par == upper)
  
  res
    
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
  p_ordered <- sort(p_vals)
  
  J <- length(steps)
  
  spread_to <- which.max(n_count - cumsum(pmax(k_min - n_count, 0)) > k_min)
  
  if (spread_to > 1) {
    for (j in 1:(spread_to - 1)) {
      new_steps[j] <- mean(p_ordered[p_ordered > steps[j]][0:1 + k_min - n_count[[j]]])
      n_count[j + 1] <- n_count[[j + 1]] - k_min + n_count[[j]]
      n_count[j] <- k_min
    }
  }

  for (j in J:spread_to) {
    if (n_count[[1 + j]] < k_min) {
      new_steps[j] <- mean(rev(p_ordered[p_ordered < steps[j]])[0:1 + k_min - n_count[[1 + j]]])
      n_count[j] <- n_count[[j]] - k_min + n_count[[1 + j]]
      n_count[j + 1] <- k_min
    } 
  }
  
  new_steps
}

LRT_VHSM <- function(model, 
                     steps = .025, 
                     k_min = 3L, 
                     two_sided = TRUE,
                     tol = 10^-3, 
                     method = "L-BFGS-B", 
                     use_gradient = TRUE, 
                     control = list()) {
  
  if (length(steps) > 1 & !two_sided) stop("One-sided tests only available for models with a single step.")
  
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
  
  if (two_sided) {
    p_val <- pchisq(LRT, df = df, lower.tail = FALSE)

  } else {
    
    Z <- sign(VHSM_opt$par[length(VHSM_opt$par)] - 1) * sqrt(LRT)
    p_val <- pnorm(Z)
    
  }  
  
  tibble::data_frame(
    LRT = LRT,
    df = df,
    p_val = p_val,
    edge_RE = null_opt$edge,
    edge_VHSM = VHSM_opt$edge,
    RE = list(null_opt),
    VHSM = list(VHSM_opt)
  )
  
}
