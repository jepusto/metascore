r_rademacher <- function(n) sample(c(-1L, 1L), size = n, replace = TRUE)

r_mammen <- function(n) {
  sample(
    (1 + c(-1, 1) * sqrt(5)) / 2, 
    size = n, replace = TRUE, 
    prob = (sqrt(5) + 1) / (2 * sqrt(5)) * c(1, -1) + c(0, 1)
  )
}

r_webb <- function(n) sample(sqrt(c(-3, -2, -1, 1, 2, 3) / 2), size = n, replace = TRUE)
  
# stat <- function(x) VHSM_score_test(x, steps = .025, type = "robust")
# reps <- 4999L
# r_func <- r_rademacher

wild_bootstrap <- function(
  model, 
  stat, ..., 
  update = TRUE,
  reps = 999L, 
  r_func = r_rademacher, 
  plot = FALSE,
  return_boot_dist = FALSE,
  seed = NULL
) {
  
  if (!is.null(seed)) set.seed(seed)
  
  x_hat <- fitted(model)
  res <- residuals(model)
  k <- model$k
  
  S <- stat(model, ...)

  if (update) {
    
    booties <- purrr::rerun(.n = reps, {
      z <- r_func(k)
      yi <- x_hat + res * z
        mod_b <- update(model, yi = yi)
        stat(mod_b, ...)
    })
  
  } else {
    
    booties <- purrr::rerun(.n = reps, {
      z <- r_func(k)
      model$yi <- x_hat + res * z
      stat(model, ...)
    })
    
  } 
    
  booties <- unlist(booties)
  
  if (plot) {
    print(plot(density(booties), main = "Bootstrap distribution"))
    abline(v = S)
  }
  
  if (return_boot_dist) return(booties)
  
  p_val <- mean(S <= booties)
  
  data.frame(Stat = S, p_val = p_val)
  
}

n_sig <- function(yi, sei, df = Inf, step = .025) {
  p_onesided <- pt(yi / sei, df = df, lower.tail = FALSE)
  sum(p_onesided < step)
}

bootstrap_n_sig <- function(
  model, 
  step = .025,
  reps = 999L, 
  r_func = r_rademacher, 
  return_boot_dist = FALSE,
  seed = NULL
) {
  
  if (!("rma.uni" %in% class(model))) {
    return(data.frame(Stat = NA, p_val = NA))
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  x_hat <- fitted(model)
  res <- residuals(model)
  sei <- sqrt(model$vi)
  k <- model$k
  
  S <- n_sig(model$yi, sei = sei, step = step)
  
  booties <- purrr::rerun(.n = reps, {
    z <- r_func(k)
    yi <- x_hat + res * z
    n_sig(yi = yi, sei = sei, step = step)
  })
  
  booties <- unlist(booties)
  
  if (return_boot_dist) return(booties)
  
  p_val <- mean(S <= booties)
  
  data.frame(Stat = S, p_val = p_val)
  
}

bootstrap_quick_score <- function(
  model, 
  steps = .025,
  reps = 999L, 
  r_func = r_rademacher, 
  return_boot_dist = FALSE,
  seed = NULL
) {

  if (!("rma.uni" %in% class(model))) {
    return(data.frame(Stat = NA, p_val = NA))
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  
  beta <- as.vector(model$beta)
  tau_sq <- model$tau2
  y <- as.vector(model$yi)
  s <- sqrt(model$vi)
  X <- model$X
  k <- model$k
  q <- length(steps)
  
  x_hat <- fitted(model)
  res <- residuals(model)
  
  S <- quick_score_Q(beta = beta, tau_sq = tau_sq, steps = steps, y = y, s = s, X = X, q = q)
  
  booties <- purrr::rerun(.n = reps, {
    z <- r_func(k)
    y_boot <- x_hat + res * z
    quick_score_Q(beta = beta, tau_sq = tau_sq, steps = steps, y = y_boot, s = s, X = X, q = q)
  })
  
  booties <- unlist(booties)
  
  if (return_boot_dist) return(booties)
  
  p_val <- mean(S <= booties)
  
  data.frame(Stat = S, p_val = p_val)
  
}
