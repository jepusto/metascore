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
    
  booties <- purrr::rerun(.n = reps, {
    z <- r_func(k)
    yi <- x_hat + res * z
    mod_b <- update(model, yi = yi)
    stat(mod_b, ...)
  })
  
  booties <- unlist(booties)
  
  if (plot) {
    print(plot(density(booties), main = "Bootstrap distribution"))
    abline(v = S)
  }
  
  if (return_boot_dist) return(booties)
  
  p_val <- mean(S <= booties)
  
  data.frame(Stat = S, p_val = p_val)
  
}
