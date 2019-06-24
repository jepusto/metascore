
score_CR <- function(mod, cluster, alpha = .025) {
  
  # pull required quantities from model
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
  
  # Fisher information entries
  FI_beta <- 1 / mod$se^2
  FI_tausq <- 1 / mod$se.tau2^2
  FI_pi <- sum(Pwr_i * (1 - Pwr_i))
  
  dnorm_c_i <- dnorm(c_i)
  FI_pi_beta <- -sum(dnorm_c_i * sqrt(inv_var_i))
  FI_pi_tausq <- -sum(c_i * dnorm_c_i * inv_var_i) / 2

  # Score contributions
  S_mu_i <- w_i * e_i
  
  S_tausq_i <- if (mod$method == "ML") {
    inv_var_i * (e_i^2 * inv_var_i - 1) / 2
  } else if (mod$method == "REML") { 
    inv_var_i * (e_i^2 * inv_var_i - 1 + inv_var_i / sum(inv_var_i)) / 2
  } else {
    0
  }
  
  # Efficient score
  F_i <- Pwr_i - O_i - S_mu_i * FI_pi_beta / FI_beta - S_tausq_i * FI_pi_tausq / FI_tausq
  
  # Clustered efficient score
  F_j <- tapply(F_i, cluster, sum)
  
  # Cluster-robust score statistic
  V_CR <- sum(F_j^2)
  Z_CR = Score_pi / sqrt(V_CR)
  
  data.frame(
    Z = Z_CR,
    p_val = pnorm(Z_CR)
  )

} 
