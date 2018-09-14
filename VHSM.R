N <- 100
p <- 3
X <- cbind(rep(1, N), matrix(rnorm(N * (p - 1)), N, p - 1))
beta <- rnorm(p, mean = 0.2, sd = 0.1)
tau_sq <- 0.1^2
omega <- c(0.5, 0.4)
steps <- c(0.025, 0.5)
s <- 1 / rchisq(N, df = 5)
y <- as.vector(X %*% beta) + rnorm(N, sd = sqrt(tau_sq)) + rnorm(N, sd = s)


selection_weight <- function(y, s, a, omega) {
  p_vals <- 1 - pnorm(y / s)
  cats <- cut(p_vals, breaks = a, include.lowest = TRUE)
  weight_vec <- omega[cats]
  cj <- as.vector(table(cats))
  list(weight_vec = weight_vec, cj = cj)
}

B_matrix <- function(mu_vec, steps, s, eta) {
  b_mat <- tcrossprod(s, -qnorm(steps)) - mu_vec
  N_mat <- pnorm(b_mat / eta)
  cbind(1, N_mat) - cbind(N_mat, 0)
}
  
VHSM_loglik <- function(beta, tau_sq, omega, steps, y, s, X) {
 
  a <- c(0, steps, 1)
  omega_vec <- c(1, omega)
  mu_vec <- as.vector(X %*% beta)
  eta <- sqrt(tau_sq + s^2)
  
  weight_vec <- selection_weight(y, s, a = a, omega = omega_vec)$weight_vec
  B_ij <- B_matrix(mu_vec = mu_vec, steps = steps, s = s, eta = eta)
  Ai <- as.vector(B_ij %*% omega_vec)
  
  sum(log(weight_vec)) - sum((y - mu_vec)^2 / V) / 2 - sum(log(V)) / 2 - sum(log(Ai))
  
}

VHSM_loglik(beta, tau_sq, omega, steps, y, s, X)

VHSM_score <- function(beta, tau_sq, omega, steps, y, s, X) {

  a <- c(0, steps, 1)
  omega_vec <- c(1, omega)
  mu_vec <- as.vector(X %*% beta)
  er <- y - mu_vec
  eta <- sqrt(tau_sq + s^2)
  
  N_mat <- pnorm(b_mat / eta)
  B_mat <- cbind(1, N_mat) - cbind(N_mat, 0)
  Ai <- as.vector(B_mat %*% omega_vec)
  
  s_weights <- selection_weight(y, s, a = a, omega = omega_vec)
  weight_vec <- s_weights$weight_vec
  cj <- s_weights$cj
  
  b_mat <- tcrossprod(s, -qnorm(steps)) - mu_vec
  n_mat <- dnorm(b_mat / eta)
  dA_dbeta <- X * (as.vector((cbind(n_mat, 0) - cbind(0, n_mat)) %*% omega_vec) / eta)

  t_mat <- b_mat * n_mat
  dA_dtau_sq <- (as.vector((cbind(t_mat, 0) - cbind(0, t_mat)) %*% omega_vec) / eta^3)

  dl_dbeta <- colSums(X * er / eta^2) - colSums(dA_dbeta / Ai)
  dl_dtau_sq <- sum(er^2 / eta^4) / 2 - sum(1 / eta^2) / 2 - sum(dA_dtau_sq / Ai)
  dl_domega <- cj[-1] / omega - colSums(B_mat[,-1] / Ai)
}

VHSM_Hessian <- function(beta, tau_sq, omega, steps, y, s, X) {
  
}
