#---------------------------------------
# Vevea-Hedges selection model 
# log likelihood 
#---------------------------------------

VHSM_loglik <- function(beta, tau_sq, omega, steps, y, s, X = matrix(rep(1, length(y)))) {
 
  omega_vec <- c(1, omega)
  mu_vec <- as.vector(X %*% beta)
  eta <- sqrt(tau_sq + s^2)
  
  c_mat <- (tcrossprod(s, -qnorm(steps)) - mu_vec) / eta
  N_mat <- pnorm(c_mat)
  B_mat <- cbind(1, N_mat) - cbind(N_mat, 0)
  
  p_vals <- 1 - pnorm(y / s)
  cats <- cut(p_vals, breaks = c(0, steps, 1), include.lowest = TRUE)
  weight_vec <- omega_vec[cats]
  
  Ai <- as.vector(B_mat %*% omega_vec)
  
  sum(log(weight_vec)) - sum(((y - mu_vec) / eta)^2) / 2 - sum(log(eta)) - sum(log(Ai))
}

VHSM_negloglik_theta <- function(theta, steps, y, s, X = matrix(rep(1, length(y)))) {

  p <- NCOL(X)
  
  if (is.null(steps)) {
    steps <- 0.5
    omega <- 1L
  } else {
    q <- length(steps)
    omega <- theta[p + 1 + 1:q]
  }

  -1 * VHSM_loglik(
    beta = theta[1:p], 
    tau_sq = theta[p+1], 
    omega = omega,
    steps = steps,
    y = y,
    s = s,
    X = X)
} 

#---------------------------------------
# Vevea-Hedges selection model 
# score and Hessian
#---------------------------------------

VHSM_prep <- function(beta, tau_sq, omega, steps, y, s, X) {
  
  omega_vec <- c(1, omega)
  mu_vec <- as.vector(X %*% beta)
  er <- y - mu_vec
  eta <- sqrt(tau_sq + s^2)
  
  c_mat <- (tcrossprod(s, -qnorm(steps)) - mu_vec) / eta
  
  N_mat <- pnorm(c_mat)
  B_mat <- cbind(1, N_mat) - cbind(N_mat, 0)
  Ai <- as.vector(B_mat %*% omega_vec)
  
  p_vals <- 1 - pnorm(y / s)
  cats <- cut(p_vals, breaks = c(0, steps, 1), include.lowest = TRUE)
  n_s <- as.vector(table(cats))
  
  d1_mat <- dnorm(c_mat)
  d2_mat <- c_mat * d1_mat
  
  dB_dmu <- (cbind(d1_mat, 0) - cbind(0, d1_mat)) / eta
  dA_dmu <- as.vector(dB_dmu %*% omega_vec)
  
  dB_dtausq <- (cbind(d2_mat, 0) - cbind(0, d2_mat)) / (2 * eta^2)
  dA_dtausq <- as.vector(dB_dtausq %*% omega_vec)
  
  dl_dbeta <- colSums((er / eta^2 - dA_dmu / Ai) * X)
  dl_dtausq <- sum(er^2 / eta^4) / 2 - sum(1 / eta^2) / 2 - sum(dA_dtausq / Ai)
  dl_domega <- n_s[-1] / omega - colSums(B_mat[,-1,drop=FALSE] / Ai)
  
  list(
    omega_vec = omega_vec,
    mu_vec = mu_vec,
    er = er,
    eta = eta,
    c_mat = c_mat,
    B_mat = B_mat,
    Ai = Ai,
    n_s = n_s,
    d1_mat = d1_mat,
    d2_mat = d2_mat,
    dB_dmu = dB_dmu,
    dB_dtausq = dB_dtausq,
    dA_dmu = dA_dmu,
    dA_dtausq = dA_dtausq,
    dl_dbeta = dl_dbeta,
    dl_dtausq = dl_dtausq,
    dl_domega = dl_domega
  )
}

VHSM_score <- function(beta, tau_sq, omega, steps, y, s, X, prep = NULL) {
  
  if (is.null(prep)) prep <- VHSM_prep(beta, tau_sq, omega, steps, y, s, X)
  
  with(prep, c(dl_dbeta, dl_dtausq, dl_domega))
}

VHSM_neg_score_theta <- function(theta, steps, y, s, X) {
  
  p <- NCOL(X)
  
  if (is.null(steps)) {
    
    -1 * null_score(
      beta = theta[1:p],
      tau_sq = theta[p+1],
      steps = 0.5, 
      y = y,
      s = s,
      X = X)[1:(p+1)]
    
  } else {
    
    q <- length(steps)
    
    -1 * VHSM_score(
      beta = theta[1:p], 
      tau_sq = theta[p+1], 
      omega = theta[p + 1 + 1:q],
      steps = steps,
      y = y,
      s = s,
      X = X)
  }
  
} 

VHSM_Info <- function(beta, tau_sq, omega, steps, y, s, X, prep = NULL) {
  
  if (is.null(prep)) prep <- VHSM_prep(beta, tau_sq, omega, steps, y, s, X)
  
  d3_mat <- with(prep, (c_mat^2 - 3) * d2_mat)
  d4_mat <- with(prep, (c_mat^2 - 1) * d1_mat)
  
  d2A_dmu_dmu <- 2 * prep$dA_dtausq
  d2A_dtausq_dtausq <- as.vector((cbind(d3_mat, 0) - cbind(0, d3_mat)) %*% prep$omega_vec) / (4 * prep$eta^4)
  d2A_dmu_dtausq <- as.vector((cbind(d4_mat, 0) - cbind(0, d4_mat)) %*% prep$omega_vec) / (2 * prep$eta^3)

  d2l_dbeta_dbeta <- with(prep, crossprod(X, (1 / eta^2 + d2A_dmu_dmu / Ai - (dA_dmu / Ai)^2) * X))
  d2l_dbeta_dtausq <- with(prep, colSums((d2A_dmu_dtausq / Ai + er / eta^4 - dA_dmu * dA_dtausq / Ai^2) * X))
  d2l_dbeta_domega <- with(prep, crossprod(X, (dB_dmu[,-1, drop = FALSE] - (dA_dmu / Ai) * B_mat[,-1,drop=FALSE]) / Ai))
  d2l_dtausq_dtausq <- with(prep, sum(er^2 / eta^6) + sum(d2A_dtausq_dtausq / Ai) - sum(1 / eta^4) / 2 - sum((dA_dtausq / Ai)^2))
  d2l_dtausq_domega <- with(prep, colSums(dB_dtausq[,-1] / Ai - (dA_dtausq / Ai^2) * B_mat[,-1]))
  d2l_domega_domega <- with(prep, diag(n_s[-1] / omega^2, nrow = length(n_s) - 1L) - crossprod(B_mat[,-1] / Ai))
  
  I_mat <- rbind(
    cbind(d2l_dbeta_dbeta, d2l_dbeta_dtausq, d2l_dbeta_domega),
    cbind(t(d2l_dbeta_dtausq), d2l_dtausq_dtausq, t(d2l_dtausq_domega)),
    cbind(t(d2l_dbeta_domega), d2l_dtausq_domega, d2l_domega_domega)
  )
  
  dimnames(I_mat) <- NULL
  
  I_mat
}

#---------------------------------------
# Vevea-Hedges selection model 
# score and Hessian 
# under null of no selection
#---------------------------------------

null_prep <- function(beta, tau_sq, steps, y, s, X) {
  
  mu_vec <- as.vector(X %*% beta)
  er <- y - mu_vec
  eta <- sqrt(tau_sq + s^2)
  
  c_mat <- (tcrossprod(s, -qnorm(steps)) - mu_vec) / eta
  
  N_mat <- pnorm(c_mat)
  B_mat <- cbind(1, N_mat) - cbind(N_mat, 0)
  
  p_vals <- 1 - pnorm(y / s)
  cats <- cut(p_vals, breaks = c(0, steps, 1), include.lowest = TRUE)
  n_s <- as.vector(table(cats))
  
  list(
    mu_vec = mu_vec,
    er = er,
    eta = eta,
    c_mat = c_mat,
    B_mat = B_mat,
    cats = cats,
    n_s = n_s
  )
}

null_score_matrix <- function(beta, tau_sq, steps, y, s, X, prep = NULL) {
  
  if (is.null(prep)) prep <- null_prep(beta, tau_sq, steps, y, s, X)
  
  obs_matrix <- model.matrix(~ 0 + prep$cats)
    
  dl_dbeta <- with(prep, (er / eta^2) * X)
  dl_dtausq <- with(prep, er^2 / eta^4 - 1 / eta^2) / 2
  dl_domega <- with(prep, obs_matrix[,-1] - B_mat[,-1])
  
  cbind(dl_dbeta, dl_dtausq, dl_domega)
}

null_score <- function(beta, tau_sq, steps, y, s, X, prep = NULL) {
  
  if (is.null(prep)) prep <- null_prep(beta, tau_sq, steps, y, s, X)
  
  dl_dbeta <- with(prep, colSums((er / eta^2) * X))
  dl_dtausq <- with(prep, sum(er^2 / eta^4) / 2 - sum(1 / eta^2) / 2)
  dl_domega <- with(prep, n_s[-1] - colSums(B_mat[,-1,drop=FALSE]))
  
  c(dl_dbeta, dl_dtausq, dl_domega)
}

null_Info <- function(beta, tau_sq, steps, y, s, X, prep = NULL, info = "observed") {
  
  if (is.null(prep)) prep <- null_prep(beta, tau_sq, steps, y, s, X)
  
  d1_mat <- dnorm(prep$c_mat)
  d2_mat <- prep$c_mat * d1_mat
  
  dB_dmu <- (cbind(d1_mat, 0) - cbind(0, d1_mat)) / prep$eta
  dB_dtausq <- (cbind(d2_mat, 0) - cbind(0, d2_mat)) / (2 * prep$eta^2)
  
  d2l_dbeta_dbeta <- with(prep, crossprod(X / eta))
  d2l_dbeta_domega <- crossprod(X, dB_dmu[,-1])
  d2l_dtausq_domega <- colSums(dB_dtausq[,-1,drop=FALSE])
  
  if (info == "observed") {
    d2l_dbeta_dtausq <- with(prep, colSums((er / eta^4) * X))
    d2l_dtausq_dtausq <- with(prep, sum(er^2 / eta^6) - sum(1 / eta^4) / 2)
    d2l_domega_domega <- with(prep, diag(n_s[-1], nrow = length(n_s) - 1L) - crossprod(B_mat[,-1]))
  } else {
    d2l_dbeta_dtausq <- rep(0, NCOL(X))
    d2l_dtausq_dtausq <- sum(1 / prep$eta^4) / 2
    d2l_domega_domega <- with(prep, diag(colSums(B_mat[,-1,drop=FALSE]), nrow = length(n_s) - 1L) - crossprod(B_mat[,-1]))
  }
  
  I_mat <- rbind(
    cbind(d2l_dbeta_dbeta, d2l_dbeta_dtausq, d2l_dbeta_domega),
    cbind(t(d2l_dbeta_dtausq), d2l_dtausq_dtausq, t(d2l_dtausq_domega)),
    cbind(t(d2l_dbeta_domega), d2l_dtausq_domega, d2l_domega_domega)
  )
  dimnames(I_mat) <- NULL  
  
  I_mat  
}
