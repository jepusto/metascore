try_inverse <- function(x) tryCatch(chol2inv(chol(x)), error = function(e) NULL)

VHSM_score_test <- function(
  model, 
  steps, 
  type = "parametric", 
  info = "expected",
  prior_mass = 0L,
  return = "p-val",
  bootstrap = FALSE,
  ...
) {
  
  if (!("rma.uni" %in% class(model))) {
    return(data.frame(non_sig = NA, Q_score = NA, df = NA, p_val = NA))
  }

  beta <- as.vector(model$beta)
  tau_sq <- model$tau2
  y <- as.vector(model$yi)
  s <- sqrt(model$vi)
  X <- model$X
  k <- model$k
  
  prep <- null_prep(beta, tau_sq, steps, y, s, X)
  q <- length(steps)
  
  if (bootstrap) {
    BS <- wild_bootstrap(model, 
                         stat = VHSM_score_test, 
                         steps = steps, type = type, info = info,
                         return = "Q", 
                         ...)
    res <- data.frame(
      Q_score = BS$Stat, 
      df = q, 
      p_val = BS$p_val
    ) 
    return(res)
  }
  
  I_mat <- null_Info(beta, tau_sq, steps, y, s, X, prep = prep, info = info) / k
  
  if (type == "parametric") { 
    
    S_vec <- null_score(beta, tau_sq, steps, y, s, X, prep = prep)
    
    I_mat_inv <- try_inverse(I_mat)
    
    Q <- if (is.null(I_mat_inv)) NA else sum(I_mat_inv * tcrossprod(S_vec)) / k
    
  } else if (type == "subscore") {
    
    S_vec <- null_score(beta, tau_sq, steps, y, s, X, prep = prep)
    omega_index <- length(beta) + 1 + 1:length(steps)
    
    I_sub_inverse <- try_inverse(I_mat[omega_index, omega_index])
    
    Q <- if (is.null(I_sub_inverse)) NA else sum(I_sub_inverse * tcrossprod(S_vec[omega_index])) / k
    
  } else if (type == "robust") {
    
    S_mat <- null_score_matrix(beta, tau_sq, steps, y, s, X, prep = prep)
    S_vec <- colSums(S_mat)
    
    omega_index <- length(beta) + 1 + 1:length(steps)
    S_omega <- S_vec[omega_index] + prior_mass
    I_model_inv <- try_inverse(I_mat[-omega_index, -omega_index])
    
    if (is.null(I_model_inv)) {
      
      Q <- NA 
      
    } else {
      
      I_model_omega <- I_mat[omega_index, -omega_index]
      Bread <- cbind(- I_model_omega %*% I_model_inv, diag(1L, nrow = q))
      Meat <- crossprod(S_mat)
      
      V_mat <- Bread %*% Meat %*% t(Bread) / k
      
      V_inv <- try_inverse(V_mat)
      
      Q <- if (is.null(V_inv)) NA else sum(V_inv * tcrossprod(S_omega)) / k
    }
    
  } else {
    
    stop("Type must be one of 'parametric', 'subscore', or 'robust'.")
    
  }
  
  if (return == "Q") {
    return(Q)
  } else {
    p_val <- pchisq(Q, df = q, lower.tail = FALSE)
    
    data.frame(
      Q_score = Q, 
      df = q, 
      p_val = p_val
    ) %>%
      return()
      
  }

}


quick_score_Q <- function(beta, tau_sq, steps, y, s, X, prep = NULL, q = length(steps)) {
  
  if (is.null(prep)) prep <- null_prep(beta, tau_sq, steps, y, s, X)
  S_mat <- null_score_matrix(beta, tau_sq, steps, y, s, X, prep = prep)
  S_vec <- colSums(S_mat)
  
  S_cov <- crossprod(S_mat)
  
  S_cov_inv <- try_inverse(S_cov)
  
  Q <- if (is.null(S_cov_inv)) NA else sum(S_cov_inv * tcrossprod(S_vec)) / k
  p_val <- pchisq(Q, df = q, lower.tail = FALSE)
  
  data.frame(
    Q_score = Q, 
    df = q, 
    p_val = p_val
  )
  
}
