
VHSM_score_test <- function(model, steps, type = "parametric", info = "expected") {
  
  if (!("rma.uni" %in% class(model))) {
    return(data.frame(Q_score = NA, df = NA, p_val = NA))
  }

  beta <- as.vector(model$beta)
  tau_sq <- model$tau2
  y <- as.vector(model$yi)
  s <- sqrt(model$vi)
  X <- model$X
  
  prep <- null_prep(beta, tau_sq, steps, y, s, X)
  
  
  I_mat <- null_Info(beta, tau_sq, steps, y, s, X, prep = prep, info = info)
  
  q <- length(steps)
  
  if (type == "parametric") { 
    
    S_vec <- null_score(beta, tau_sq, steps, y, s, X, prep = prep)
    
    Q <- sum(chol2inv(chol(I_mat)) * tcrossprod(S_vec))
    
  } else {
    
    S_mat <- null_score_matrix(beta, tau_sq, steps, y, s, X, prep = prep)
    S_vec <- colSums(S_mat)
    
    omega_index <- length(beta) + 1 + 1:length(steps)
    S_omega <- S_vec[omega_index]
    I_model <- I_mat[-omega_index, -omega_index]
    I_model_omega <- I_mat[omega_index, -omega_index]
    Bread <- cbind(- I_model_omega %*% chol2inv(chol(I_model)), diag(1L, nrow = q))
    Meat <- crossprod(S_mat)
    
    V_mat <- Bread %*% Meat %*% t(Bread)
    
    Q <- sum(chol2inv(chol(V_mat)) * tcrossprod(S_omega))
    
  }
  
  p_val <- pchisq(Q, df = q, lower.tail = FALSE)

  data.frame(Q_score = Q, df = q, p_val = p_val)
}


