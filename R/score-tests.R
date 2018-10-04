try_inverse <- function(x) tryCatch(chol2inv(chol(x)), error = function(e) NULL)

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
    
    I_mat_inv <- try_inverse(I_mat)
    
    Q <- if (is.null(I_mat_inv)) NA else sum(I_mat_inv * tcrossprod(S_vec))
    
  } else if (type == "subscore") {
    
    S_vec <- null_score(beta, tau_sq, steps, y, s, X, prep = prep)
    omega_index <- length(beta) + 1 + 1:length(steps)
    
    I_sub_inverse <- try_inverse(I_mat[omega_index, omega_index])
    
    Q <- if (is.null(I_sub_inverse)) NA else sum(I_sub_inverse * tcrossprod(S_vec[omega_index]))
    
  } else if (type == "robust") {
    
    S_mat <- null_score_matrix(beta, tau_sq, steps, y, s, X, prep = prep)
    S_vec <- colSums(S_mat)
    
    omega_index <- length(beta) + 1 + 1:length(steps)
    S_omega <- S_vec[omega_index]
    I_model_inv <- try_inverse(I_mat[-omega_index, -omega_index])
    
    if (is.null(I_model_inv)) {
      Q <- NA 
    } else {
      
      I_model_omega <- I_mat[omega_index, -omega_index]
      Bread <- cbind(- I_model_omega %*% I_model_inv, diag(1L, nrow = q))
      Meat <- crossprod(S_mat)
      
      V_mat <- Bread %*% Meat %*% t(Bread)
      
      V_inv <- try_inverse(V_mat)
      
      Q <- if (is.null(V_inv)) NA else sum(V_inv * tcrossprod(S_omega))
    }
    
  } else {
    
    stop("Type must be one of 'parametric', 'subscore', or 'robust'.")
    
  }
  
  p_val <- pchisq(Q, df = q, lower.tail = FALSE)

  data.frame(Q_score = Q, df = q, p_val = p_val)
}


