context("VHSM log likelihood, score, and Hessian")

rm(list=ls())
N <- 200
p <- 3
X <- cbind(rep(1, N), matrix(rnorm(N * (p - 1)), N, p - 1))
beta <- rnorm(p, mean = 0.2, sd = 0.1)
tau_sq <- 0.2^2
omega <- c(0.5, 0.4)
steps <- c(0.025, 0.5)
s <- 1 / rchisq(N, df = 5)
y <- as.vector(X %*% beta) + rnorm(N, sd = sqrt(tau_sq)) + rnorm(N, sd = s)
dat <- data.frame(y = y, s = s, X)

intercept <- matrix(rep(1, N))

library(weightr)

weightr_est <- with(dat, 
  weightfunct(effect = y, v = s^2, 
              steps = c(steps, 1), 
              table = TRUE)
) 

weightr_reg <- with(dat, 
  weightfunct(effect = y, v = s^2, 
              steps = c(steps, 1), 
              mods = ~ X2 + X3, 
              table = TRUE)
) 

test_that("likelihood ratio test agrees with weightr.", {
  
  meta_unadj <- VHSM_loglik(beta = weightr_est[[1]]$par[-1], 
                            tau_sq = weightr_est[[1]]$par[1],
                            omega = c(1,1),
                            steps = steps,
                            y = y, s = s)

  meta_adj <- VHSM_loglik(beta = weightr_est[[2]]$par[2], 
                          tau_sq = weightr_est[[2]]$par[1],
                          omega = weightr_est[[2]]$par[3:4],
                          steps = steps,
                          y = y, s = s)
  
  
  LRT_meta <- 2 * (meta_adj - meta_unadj)
  
  expect_equal(2 * (weightr_est[[1]]$value - weightr_est[[2]]$value), LRT_meta)
  
  reg_unadj <- VHSM_loglik(beta = weightr_reg[[1]]$par[-1], 
                            tau_sq = weightr_reg[[1]]$par[1],
                            omega = c(1,1),
                            steps = steps,
                            y = y, s = s, X = X)
  
  reg_adj <- VHSM_loglik(beta = weightr_reg[[2]]$par[2:4], 
                         tau_sq = weightr_reg[[2]]$par[1],
                         omega = weightr_reg[[2]]$par[5:6],
                         steps = steps,
                         y = y, s = s, X = X)
  
  
  LRT_reg <- 2 * (reg_adj - reg_unadj)
  
  expect_equal(2 * (weightr_reg[[1]]$value - weightr_reg[[2]]$value), LRT_reg)
  
})

test_that("Null score is consistent with score.", {
  
  S_VHSM_meta <- VHSM_score(beta = weightr_est[[1]]$par[-1], 
                            tau_sq = weightr_est[[1]]$par[1],
                            omega = c(1,1),
                            steps = steps,
                            y = y, s = s, X = intercept)
  
  S_null_meta <- null_score(beta = weightr_est[[1]]$par[-1], 
                            tau_sq = weightr_est[[1]]$par[1],
                            steps = steps,
                            y = y, s = s, X = intercept)
  
  expect_equal(S_VHSM_meta, S_null_meta)
  
  
  S_VHSM_reg <- VHSM_score(beta = weightr_reg[[1]]$par[-1], 
                           tau_sq = weightr_reg[[1]]$par[1],
                           omega = c(1,1),
                           steps = steps,
                           y = y, s = s, X = X)
  
  S_null_reg <- null_score(beta = weightr_reg[[1]]$par[-1], 
                       tau_sq = weightr_reg[[1]]$par[1],
                       steps = steps,
                       y = y, s = s, X = X)
  
  expect_equal(S_VHSM_reg, S_null_reg)
  
})

test_that("Null expected information is consistent with expected information", {
  
  I_VHSM_meta <- VHSM_Info(beta = weightr_est[[1]]$par[-1], 
                           tau_sq = weightr_est[[1]]$par[1],
                           omega = c(1,1),
                           steps = steps,
                           y = y, s = s, X = intercept)
  
  IO_null_meta <- null_Info(beta = weightr_est[[1]]$par[-1], 
                            tau_sq = weightr_est[[1]]$par[1],
                            steps = steps,
                            y = y, s = s, X = intercept, info = "observed")
  
  IE_null_meta <- null_Info(beta = weightr_est[[1]]$par[-1], 
                            tau_sq = weightr_est[[1]]$par[1],
                            steps = steps,
                            y = y, s = s, X = intercept, info = "expected")
  
  expect_equal(I_VHSM_meta, IO_null_meta)
  expect_type(all.equal(I_VHSM_meta, IE_null_meta), "character")
  
  
  I_VHSM <- VHSM_Info(beta = weightr_reg[[1]]$par[-1], 
                         tau_sq = weightr_reg[[1]]$par[1],
                         omega = c(1,1),
                         steps = steps,
                         y = y, s = s, X = X)
  
  IO_null <- null_Info(beta = weightr_reg[[1]]$par[-1], 
                         tau_sq = weightr_reg[[1]]$par[1],
                         steps = steps,
                         y = y, s = s, X = X, info = "observed")
  
  IE_null <- null_Info(beta = weightr_reg[[1]]$par[-1], 
                       tau_sq = weightr_reg[[1]]$par[1],
                       steps = steps,
                       y = y, s = s, X = X, info = "expected")
  
  expect_equal(I_VHSM, IO_null)
  expect_type(all.equal(I_VHSM, IE_null), "character")
  
})

test_that("score function is consistent with weightr.", {
  
  S_unadj <- VHSM_score(beta = weightr_est[[1]]$par[-1], 
                        tau_sq = weightr_est[[1]]$par[1],
                        omega = c(1,1),
                        steps = steps,
                        y = y, s = s, X = intercept)
  
  expect_true(all(abs(S_unadj[1:2]) < 10^-2))
  
  S_adj <- VHSM_score(beta = weightr_est[[2]]$par[2], 
                      tau_sq = weightr_est[[2]]$par[1],
                      omega = weightr_est[[2]]$par[3:4],
                      steps = steps,
                      y = y, s = s, X = intercept)
  
  expect_true(all(abs(S_adj) < 10^-2))
  
})

test_that("Hessian is consistent with weightr.", {
  
  H_unadj <- VHSM_Info(beta = weightr_est[[1]]$par[-1], 
                       tau_sq = weightr_est[[1]]$par[1],
                       omega = c(1,1),
                       steps = steps,
                       y = y, s = s, X = intercept)
  attr(H_unadj, "dimnames") <- NULL
  
  H_adj <- VHSM_Info(beta = weightr_est[[2]]$par[2], 
                     tau_sq = weightr_est[[2]]$par[1],
                     omega = weightr_est[[2]]$par[3:4],
                     steps = steps,
                     y = y, s = s, X = intercept)
  attr(H_adj, "dimnames") <- NULL
  
  expect_equal(weightr_est[[1]]$hessian[2:1,2:1], H_unadj[1:2, 1:2], tol = 10^-3)
  expect_equal(weightr_est[[2]]$hessian[c(2,1,3,4), c(2,1,3,4)], H_adj, tol = 10^-3)
    
})
