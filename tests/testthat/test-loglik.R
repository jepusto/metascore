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

library(weightr)

weightr_est <- with(dat, 
  weightfunct(effect = y, v = s^2, 
              steps = c(steps, 1), 
              mods = ~ X2 + X3, 
              table = TRUE)
) 

test_that("likelihood ratio test agrees with weightr.", {
  
  LRT_weightr <- 2 * (weightr_est[[1]]$value - weightr_est[[2]]$value)
  
  ll_unadj <- VHSM_loglik(beta = weightr_est[[1]]$`par`[-1], 
                          tau_sq = weightr_est[[1]]$`par`[1],
                          omega = c(1,1),
                          steps = steps,
                          y = y, s = s, X = X)

  ll_adj <- VHSM_loglik(beta = weightr_est[[2]]$`par`[2:4], 
                        tau_sq = weightr_est[[2]]$`par`[1],
                        omega = weightr_est[[2]]$`par`[5:6],
                        steps = steps,
                        y = y, s = s, X = X)
  
  
  LRT_VHSM <- 2 * (ll_adj - ll_unadj)
  
  expect_equal(LRT_weightr, LRT_VHSM)
  
})

test_that("Null score is consistent with score.", {
  
  S_VHSM <- VHSM_score(beta = weightr_est[[1]]$par[-1], 
                       tau_sq = weightr_est[[1]]$par[1],
                       omega = c(1,1),
                       steps = steps,
                       y = y, s = s, X = X)
  
  S_null <- null_score(beta = weightr_est[[1]]$par[-1], 
                       tau_sq = weightr_est[[1]]$par[1],
                       steps = steps,
                       y = y, s = s, X = X)
  
  expect_equal(S_VHSM, S_null)
  
})

test_that("Null Hessian is consistent with score.", {
  
  H_VHSM <- VHSM_Hessian(beta = weightr_est[[1]]$par[-1], 
                         tau_sq = weightr_est[[1]]$par[1],
                         omega = c(1,1),
                         steps = steps,
                         y = y, s = s, X = X)
  
  H_null <- null_Hessian(beta = weightr_est[[1]]$par[-1], 
                         tau_sq = weightr_est[[1]]$par[1],
                         steps = steps,
                         y = y, s = s, X = X)
  
  expect_equal(H_VHSM, H_null)
  
})

test_that("score function is consistent with weightr.", {
  
  skip("Score function")
  
  S_unadj <- VHSM_score(beta = weightr_est[[1]]$par[-1], 
                        tau_sq = weightr_est[[1]]$par[1],
                        omega = c(1,1),
                        steps = steps,
                        y = y, s = s, X = X)
  
  expect_true(all(abs(S_unadj[1:4]) < 10^-2))
  
  S_adj <- VHSM_score(beta = weightr_est[[2]]$par[2:4], 
                      tau_sq = weightr_est[[2]]$par[1],
                      omega = weightr_est[[2]]$par[5:6],
                      steps = steps,
                      y = y, s = s, X = X)
  
  expect_true(all(abs(S_adj) < 10^-2))
  
})

test_that("Hessian is consistent with weightr.", {
  
  skip("Hessian")
  
  H_unadj <- VHSM_Hessian(beta = weightr_est[[1]]$`par`[-1], 
                          tau_sq = weightr_est[[1]]$`par`[1],
                          omega = c(1,1),
                          steps = steps,
                          y = y, s = s, X = X)
  
  H_adj <- VHSM_Hessian(beta = weightr_est[[2]]$`par`[2:4], 
                        tau_sq = weightr_est[[2]]$`par`[1],
                        omega = weightr_est[[2]]$`par`[5:6],
                        steps = steps,
                        y = y, s = s, X = X)
    
})