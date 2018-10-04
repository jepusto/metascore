context("Viechtbauer (2005).")
library(metafor)
library(purrr)
library(dplyr)

dat <- r_SMD(studies = 100, mean_effect = 0.1, sd_effect = 0.1,
             n_sim = n_beta(10, 50, na = 3, nb = 3),
             p_thresholds = .025, p_RR = 0.5)

meta_fit <- rma(g, sei = sda, data = dat, method = "ML")

mu <- as.vector(meta_fit$b)
tau_sq <- meta_fit$tau2
intercept <- matrix(rep(1, nrow(dat)))
steps <- c(.05, .5)

test_that("Null observed information agrees with Viechtbauer (2005).", {
  
  I_VHSM_meta <- null_Info(beta = mu, tau_sq = tau_sq, steps = steps,
                           y = dat$g, s = dat$sda, X = intercept)
  
  w <- 1 / (meta_fit$tau2 + meta_fit$vi)
  e <- meta_fit$yi - mu
    
  I_mu_mu <- sum(w)
  I_mu_tau <- sum(w^2 * e)
  I_tau_tau <- sum(w^3 *e^2) - sum(w^2) / 2
  I_ML <- matrix(c(I_mu_mu, I_mu_tau, I_mu_tau, I_tau_tau), 2, 2)

  V_ML <- 2 / sum(w^2)
  
  expect_equal(I_VHSM_meta[1:2,1:2], I_ML)
  
})
