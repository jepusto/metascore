context("VHSM score tests")

rm(list=ls())
N <- 200
p <- 3
X <- matrix(rnorm(N * (p - 1)), N, p - 1)
beta <- rnorm(p, mean = 0.2, sd = 0.1)
tau_sq <- 0.2^2
omega <- c(0.5, 0.4)
steps <- c(0.025, 0.5)
s <- 1 / rchisq(N, df = 5)
y <- as.vector(cbind(rep(1, N), X) %*% beta) + rnorm(N, sd = sqrt(tau_sq)) + rnorm(N, sd = s)
dat <- data.frame(y = y, s = s, X)

library(metafor)
model <- rma(y ~ X1 + X2, sei = s, data = dat)

beta <- as.vector(model$beta)
tau_sq <- model$tau2

VHSM_score_test(model, steps = steps, type = "parametric", info = "expected")
VHSM_score_test(model, steps = steps, type = "parametric", info = "observed")
VHSM_score_test(model, steps = steps, type = "robust", info = "expected")
VHSM_score_test(model, steps = steps, type = "robust", info = "observed")
