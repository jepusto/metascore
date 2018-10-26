library(tidyverse)
library(metafor)
library(weightr)
devtools::load_all()
rm(list=ls())

studies = 50
mean_effect = 0.9
sd_effect = 0.2
n_sim = n_beta(n_min = 20, n_max = 120, na = 1, nb = 3)
p_thresholds = .025
p_RR = 1

dat <- r_SMD(studies, mean_effect, sd_effect, n_sim, p_thresholds = p_thresholds, p_RR = p_RR)
mean(dat$p > .05)
model <- fit_meta_ML(dat)

(weightr_res <- with(dat, weightfunct(effect = dat$g, v = dat$Vg, steps = c(steps, 1))))
LRT_VHSM(model, steps = steps)

y <- model$yi
s <- sqrt(model$vi)
X <- model$X
steps = c(0.01, .025, .5)
k_min = 2
tol = 10^-3
use_gradient = TRUE
method = "L-BFGS-B"
control = list()


RE_nogr <- fit_VHSM(y = y, s = s, X = X, steps = NULL, use_gradient = FALSE)
RE_grad <- fit_VHSM(y = y, s = s, X = X, steps = NULL)
VHSM_nogr <- fit_VHSM(y = y, s = s, X = X, steps = steps, use_gradient = FALSE)
VHSM_grad <- fit_VHSM(y = y, s = s, X = X, steps = steps)

# RE fits

weightr_res[[1]]$value
VHSM_negloglik_theta(theta = weightr_res[[1]]$par, steps = NULL, 
                     y =y, s = s)
RE_nogr$value
RE_grad$value

c(model$tau2, as.vector(model$b))
weightr_res[[1]]$par
RE_nogr$par
RE_grad$par

# VHSM fits

weightr_res[[2]]$value
VHSM_negloglik_theta(theta = weightr_res[[2]]$par, steps = steps, 
                     y = model$yi, s = sqrt(model$vi))
VHSM_nogr$value
VHSM_grad$value

weightr_res[[2]]$par
VHSM_nogr$par
VHSM_grad$par
VHSM_neg_score_theta(theta = weightr_res[[2]]$par, steps = steps, 
                     y = y, s = s, X = X)
VHSM_neg_score_theta(theta = VHSM_grad$par, steps = steps, 
                     y = y, s = s, X = X)


prep <- null_prep(beta = as.vector(model$b), 
                  tau_sq = model$tau2, 
                  steps = steps,
                  y = y, s = s, X = X)

with(prep, as.vector(colSums((er / eta^2) * X)))
with(prep, sum(er^2 / eta^4) / 2 - sum(1 / eta^2) / 2)
with(prep, n_s[-1] - colSums(B_mat[,-1,drop=FALSE]))
