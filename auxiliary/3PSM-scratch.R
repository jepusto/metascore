library(tidyverse)
library(metafor)
library(weightr)
devtools::load_all()
rm(list=ls())

studies = 50
mean_effect = 1.0
sd_effect = 0.1
n_sim = n_beta(n_min = 20, n_max = 120, na = 1, nb = 3)
p_thresholds = .025
p_RR = 0.1
steps = .025


dat <- r_SMD(studies, mean_effect, sd_effect, n_sim, p_thresholds = p_thresholds, p_RR = p_RR)
sum(dat$p > .05)
model <- fit_meta_ML(dat)
tol = 10^-3
use_gradient = TRUE
method = "L-BFGS-B"
control = list()
k_min = 2L

(weightr_res <- with(dat, weightfunct(effect = dat$g, v = dat$Vg, steps = c(steps, 1), table = TRUE)))
(VHSM_res <- LRT_VHSM(model, steps = steps, k_min = k_min))

LRT_VHSM(model, steps = steps, k_min = 0)
LRT_VHSM(model, steps = steps, k_min = 1)$VHSM[[1]]$par
(one_sided <- LRT_VHSM(model, steps = steps, k_min = k_min, two_sided = FALSE))
one_sided$VHSM[[1]]$par

y <- model$yi
s <- sqrt(model$vi)
X <- model$X

p_vals <- pnorm(y / s, lower.tail = FALSE)
p_cat(p_vals, steps)
new_steps <- find_new_steps(p_vals = p_vals, steps = steps, k_min = k_min)
p_cat(p_vals, new_steps)

RE_nogr <- fit_VHSM(y = y, s = s, X = X, steps = NULL, use_gradient = FALSE)
RE_grad <- fit_VHSM(y = y, s = s, X = X, steps = NULL)
VHSM_nogr <- fit_VHSM(y = y, s = s, X = X, steps = steps, use_gradient = FALSE)
VHSM_grad <- fit_VHSM(y = y, s = s, X = X, steps = steps)

VHSM_res$VHSM[[1]]$par
VHSM_grad$par

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
