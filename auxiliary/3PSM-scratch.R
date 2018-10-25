library(tidyverse)
library(metafor)
library(weightr)
devtools::load_all()
rm(list=ls())

studies = 200
mean_effect = 0.4
sd_effect = 0.2
n_sim = n_beta(n_min = 20, n_max = 120, na = 1, nb = 3)
p_thresholds = .025
p_RR = .3

dat <- r_SMD(studies, mean_effect, sd_effect, n_sim, p_thresholds = p_thresholds, p_RR = p_RR)
mean(dat$p > .05)
model <- fit_meta_ML(dat)

steps = .025
k_min = 2
use_gradient = TRUE
method = "L-BFGS-B"
control = list()

(weightr_res <- with(dat, weightfunct(effect = dat$g, v = dat$Vg)))

RE_nogr <- refit_RE(model, use_gradient = FALSE)
RE_grad <- refit_RE(model)
VHSM_nogr <- fit_VHSM(model, use_gradient = FALSE)
VHSM_grad <- fit_VHSM(model)


# RE fits
weightr_res[[1]]$value
RE_nogr$value
RE_grad$value

weightr_res[[1]]$par
RE_nogr$par
RE_grad$par
c(model$tau2, as.vector(model$b))

RE_nogr$counts
RE_grad$counts




# VHSM fits

weightr_res[[2]]$value
VHSM_nogr$value
VHSM_grad$value

weightr_res[[2]]$par
VHSM_nogr$par
VHSM_grad$par

VHSM_nogr$counts
VHSM_grad$counts


