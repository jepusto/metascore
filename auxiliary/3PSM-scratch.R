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
model <- fit_meta(dat)

step = .025
k_min = 2
method = "L-BFGS-B"
control = list()

(weightr_res <- with(dat, weightfunct(effect = dat$g, v = dat$Vg)))
PSM_fit <- fit_3PSM(model)

weightr_res[[1]]$par
c(model$tau2, as.vector(model$b))
PSM_fit$RE$par[2:1]

weightr_res[[2]]$par[c(2,1,3)]
PSM_fit$`3PSM`$par

weightr_res
PSM_fit$LRT
PSM_fit$p_val
