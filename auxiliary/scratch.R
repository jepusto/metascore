library(tidyverse)
devtools::load_all()
rm(list=ls())

test_types <- 
  list(type = c("parametric","subscore","robust"), info = c("expected","observed")) %>%
  cross()

n_beta(20, 120, 1, 3)(1000) %>% hist()

runSim(reps = 4000, studies = 80, mean_effect = 0.1, sd_effect = 0.0,
       n_sim = n_beta(20, 120, 1, 3), 
       p_thresholds = .025, p_RR = 1, test_types = test_types)


reps <- 1000
studies = 100
mean_effect = 0.1
sd_effect = 0.0
n_sim = n_beta(n_min = 20, n_max = 120, na = 1, nb = 3)
p_thresholds = .025
p_RR = 1

plot(density(n_sim(1000)))
#--------------------------------
# runSim innards 
#--------------------------------

meta_dat <- rerun(reps, r_SMD(studies, mean_effect, sd_effect, n_sim, p_thresholds = p_thresholds, p_RR = p_RR))

test_results <- map_df(meta_dat, estimate_effects, test_types = test_types)

test_results %>%
  group_by(type, info) %>% 
  summarise(
    pct_NA = mean(is.na(p_val)),
    reject_025 = mean(p_val < .025, na.rm = TRUE),
    reject_050 = mean(p_val < .050, na.rm = TRUE),
    reject_100 = mean(p_val < .100, na.rm = TRUE)
  )

for (i in 1:reps) {
  estimate_effects(dat = meta_dat[[i]], test_types = test_types)
}

#--------------------------------
# estimate_effects innards 
#--------------------------------

dat <- meta_dat[[i]]
test_steps <- .025
max_iter <- 100
step_adj <- 1
tau2_min <- -min(dat$Va)


rma_ML <- rma(yi = g, vi = Va, data = dat, method = "ML", 
              control = list(maxiter = max_iter, stepadj = step_adj, tau2.min = tau2_min))
rma_ML

score_tests <- map_dfr(test_types, .f = ~ VHSM_score_test(rma_ML, steps = test_steps, 
                                                          type = .$type, info = .$info))
map_dfr(test_types, as_data_frame) %>%
  bind_cols(score_tests)

#--------------------------------
# VHSM_score_test innards
#--------------------------------

model <- rma_ML
info <- "observed"
steps <- test_steps

beta <- as.vector(model$beta)
tau_sq <- model$tau2
y <- as.vector(model$yi)
s <- sqrt(model$vi)
X <- model$X

prep <- null_prep(beta, tau_sq, steps, y, s, X)

I_mat <- null_Info(beta, tau_sq, steps, y, s, X, prep = prep, info = info)

q <- length(steps)

# parametric

S_vec <- null_score(beta, tau_sq, steps, y, s, X, prep = prep)

I_mat_inv <- try_inverse(I_mat)

Q <- if (is.null(I_mat_inv)) NA else sum(I_mat_inv * tcrossprod(S_vec))

# subscore

S_vec <- null_score(beta, tau_sq, steps, y, s, X, prep = prep)
omega_index <- length(beta) + 1 + 1:length(steps)

I_sub_inverse <- try_inverse(I_mat[omega_index, omega_index])

Q <- if (is.null(I_sub_inverse)) NA else sum(I_sub_inverse * tcrossprod(S_vec[omega_index]))

# robust

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

