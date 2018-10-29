library(rlang)
library(purrr)

generate_data <- function(n) {
  data.frame(y = rnorm(n), x = rpois(n, lambda = 3))
}

generate_data(5)

lotsa_data <- rerun(10, generate_data(10))

fit_model <- function(dat) {
  dat_expr <- enexpr(dat)
  lm_call <- expr(lm(y ~ x, data = !!dat_expr))
  eval_bare(lm_call, caller_env())
}

some_data <- generate_data(10)
fit_model(some_data)

lotsa_models <- map(lotsa_data, fit_model)

update_model <- function(mod) {
  update(mod, formula = . ~ . + I(x^2))
}

some_mod <- fit_model(some_data)
update_model(some_mod)

map(lotsa_models, update_model)
