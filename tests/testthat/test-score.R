context("VHSM score tests")

dat <- r_SMD(studies = 100, mean_effect = 0.1, sd_effect = 0.1,
             n_sim = n_beta(10, 50, na = 3, nb = 3),
             p_thresholds = .025, p_RR = 0.5)

dat$X1 <- rnorm(100)
dat$X2 <- sample(LETTERS[1:4], size = 100, replace = TRUE)

test_types <- 
  list(
    type = c("parametric","robust"), 
    info = c("expected","observed"),
    steps = list(one = .025, two = c(.025, .500))
  ) %>%
  purrr:cross()


test_that("Score tests work for intercept only models.", {
  
  meta_fit <- rma(g, sei = sda, data = dat, method = "ML")
  
  score_tests <- map_dfr(test_types, .f = ~ VHSM_score_test(meta_fit, steps = .$steps, type = .$type, info = .$info))
  
  expect_true(is.data.frame(score_tests))
})

test_that("Score tests work for meta-regression models.", {
  
  meta_reg <- rma(g ~ X1 + X2, sei = sda, data = dat, method = "ML")
  
  score_tests <- map_dfr(test_types, .f = ~ VHSM_score_test(meta_reg, steps = .$steps, type = .$type, info = .$info))
  
  expect_true(is.data.frame(score_tests))
})

