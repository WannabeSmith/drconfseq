test_that("pseudo_outcome_sequential sample splitting works", {
  n <- 1000
  alpha <- 0.05
  t_opt <- 10
  y <- rnorm(n, mean = 1, sd = 2)
  treatment <- rbinom(n = n, size = 1, prob = 1/2)
  times <- 1:n

  train_idx = rbinom(n, 1, 1/2)
  manual_sample_split_cs <-
    c(pseudo_outcome_sequential(y,
                              X = as.data.frame(matrix(rnorm(n*1), ncol=1)),
                              treatment = treatment,
                              regression_fn_1 = function(y, X, newX){0},
                              propensity_score_fn = function(y, X, newX){1/2},
                              train_idx = train_idx,
                              cross_fit = FALSE,
                              times = 1:n)[[n]],
      pseudo_outcome_sequential(y,
                              X = as.data.frame(matrix(rnorm(n*1), ncol=1)),
                              treatment = treatment,
                              regression_fn_1 = function(y, X, newX){0},
                              propensity_score_fn = function(y, X, newX){1/2},
                              train_idx = 1-train_idx,
                              cross_fit = FALSE,
                              times = 1:n)[[n]])

  auto_sample_split_cs <-
    pseudo_outcome_sequential(y,
                              X = as.data.frame(matrix(rnorm(n*1), ncol=1)),
                              treatment = treatment,
                              regression_fn_1 = function(y, X, newX){0},
                              propensity_score_fn = function(y, X, newX){1/2},
                              train_idx = train_idx,
                              cross_fit = TRUE,
                              times = 1:n)[[n]]

  expect_equal(manual_sample_split_cs, auto_sample_split_cs)
})
