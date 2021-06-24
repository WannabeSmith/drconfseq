test_that("confseq_ate coincides with confseq_ate_unadjusted when train_idx set appropriately", {
  n <- 1000
  alpha <- 0.05
  t_opt <- 10
  y <- rnorm(n, mean = 1, sd = 2)
  treatment <- rbinom(n = n, size = 1, prob = 1/2)
  times <- (1:n)[rbinom(n, 1, 0.5) == 1]

  train_idx = rep(FALSE, n)
  adjusted <-
    confseq_ate(y, X = as.data.frame(matrix(rnorm(n*1), ncol=1)),
                treatment = treatment,
                regression_fn_1 = function(y, X, newX){0},
                propensity_score_fn = function(y, X, newX){1/2},
                t_opt = t_opt,
                train_idx = train_idx,
                times = times,
                n_cores = 1,
                alpha = 0.05,
                cross_fit = FALSE)

  unadjusted <-
    confseq_ate_unadjusted(y,
                           treatment = treatment,
                           propensity_score = 1/2,
                           t_opt = t_opt,
                           alpha = 0.05,
                           times = times)

  expect_equal(adjusted, unadjusted)

  # Should also work if times <- 1:n
  times <- 1:n

  train_idx = rep(FALSE, n)
  adjusted <-
    confseq_ate(y, X = as.data.frame(matrix(rnorm(n*1), ncol=1)),
                treatment = treatment,
                regression_fn_1 = function(y, X, newX){0},
                propensity_score_fn = function(y, X, newX){1/2},
                t_opt = t_opt,
                train_idx = train_idx,
                times = times,
                n_cores = 1,
                alpha = 0.05,
                cross_fit = FALSE)

  unadjusted <-
    confseq_ate_unadjusted(y,
                           treatment = treatment,
                           propensity_score = 1/2,
                           t_opt = t_opt,
                           alpha = 0.05,
                           times = times)

  expect_equal(adjusted, unadjusted)
})
