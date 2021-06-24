test_that("pseudo_outcome_sequential sample splitting works", {
  n <- 1000
  alpha <- 0.05
  t_opt <- 10
  y <- rnorm(n, mean = 1, sd = 2)
  treatment <- rbinom(n = n, size = 1, prob = 1/2)
  times <- 1:n

  train_idx = sample(c(rep(1, 0.4*n), rep(0, 0.6*n)), replace = FALSE)

  # Create pseudo-outcomes using the above train_idx split
  pseudo_outcomes_trainsplit <-
    pseudo_outcome_sequential(y,
                              X = as.data.frame(matrix(rnorm(n*1), ncol=1)),
                              treatment = treatment,
                              regression_fn_1 = function(y, X, newX){0},
                              propensity_score_fn = function(y, X, newX){1/2},
                              train_idx = train_idx,
                              cross_fit = FALSE,
                              times = 1:n)
  # Create the pseudo-outcomes using the opposite of the above train_idx split
  pseudo_outcomes_evalsplit <-
    pseudo_outcome_sequential(y,
                              X = as.data.frame(matrix(rnorm(n*1), ncol=1)),
                              treatment = treatment,
                              regression_fn_1 = function(y, X, newX){0},
                              propensity_score_fn = function(y, X, newX){1/2},
                              train_idx = 1-train_idx,
                              cross_fit = FALSE,
                              times = 1:n)

  # pseudo_outcomes_trainsplit at the final time should have exactly the number
  # of *evaluation* split subjects.
  expect_equal(length(pseudo_outcomes_trainsplit[[n]]), sum(1-train_idx))
  # Opposite should be true for pseudo_outcomes_eval.
  expect_equal(length(pseudo_outcomes_evalsplit[[n]]), sum(train_idx))

  # Create pseudo-outcomes using the above
  # train_idx split butwith cross-fitting
  auto_sample_splits <-
    pseudo_outcome_sequential(y,
                              X = as.data.frame(matrix(rnorm(n*1), ncol=1)),
                              treatment = treatment,
                              regression_fn_1 = function(y, X, newX){0},
                              propensity_score_fn = function(y, X, newX){1/2},
                              train_idx = train_idx,
                              cross_fit = TRUE,
                              times = 1:n)

  auto_sample_split_final <- auto_sample_splits[[n]]

  manual_sample_split_final <-
    c(pseudo_outcomes_trainsplit[[n]], pseudo_outcomes_evalsplit[[n]])

  # Splitting manually and concatenating training and evaluation sets should be
  # equivalent to performing the sample split automatically.
  expect_equal(manual_sample_split_final, auto_sample_split_final)

  lengths <- unname(unlist(purrr::map(auto_sample_splits, length)))
  # At time t, the cross-fit set
  # should have exactly t pseudooutcomes.
  expect_equal(lengths, 1:n)
})
