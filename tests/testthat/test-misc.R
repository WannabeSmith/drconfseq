test_that("cumul_mean works", {
  x <- rnorm(100)
  expect_equal(cumul_mean(x)[length(x)], mean(x))
  expect_equal(cumul_mean(x)[1], x[1])
  expect_equal(length(cumul_mean(x)), length(x))
})

test_that("cumul_var works", {
  x <- rnorm(100)
  expect_equal(cumul_var(x)[length(x)], var(x))
  expect_equal(cumul_var(x)[1], NaN)
  expect_equal(length(cumul_var(x)), length(x))
})

test_that("naive_confidence_intervals works", {
  n <- 100
  alpha <- 0.1
  x <- rnorm(n)
  ci_all_times <- naive_confidence_intervals(x, alpha = alpha,
                                             var = 1, return_all_times = TRUE)
  ci_final_time <- naive_confidence_intervals(x, alpha = alpha,
                                             var = 1, return_all_times = FALSE)
  expect_equal(ci_all_times[n], ci_final_time[n])

  expect_equal(ci_final_time$l, mean(x) - qnorm(1-alpha/2) / sqrt(n))
})
