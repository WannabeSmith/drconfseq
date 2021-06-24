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
  ci_all_times <- naive_confidence_intervals(x,
                                             alpha = alpha,
                                             var = 1,
                                             return_all_times = TRUE)
  ci_final_time <- naive_confidence_intervals(x,
                                              alpha = alpha,
                                              var = 1,
                                              return_all_times = FALSE)
  expect_equal(ci_all_times[n], ci_final_time[n])

  expect_equal(ci_final_time$l, mean(x) - qnorm(1 - alpha / 2) / sqrt(n))
})

test_that("get_cumul_miscoverage_rate works", {
  n <- 1000
  num_repeats <- 500
  var <- 10
  alpha <- 0.1
  miscoverage_rate <-
    get_cumul_miscoverage_rate(
      data_generator_fn = function() {
        rnorm(n, sd = sqrt(var))
      },
      conf_set_fn = function(x) {
        asymptotic_confseq(x,
                           t_opt = 10,
                           alpha = alpha,
                           var = var)
      },
      times = 1:n,
      num_repeats = num_repeats,
      mu = 0,
      n_cores = 1
    )

  # Check to make sure that miscoverage is nondecreasing
  expect_true(all(miscoverage_rate - c(0, miscoverage_rate[0:(n - 1)]) >= 0))

  # Check to make sure miscoverage is nonzero at the final time.
  # This doesn't *always* need to be true, but in almost all cases, it will.
  expect_true(miscoverage_rate[n] > 0)

  # Calculate confidence interval for miscoverage rate.
  # This test will fail at most 1% of the time.
  coverage_test <-
    binom.test(
      miscoverage_rate[num_repeats] * num_repeats,
      n = num_repeats,
      p = alpha,
      alternative = "greater",
      conf.level = 0.99
    )
  expect_lt(coverage_test$conf.int[1], alpha)
})
