test_that("std_conjmix_margin fails reasonably", {
  expect_error(std_conjmix_margin(t = 10:1000, rho2 = 0, alpha = 0.05))
  expect_error(std_conjmix_margin(t = c(-1, 10:1000), rho2 = 0.1, alpha = 0.05))
  expect_error(std_conjmix_margin(t = 10:1000, rho2 = 0.1, alpha = 0))
  expect_error(std_conjmix_margin(t = 10:1000, rho2 = 0.1, alpha = 1))
  expect_silent(std_conjmix_margin(t = 10:1000, rho2 = 0.1, alpha = 0.05))
})

test_that("best_rho2_exact fails reasonably", {
  expect_silent(best_rho2_exact(t_opt = 100,
                               alpha_opt = sqrt(lamW::lambertW0(1)) -
                                 .Machine$double.eps))
  expect_error(best_rho2_exact(t_opt = 100,
                               alpha_opt = sqrt(lamW::lambertW0(1))))
})

test_that("lambertWm1_approx is close to lambertWm1", {
  expect_lte(abs(lambertWm1_approx(-10^(-20)) - lambertWm1(-10^(-20))), 0.1)
  expect_lte(abs(lambertWm1_approx(-10^(-50)) - lambertWm1(-10^(-50))), 0.05)
  expect_lte(abs(lambertWm1_approx(-10^(-100)) - lambertWm1(-10^(-100))), 0.03)
})

test_that("best_rho2_approx is close to best_rho2_exact", {
  expect_lt(abs(best_rho2_approx(t_opt = 1, alpha_opt = 0.05) -
                  best_rho2_exact(t_opt = 1, alpha_opt = 0.05)), 0.5)
  expect_lt(abs(best_rho2_approx(t_opt = 1, alpha_opt = 0.001) -
                  best_rho2_exact(t_opt = 1, alpha_opt = 0.001)), 0.2)
  expect_lt(abs(best_rho2_approx(t_opt = 1, alpha_opt = 0.00001) -
                  best_rho2_exact(t_opt = 1, alpha_opt = 0.00001)), 0.13)
  expect_lt(abs(best_rho2_approx(t_opt = 1, alpha_opt = 10^(-10)) -
                  best_rho2_exact(t_opt = 1, alpha_opt = 10^(-10))), 0.09)
})

test_that("asymptotic_confseq has valid coverage", {
  n <- 10000
  repeats <- 500
  alpha <- 0.1
  # Generate standard Gaussian data

  miscoverages <- purrr::map(1:repeats, function(i){
    x <- rnorm(n)
    acs <- asymptotic_confseq(x, t_opt = 10, alpha = alpha, var = 1)
    l = acs$l
    u = acs$u

    return(any(l > 0 | u < 0))
  }) %>%
    unlist() %>%
    sum()

  # Need to have a small number of miscoverages.
  expect_lt(miscoverages, repeats)
  # Miscoverage shouldn't be near 0,
  # otherwise our CSs are overly conservative
  expect_gt(miscoverages, 0)

  coverage_test <- binom.test(miscoverages,
                              n = repeats,
                              p = 0.1,
                              alternative = "greater",
                              # This random test will fail at
                              # most 1% of the time
                              conf.level = 0.99)
  # Testing to see whether the lower binomial
  # confidence interval is less than alpha
  expect_lt(coverage_test$conf.int[1], alpha)
})

test_that("plot_cs_shape does not throw an error", {
  expect_silent(plot_cs_shape(t_opts = c(1, 10, 50, 500, 1000),
                              t = 1:10000, alpha = 0.05, log_scale = FALSE))
  # Should work even if only one t_opts is passed
  expect_silent(plot_cs_shape(t_opts = c(1),
                              t = 1:10000, alpha = 0.05, log_scale = FALSE))
  # Should work even if only one t is passed
  expect_silent(plot_cs_shape(t_opts = c(1, 10, 50, 500, 1000),
                              t = 1, alpha = 0.05, log_scale = FALSE))
  # Make sure log_scale = TRUE also works
  expect_silent(plot_cs_shape(t_opts = c(1, 10, 50, 500, 1000),
                              t = 1:10000, alpha = 0.05, log_scale = TRUE))
})
