test_that("std_conjmix_margin fails reasonably", {
  expect_error(std_conjmix_margin(t = 10:1000, rho2 = 0, alpha = 0.05))
  expect_error(std_conjmix_margin(t = c(-1, 10:1000), rho2 = 0.1, alpha = 0.05))
  expect_error(std_conjmix_margin(t = 10:1000, rho2 = 0.1, alpha = 0))
  expect_error(std_conjmix_margin(t = 10:1000, rho2 = 0.1, alpha = 1))
  expect_silent(std_conjmix_margin(t = 10:1000, rho2 = 0.1, alpha = 0.05))
})

