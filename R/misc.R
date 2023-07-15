#' Cumulative mean
#'
#' @param x The real-valued data
#' @param regularizer_obs How many fake "regularization" observations to add.
#'                        Setting this to 0 corresponds to the usual cumulative
#'                        sample mean, while larger values lead to more
#'                        regularization.
#' @param regularizer_mean The value of the fake regularization observations.
#'                         If `regularizer_obs` is 0, then this parameter has
#'                         no effect on the output of `cumul_mean`.
#'
#' @return The sample mean of `x` at each time
#' @export
cumul_mean <- function(x,
                       regularizer_obs = 0,
                       regularizer_mean = 1 / 2) {
  t <- seq(1, length(x))
  return((cumsum(x) + regularizer_obs * regularizer_mean) / (t + regularizer_obs))
}

#' Cumulative variance
#'
#' @param x The real-valued data
#'
#' @return The sample variance of `x` at each time
#' @export
cumul_var <- function(x) {
  t <- 1:length(x)
  sigma2_t <- (cumul_mean(x^2) - cumul_mean(x)^2) * t / (t - 1)

  return(sigma2_t)
}

#' Confidence interval margin for the mean of unit variance Gaussians
#'
#' @importFrom stats qnorm
#' @param t The times at which the margin will be evaluated
#' @param alpha The (0, 1)-valued confidence level
#'
#' @return A vector containing the margin for each $t$
#' @export
naive_std_margin <- function(t, alpha) {
  return(qnorm(p = 1 - alpha) / sqrt(t))
}

#' Naive confidence interval
#'
#' @param x The observed data points (a real vector).
#' @param alpha The significance level (a (0, 1)-valued real).
#' @param var The known or estimated variance of the observations, `x`.
#'            If left `NULL`, then the empirical variance of `x` will be
#'            taken (positive real-valued vector or `NULL`).
#' @param return_all_times Should the confidence sequence be returned at each
#'                         time? (boolean)
#' @return A list containing the following vectors: \cr
#' \item{l}{The lower confidence interval (a real vector).}
#' \item{u}{The upper confidence interval (a real vector).}
#' @export
naive_confidence_intervals <- function(x,
                                       alpha = 0.05,
                                       var = NULL,
                                       return_all_times = FALSE) {
  # If the user wants results at each time, use time from 1 to n.
  # Otherwise, just use time n.
  if (return_all_times) {
    t <- seq(1, length(x))
    mu_hat_t <- cumul_mean(x)
    if (is.null(var)) {
      var <- cumul_var(x)
    }
  } else {
    t <- length(x)
    mu_hat_t <- mean(x)
    if (is.null(var)) {
      var <- var(x)
    }
  }

  std_margin <- naive_std_margin(t = t, alpha = alpha / 2)

  margin <- sqrt(var) * std_margin
  # When the margin is NA, such as on the first observation, just return
  # infinity
  margin[is.na(margin)] <- Inf

  return(list(
    "l" = mu_hat_t - margin,
    "u" = mu_hat_t + margin
  ))
}

#' Get the cumulative empirical miscoverage rate of a confidence sequence
#'
#' @param data_generator_fn A function which generates a vector of data
#'                          (() -> real)
#' @param conf_set_fn A function which takes a vector of data and produces a
#'                    sequence of confidence intervals or a confidence sequence
#'                    (real -> (real, real))
#' @param times The times for which the miscoverage should be considered.
#'              For example, we might want to start the experiment late since
#'              early start times can lead to poor asymptotic approximations
#'              (a vector of positive integers)
#' @param num_repeats How many repeats to perform. Higher values lead to
#'                    better estimates of the miscoverage rate (positive int).
#' @param mu The parameter which `conf_set_fn` should be covering (a real).
#' @param n_cores How many cores to parallel process on (positive int).
#' @return A vector of (increasing) miscoverage rates at each time
#'         (a (0, 1)-valued vector).
#' @export
get_cumul_miscoverage_rate <-
  function(data_generator_fn,
           conf_set_fn,
           times,
           num_repeats,
           mu = 0,
           n_cores = 1) {
    miscoverage_list <- mclapply(1:num_repeats, function(i) {
      x <- data_generator_fn()
      conf_sets <- conf_set_fn(x)
      l <- conf_sets$l[times]
      u <- conf_sets$u[times]
      miscoverage <- cummax(l > mu | u < mu)
      stopifnot(all(!is.na(miscoverage)))

      if (i %% ceiling(0.1 * num_repeats) == 0) {
        print(paste("Finished simulation", i, "out of", num_repeats, sep = " "))
      }
      return(miscoverage)
    }, mc.cores = n_cores)

    miscoverage_rate <- colMeans(do.call(rbind, miscoverage_list))

    return(miscoverage_rate)
  }

#' Get the empirical width of a confidence sequence
#'
#' @importFrom parallel mclapply
#' @param data_generator_fn A function which generates a vector of data
#'                          (() -> real)
#' @param conf_set_fn A function which takes a vector of data and produces a
#'                    sequence of confidence intervals or a confidence sequence
#'                    (real -> (real, real))
#' @param times The times for which the miscoverage should be considered.
#'              For example, we might want to start the experiment late since
#'              early start times can lead to poor asymptotic approximations
#'              (a vector of positive integers)
#' @param num_repeats How many repeats to perform. Higher values lead to
#'                    better estimates of the miscoverage rate (positive int).
#' @param n_cores How many cores to parallel process on (positive int).
#' @return A vector of (increasing) miscoverage rates at each time
#'         (a (0, 1)-valued vector).
#' @export
get_avg_width <-
  function(data_generator_fn,
           conf_set_fn,
           times,
           num_repeats,
           n_cores = 1) {
    width_list <- mclapply(1:num_repeats, function(i) {
      x <- data_generator_fn()
      conf_sets <- conf_set_fn(x)
      l <- conf_sets$l[times]
      u <- conf_sets$u[times]
      width <- u - l
      stopifnot(all(!is.na(width)))

      if (i %% ceiling(0.1 * num_repeats) == 0) {
        print(paste("Finished simulation", i, "out of", num_repeats, sep = " "))
      }
      return(width)
    }, mc.cores = n_cores)

    avg_width <- colMeans(do.call(rbind, width_list))

    return(avg_width)
  }


#' Robbins' sub-Gaussian confidence sequence
#'
#' @param x A vector of data (a real vector).
#' @param alpha The desired type-I error level (a real in (0, 1)).
#' @param t_opt The time to optimize the CS for (a positive integer).
#' @param sigma2 The sub-Gaussian parameter (a real greater than 0).
#' @return A list with two elements, `l` and `u`, which are vectors of lower
#'        and upper confidence bounds, respectively.
#' @export
robbins_subGaussian_cs <- function(x, alpha, t_opt, sigma2) {
  # Yes, we're using the code of AsympCSs but when plugging in
  # the true sub-Gaussian parameter, this becomes the nonasymptotic
  # sub-Gaussian CS of Robbins.
  cs <- asymptotic_confseq(x = x, t_opt = t_opt, alpha = alpha, var = sigma2, return_all_times = TRUE)
  return(list(
    "l" = cs$l,
    "u" = cs$u
  ))
}

#' Howard et al. sub-exponential confidence sequence
#'
#' @importFrom reticulate import
#' @param x A vector of data in [0, 1] (a real vector).
#' @param alpha The desired type-I error level (a real in (0, 1)).
#' @param v_opt The time to optimize the CS for (a positive integer).
#' @return A list with two elements, `l` and `u`, which are vectors of lower
#'        and upper confidence bounds, respectively.
#' @export
howard_subexponential_cs <- function(x, alpha, v_opt, running_intersection = FALSE, lower_bd = 0, upper_bd = 1) {
  confseq <- import("confseq")

  cs <- confseq$conjmix_bounded$conjmix_empbern_twosided_cs(
    x = x, v_opt = v_opt,
    alpha = alpha,
    running_intersection = running_intersection,
    lower_bd = lower_bd, upper_bd = upper_bd
  )
  ## l <- confseq$conjmix_bounded$conjmix_empbern_lower_cs(
  ##   x = x, v_opt = v_opt, alpha = alpha, running_intersection = running_intersection
  ## )

  ## u <- 1 - confseq$conjmix_bounded$conjmix_empbern_lower_cs(
  ##   x = 1 - x, v_opt = v_opt, alpha = alpha, running_intersection = running_intersection
  ## )

  return(list("l" = cs[[1]], "u" = cs[[2]]))
}


#' Betting confidence sequence
#'
#' @param x A vector of data (a real vector).
#' @param alpha The desired type-I error level (a real in (0, 1)).
#' @param breaks The number of breaks to use for the CS.
#' @param running_intersection Whether to use the running intersection
#' @param parallel Whether to use parallel processing
#' @param convex_comb Whether to use convex combination
#' @param theta The convex combination parameter
#' @param trunc_scale The truncation scale
#' @param m_trunc Whether to scale truncation based on m
#' @return A list with two elements, `l` and `u`, which are vectors of lower
#'        and upper confidence bounds, respectively.
#' @export
betting_cs <- function(x, alpha, breaks = 1000,
                       running_intersection = FALSE, parallel = FALSE,
                       convex_comb = FALSE,
                       theta = 1 / 2, trunc_scale = 1 / 2, m_trunc = TRUE) {
  confseq <- import("confseq")
  cs <- confseq$betting$betting_cs(
    x = x,
    alpha = alpha
  )
  return(list(
    "l" = cs[[1]],
    "u" = cs[[2]]
  ))
}

#' Get a sequence of logarithmically-spaced numbers
#'
#' @param from The starting value of the sequence
#' @param to The final value of the sequence
#' @param n How many values of the sequence to output
#' @return A sequence of logarithmically-spaced numbers from `from` to `to`.
#' @export
logseq <- function(from, to, n) {
  exp(seq(log(from), log(to), length.out = n))
}
