# This file contains functions related to asymptotic confidence sequences.

#' Conjugate mixture standard (unit-variance) margin
#'
#' Function signature: (int, real, real) -> real
#'
#' @param t The times at which to product the margin
#'          (a positive integer vector).
#' @param rho2 The tuning parameter (a positive real).
#' @param alpha The significance level (a (0, 1)-valued real).
#' @return The standard conjugate mixture margins (a real vector).
#' @export
std_conjmix_margin <- function(t, rho2, alpha=0.05/2)
{
  return(
    sqrt(
      2*(t*rho2 + 1)*
        log(sqrt(t*rho2 + 1) / alpha) / (t^2 * rho2)
    )
  )
}

#' Get the best value of $\rho^2$ for a time t_opt (exact optimization)
#'
#' Function signature: (int, real) -> real
#'
#' @param t_opt The time for which $\rho^2$
#'              should be optimized (a positive integer)
#' @param alpha The significance level (a (0, 1)-valued real).
#' @return The corresponding value of $\rho^2$ (positive real)
#' @export
best_rho2_exact <- function(t_opt, alpha_opt=0.05/2)
{
  # Need to adjust this hard-coded upper search bound of 10
  optimize(function(rho2)
  {std_conjmix_margin(t=t_opt, rho2=rho2, alpha=alpha_opt)},
  interval=c(0, 10))$minimum
}

#' Get the best value of $\rho^2$ for a time t_opt (approximate optimization)
#'
#' Function signature: (int, real) -> real
#'
#' @param t_opt The time for which $\rho^2$
#'              should be optimized (a positive integer)
#' @param alpha The significance level (a (0, 1)-valued real).
#' @return The corresponding value of $\rho^2$ (positive real)
#' @export
best_rho2_approx <- function(t_opt, alpha_opt=0.05/2)
{
  (2*log(1/alpha_opt) +
     log(1 + 2*log(1/alpha_opt))) / t_opt
}

#' LIL standard (unit-variance) margin
#'
#' Function signature: (int, real) -> real
#'
#' @param t The times at which to product the margin
#'          (a positive integer vector).
#' @param alpha The significance level (a (0, 1)-valued real).
#' @return The standard LIL margins (a real vector).
#' @export
std_LIL_margin <- function(t, alpha=0.05/2)
{
  return(
    1.7 *
      sqrt(log(log(2*t)) + 0.72*log(5.2 / alpha)) /
      sqrt(t)
  )
}

#' Asymptotic confidence sequence
#'
#' Function signature: (real, real, real, real, boolean) -> (real, real)
#'
#' @param x The observed data points (a real vector).
#' @param alpha The significance level (a (0, 1)-valued real).
#' @param t_opt The time for which the confidence sequence should be tightest.
#' @param var The known or estimated variance of the observations, `x`.
#'            If left `NULL`, then the empirical variance of `x` will be
#'            taken (positive real-valued vector or `NULL`).
#' @param return_all_times Should the CS be returned at
#'                         every time point? (boolean)
#' @return A list containing the following vectors: \cr
#' \item{l}{The lower confidence sequence (a real vector).}
#' \item{u}{The upper confidence sequence (a real vector).}
#' @export
asymptotic_confseq <- function(x, t_opt, alpha=0.05,
                               var=NULL, LIL=FALSE,
                               return_all_times=TRUE)
{
  # If the user wants results for each time, use time from 1 to n.
  # Otherwise, just use time n.
  if(return_all_times)
  {
    t = seq(1, length(x))
    mu_hat_t = cumul_mean(x)
    if (is.null(var)) var <- cumul_var(x)
  } else
  {
    t = length(x)
    mu_hat_t = mean(x)
    if (is.null(var)) var <- var(x)
  }

  if (LIL)
  {
    std_margin <- std_LIL_margin(t=t, alpha=alpha/2)
  } else
  {
    rho2 <- best_rho2_approx(t_opt=t_opt, alpha_opt=alpha/2)
    std_margin <- std_conjmix_margin(t=t, rho2=rho2, alpha=alpha/2)
  }

  margin <- sqrt(var)*std_margin

  # When the margin is NA, such as on the first observation, just return
  # infinity
  margin[is.na(margin)] = Inf

  return(list('l' = mu_hat_t - margin,
              'u' = mu_hat_t + margin))
}

