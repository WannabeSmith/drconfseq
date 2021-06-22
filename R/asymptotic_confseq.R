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
std_conjmix_margin <- function(t, rho2, alpha=0.05)
{
  return(
    sqrt(
      2*(t*rho2 + 1)*
        log(sqrt(t*rho2 + 1) / alpha) / (t^2 * rho2)
    )
  )
}

#' Get the best value of $rho^2$ for a time t_opt (exact optimization)
#'
#' Function signature: (int, real) -> real
#'
#' @param t_opt The time for which $rho^2$
#'              should be optimized (a positive integer)
#' @param alpha The significance level (a (0, 1)-valued real).
#' @return The corresponding value of $rho^2$ (positive real)
#' @export
best_rho2_exact <- function(t_opt, alpha_opt=0.05)
{
  (-lambertWm1(-alpha_opt^2 * exp(alpha^2 - 1)) - 1) / t_opt
}

lambertWm1_approx <- function(x)
{
  # Based on Taylor approximation of Lambert W function
  # https://en.wikipedia.org/wiki/Lambert_W_function#Asymptotic_expansions
  # https://cs.uwaterloo.ca/research/tr/1993/03/W.pdf
  log(-x) - log(-log(-x))
}

#' Get the best value of $rho^2$ for a time t_opt (approximate optimization)
#'
#' Function signature: (int, real) -> real
#'
#' @param t_opt The time for which $rho^2$
#'              should be optimized (a positive integer)
#' @param alpha The significance level (a (0, 1)-valued real).
#' @return The corresponding value of $rho^2$ (positive real)
#' @export
best_rho2_approx <- function(t_opt, alpha_opt=0.05)
{
  (-lambertWm1_approx(-alpha_opt^2 * exp(alpha^2 - 1)) - 1) / t_opt
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
    std_margin <- std_LIL_margin(t=t, alpha=alpha)
  } else
  {
    rho2 <- best_rho2_exact(t_opt=t_opt, alpha_opt=alpha)
    std_margin <- std_conjmix_margin(t=t, rho2=rho2, alpha=alpha)
  }

  margin <- sqrt(var)*std_margin

  # When the margin is NA, such as on the first observation, just return
  # infinity
  margin[is.na(margin)] = Inf

  return(list('l' = mu_hat_t - margin,
              'u' = mu_hat_t + margin))
}

#' Plot ratio of confidence sequence to confidence interval widths
#'
#' @importFrom purrr map reduce
#' @importFrom ggplot2 ggplot geom_line geom_point aes guides
#'             ylab theme_minimal theme scale_x_log10 annotation_logticks
#'             guide_legend element_text
#' @param t_opts The times for which to optimize the confidence sequence
#'               (vector of positive integers)
#' @param t The times to plot the confidence sequence
#'          (vector of positive integers)
#' @param alpha The significance level, e.g. set alpha=0.05 for 95% coverage
#'              (real number between 0 and 1)
#' @param log_scale Should the plot be returned on a log-scale? (boolean)
#'
#' @return A ggplot2 plot object
#' @export
plot_cs_shape <- function(t_opts, t, alpha = 0.05, log_scale = FALSE)
{
  alpha <- 0.05
  if(log_scale)
  {
    shape_points <- unique(round(logseq(min(t), max(t), n = 6)))
  } else
  {
    shape_points <- unique(round(seq(min(t), max(t), length.out = 6)))
  }

  plt_data_lines <- t_opts %>% map(function(t_opt){
    acs <- std_conjmix_margin(t = t,
                              rho2 = best_rho2_exact(t_opt),
                              alpha=alpha)
    naive <- naive_std_margin(t = t, alpha = alpha)
    data.frame(Time = t,
               Ratio = naive / acs,
               t_opt = as.factor(t_opt))
  }) %>% reduce(rbind)

  plt_data_points <- plt_data_lines[plt_data_lines$Time %in% shape_points, ]

  plt <-
    ggplot(plt_data_lines) +
    geom_line(aes(x = Time, y = Ratio, color = t_opt)) +
    geom_point(data = plt_data_points,
               aes(x = Time, y = Ratio, color = t_opt, shape = t_opt), size=2) +
    guides(color = guide_legend(title = "t optimized"),
           shape = guide_legend(title = "t optimized")) +

    ylab("Ratio of confidence interval widths\nto confidence sequence widths") +
    theme_minimal() +
    theme(text = element_text(family="serif"), legend.position=c(0.75, 0.3))

  if (log_scale)
  {
    plt <- plt +
      scale_x_log10(breaks = scales::trans_breaks("log10",
                                                  function(x) 10^x),
                    labels = scales::trans_format("log10",
                                                  scales::math_format(10^.x))) +
      annotation_logticks(colour = "grey", side = "b")
  }

  plt
  return(plt)
}

