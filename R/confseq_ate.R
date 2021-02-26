require(SuperLearner)
require(rpart)
require(parallel)

#' Confidence sequence for the average treatment effect.
#'
#' Function signature: (real, data.frame, int,
#'                      (real, data.frame, data.frame -> real),
#'                      (real, data.frame, data.frame -> real),
#'                      (int, data.frame, data.frame -> real),
#'                      real, real, int, int) -> data.frame
#' @importFrom parallel mclapply
#' @param y The measured outcome (a real vector).
#' @param X Measured covariates (an nxd real matrix where `n = length(y)`)
#' @param treatment Whether the subject received treatment
#'                  (a boolean vector or 0/1-valued integer vector).
#' @param regression_fn_1 A function which predicts outcomes for those who took
#'                        treatment. This function takes three arguments: `y`,
#'                        `X`, and `newX`, the training outcome, training
#'                        covariates, and evaluation covariates, respectively.
#'                        The function outputs the predicted response given the
#'                        evaluation covariates, `newX`.
#' @param regression_fn_0 The same as `regression_fn_1` but for those who did
#'                        not receive treatment. If left NULL, this will be
#'                        set to the same function as `regression_fn_1`.
#' @param propensity_score A function which predicts the propensity score for
#'                         each subject. Similar to `regression_fn_1`, this
#'                         function takes three arguments: `y`, `X`, and `newX`,
#'                         the training treatment indicator (1 if treatment, 0
#'                         if control), the training covariates, and the
#'                         evaluation covariates. The function outputs the
#'                         predicted propensity score given the evaluation
#'                         covariates, `newX`.
#' @param train_idx The indices indicating the training split for the sample
#'                  splitting algorithm. If left NULL, the training index will
#'                  be assigned randomly with probability 1/2.
#' @param t_opt Time for which the CS should be tightest
#' @param alpha Confidence level between 0 and 1 (real)
#' @param times The times for which the doubly-robust variables should be
#'              calculated. Can be a vector of times (an integer vector) or
#'              a single time (integer). If left NULL, the variables will
#'              only be computed at time n.
#' @param n_cores The number of cores to use for parallelization.
#' @return Data frame containing the lower and upper confidence sequences.
#' @export
confseq_ate <- function(y, X, treatment,
                        regression_fn_1,
                        regression_fn_0 = NULL,
                        propensity_score_fn,
                        t_opt,
                        train_idx = NULL,
                        alpha = 0.05,
                        times = NULL,
                        n_cores = 1,
                        cross_fit = TRUE)
{
  # Get the 'doubly-robust ATE' variables as a list for each time in times
  # If times is NULL, this will just compute the variables for time n.
  pseudo_outcome_list <-
    pseudo_outcome_sequential(y = y, X = X,
                              treatment = treatment,
                              regression_fn_1 = regression_fn_1,
                              regression_fn_0 = regression_fn_0,
                              propensity_score_fn = propensity_score_fn,
                              times = times,
                              n_cores = n_cores,
                              cross_fit = cross_fit)

  # Compute the confidence sequence at each of these times
  confseq_list <- mclapply(pseudo_outcome_list, function(pseudo_outcomes){
    acs <- asymptotic_confseq(pseudo_outcomes,
                              t_opt = t_opt,
                              alpha = alpha,
                              return_all_times = FALSE)
    return(c(acs$l, acs$u))
  }, mc.cores = n_cores)

  confseq <- data.frame(do.call(rbind, confseq_list))
  colnames(confseq) <- c('l', 'u')

  return(confseq)
}

#' Unadjusted estimator
#'
#' @export
confseq_ate_unadjusted <- function(y, treatment,
                                   propensity_score,
                                   t_opt,
                                   alpha = 0.05,
                                   times = NULL)
{
  pseudo_outcome <-
    pseudo_outcome_abstract(y = y, reg_1 = 0,
                            reg_0 = 0,
                            propensity_score = propensity_score,
                            treatment = treatment)

  cs <- asymptotic_confseq(pseudo_outcome, t_opt = t_opt,
                            alpha = alpha)

  df <- data.frame(l = cs$l, u = cs$u,
                   row.names = 1:length(y))

  if(all(!is.na(times)))
  {
    df <- df[times,]
  }

  return(df)
}

