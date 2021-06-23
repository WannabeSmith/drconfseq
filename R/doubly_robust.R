# This file contains functions related to causal inference,
# double robustness, double machine learning, targeted learning, etc.

#' Doubly-robust pseudo-outcome for binary treatment given nuisance function estimates.
#'
#' @param y The measured outcome (a real vector).
#' @param reg_1 The predicted outcome in the treatment group
#'                     (a real vector).
#' @param reg_0 The predicted outcome in the control group
#'                     (a real vector).
#' @param propensity_score The predicted (or known) propensity score
#'                         (a (0, 1)-valued vector)
#' @param treatment Whether the subject received treatment
#'                  (a boolean vector or 0/1-valued integer vector).
#' @return The doubly-robust variables at each time
#'         (a vector of real numbers)
#' @export
pseudo_outcome_abstract <- function(y, reg_1,
                                    reg_0,
                                    propensity_score,
                                    treatment)
{
  # Length of the experiment
  N <- length(y)

  variables <-
    # difference of regression functions
    (reg_1 - reg_0) +
    # treatments by propensity score
    (treatment/propensity_score - (1-treatment)/(1-propensity_score))*
    # y minus regression for given treatment
    (y - treatment*reg_1 - (1-treatment)*reg_0)

  return(variables)
}

#' Doubly-robust variables at many timepoints for binary treatment given nuisance function estimators.
#'
#' @importFrom stats rbinom
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
#' @param propensity_score_fn A function which predicts the propensity score for
#'                            each subject. Similar to `regression_fn_1`, this
#'                            function takes three arguments: `y`, `X`, and
#'                            `newX`,
#'                            the training treatment indicator (1 if treatment,
#'                            0 if control), the training covariates, and the
#'                            evaluation covariates. The function outputs the
#'                            predicted propensity score given the evaluation
#'                            covariates, `newX`.
#' @param train_idx The indices indicating the training split for the sample
#'                  splitting algorithm. If left NULL, the training index will
#'                  be assigned randomly with probability 1/2.
#' @param times The times for which the doubly-robust variables should be
#'              calculated. Can be a vector of times (an integer vector) or
#'              a single time (integer). If left NULL, the variables will
#'              only be computed at time n.
#' @param n_cores The number of cores to use for parallelization.
#' @param cross_fit Should cross-fitting be used? (boolean)
#' @return A list of doubly-robust variable vectors, one vector for each time
#'         requested in `times`.
#' @export
pseudo_outcome_sequential <- function(y, X, treatment,
                                      regression_fn_1,
                                      regression_fn_0 = NULL,
                                      propensity_score_fn,
                                      train_idx = NULL,
                                      times = NULL,
                                      n_cores = 1, cross_fit = FALSE)
{
  if (is.null(regression_fn_0)) regression_fn_0 <- regression_fn_1

  # Get sample splitting index if not provided by user
  if (is.null(train_idx)) train_idx <- rbinom(length(y), 1, 0.5)
  train_idx <- train_idx == TRUE

  eval_idx <- 1-train_idx == TRUE

  if (any(is.null(times)))
  {
    warning('"times" was left as null. Computing only at time n.')
    times = length(y)
  }

  if(cross_fit)
  {
    train_indices <- list(train_idx, eval_idx)
  } else
  {
    train_indices <- list(train_idx)
  }

  variables_list <- mclapply(times, function(time){
    # This unlist(lapply(...)) looks a bit strange, but it is
    # necessary for when the user wishes to perform cross-fitting.
    # See above if/else statement.
    unlist(lapply(train_indices, function(train_idx){
      train_idx_t <- train_idx == TRUE
      train_idx_t[(time+1):length(train_idx_t)] = FALSE
      eval_idx_t <- 1 - train_idx == TRUE
      eval_idx_t[(time+1):length(eval_idx_t)] = FALSE
      # y
      y_train <- y[train_idx_t]
      y_eval <- y[eval_idx_t]
      # X
      X_train <- X[train_idx_t, ]
      X_eval <- X[eval_idx_t, ]
      # treatment
      treatment_train <- treatment[train_idx_t]
      treatment_eval <- treatment[eval_idx_t]

      reg_1_est <- regression_fn_1(y = y_train[treatment_train==1],
                                   X = X_train[treatment_train==1, ],
                                   newX = X_eval)
      message(paste('Finished fitting regression model for treated at time',
                    time))
      reg_0_est <- regression_fn_0(y = y_train[treatment_train==0],
                                   X = X_train[treatment_train==0, ],
                                   newX = X_eval)
      message(paste('Finished fitting regression model for untreated at time',
                    time))
      pi_est <- propensity_score_fn(y = as.numeric(treatment_train),
                                    X = X_train,
                                    newX = X_eval)
      message(paste('Finished fitting prop. score function at time',
                    time))

      variables <- pseudo_outcome_abstract(y = y_eval,
                                           reg_1 = reg_1_est,
                                           reg_0 = reg_0_est,
                                           propensity_score = pi_est,
                                           treatment = treatment_eval)
      return(variables)
    }))
  }, mc.cores = n_cores)

  names(variables_list) <- times

  return(variables_list)
}

#' Generic doubly-robust estimator for binary treatment
#'
#' @param y The measured outcome (a real number).
#' @param reg_1 The predicted outcome in the treatment group
#'                     (a real vector).
#' @param reg_0 The predicted outcome in the control group
#'                     (a real vector).
#' @param propensity_score The predicted (or known) propensity score
#'                         (a (0, 1)-valued vector)
#' @param treatment Whether the subject received treatment
#'                  (a boolean vector or 0/1-valued integer vector).
#' @return The doubly-robust estimator at each time (a vector of real numbers)
#' @export
pseudo_outcome_estimator <- function(y, reg_1,
                                     reg_0,
                                     propensity_score,
                                     treatment)
{
  # get doubly-robust observations
  observations <- pseudo_outcome_abstract(y, reg_1,
                                          reg_0, propensity_score,
                                          treatment)
  return(cumul_mean(observations))
}

#' Variance of the doubly robust observations
#'
#' @param y The measured outcome (a real number).
#' @param reg_1 The predicted outcome in the treatment group
#'                     (a real vector).
#' @param reg_0 The predicted outcome in the control group
#'                     (a real vector).
#' @param propensity_score The predicted (or known) propensity score
#'                         (a (0, 1)-valued vector)
#' @param treatment Whether the subject received treatment
#'                  (a boolean vector or 0/1-valued integer vector).
#' @return The variance of the doubly-robust observations
#'         (a vector of real numbers)
#' @export
pseudo_outcome_variance <- function(y, reg_1, reg_0,
                                    propensity_score,
                                    treatment)
{
  # get doubly-robust observations
  observations <- pseudo_outcome_abstract(y, reg_1,
                                          reg_0, propensity_score,
                                          treatment)

  return(cumul_var(observations))
}


#' Get the Super Learner function given a library of SuperLearner algorithms.
#'
#' @importFrom SuperLearner SuperLearner
#' @importFrom stats gaussian
#' @param SL.library Character vector containing SuperLearner algorithm names.
#' @param family The family for the algorithm. If estimating the propensity
#'               score, set `family = binomial`, but if estimating the
#'               regression function, set `family = gaussian()`.
#' @return Function which takes three arguments `y`, `X`, and `newX`, the
#'         training outcome, training covariates, and evaluation covariates,
#'         respectively. This function takes `y` and `X` to train the
#'         algorithm and uses the `SuperLearner` library to predict the outcome
#'         for the covariates `newX`.
#' @export
get_SL_fn <- function(SL.library=c("SL.earth","SL.gam",
                                   "SL.glm","SL.glmnet",
                                   "SL.glm.interaction",
                                   "SL.ranger"),
                      family=gaussian())
{
  SL_fn <- function(y, X, newX)
  {
    SuperLearner(Y = y, X = X, newX = newX,
                 SL.library = SL.library,
                 family = family)$SL.predict
  }

  return(SL_fn)
}

