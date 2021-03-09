library(sequential.causal)
library(parallel)
library(pracma)

ATE <- 1
n = 10000
d = 3
X <- cbind(1, matrix(rnorm(n*d), nrow = n))

beta_mu <- c(1, -1, -2, 3)

reg_true <- function(x)
{
  beta_mu %*% c(x[1], x[2]^2, sin(x[3]), abs(x[4]))
}

prop_score_true <- function(x)
{
  logodds <- reg_true(x)
  pi = exp(logodds)/(1 + exp(logodds))
  # Ensure bounded away from 0 and 1
  pi = pi*0.6 + 0.2
}

reg_observed <- apply(X, MARGIN=1, FUN=reg_true)
p <- apply(X, MARGIN=1, FUN=prop_score_true)
treatment <- rbinom(n, 1, p)
y <- reg_observed + treatment*ATE + rt(n, df=5)

drate_variables_oracle <-
  pseudo_outcome_abstract(y = y, reg_1 = reg_observed,
                          reg_0 = reg_observed, propensity_score = p,
                          treatment = treatment)



# Get SuperLearner prediction function for $\mu^1$.
# Using default ML algorithm choices
sl_reg_1 <- get_SL_fn()

# Do the same for $\mu^0$.
sl_reg_0 <- get_SL_fn()

# Get SuperLearner prediction function for $\pi$
pi_fn <- get_SL_fn(family = binomial)

# Get GLM superlearner
glm_reg_1 = get_SL_fn(SL.library = "SL.glm")


times <- unique(round(logseq(250, 10000, n = 30)))
alpha <- 0.05
n_cores <- detectCores()

# Split the sample (if we don't do this explicitly,
# confseq_ate or drate_variables_sequential can automatically)
train_idx <- rbinom(n, p = 0.5, size = 1) == 1

# # Doing this the more abstract way so that we can plot confidence
# # sequences alongside confidence intervals without re-doing the analysis.
# # In practice this is simpler by just using confseq_ate.
# drate_variables_SL_list <-
#   drate_variables_sequential(y = y, X = X,
#                              treatment = treatment,
#                              regression_fn_1 = sl_reg_1,
#                              propensity_score_fn = function(y, X, newX){1/2},
#                              train_idx = train_idx,
#                              times = times,
#                              n_cores = n_cores)
#
# drate_variables_glm_list <-
#   drate_variables_sequential(y = y, X = X,
#                              treatment = treatment,
#                              regression_fn_1 = glm_reg_1,
#                              propensity_score_fn = function(y, X, newX){1/2},
#                              train_idx = train_idx,
#                              times = times,
#                              n_cores = n_cores)
#
# drate_variables_unadj_list <-
#   drate_variables_sequential(y = y, X = X,
#                              treatment = treatment,
#                              regression_fn_1 = function(y, X, newX){0},
#                              propensity_score_fn = function(y, X, newX){1/2},
#                              train_idx = train_idx,
#                              times = times,
#                              n_cores = n_cores)

confseq_SL <- confseq_ate(y, X, treatment, regression_fn_1 = sl_reg_1,
                          regression_fn_0 = sl_reg_1,
                          propensity_score_fn = pi_fn,
                          train_idx = train_idx, t_opt = 250, alpha=alpha,
                          times=times, n_cores = n_cores, cross_fit = TRUE)
confseq_glm <- confseq_ate(y, X, treatment, regression_fn_1 = glm_reg_1,
                           regression_fn_0 = glm_reg_1,
                           propensity_score_fn = get_SL_fn(SL.library="SL.glm",
                                                           family=binomial()),
                           train_idx = train_idx, t_opt = 250, alpha=alpha,
                           times=times, n_cores = n_cores, cross_fit = TRUE)
# confseq_unadj <- confseq_ate(y, X, treatment,
#                              regression_fn_1 = function(y, X, newX){0},
#                              propensity_score_fn = function(y, X, newX){1/2},
#                              train_idx = train_idx, t_opt = 250, alpha=alpha,
#                              times=times, n_cores = n_cores)


confseq_unadj <- confseq_ate_unadjusted(y = y, treatment = treatment,
                                        propensity_score = 1/2,
                                        t_opt = 250, alpha = alpha,
                                        times = times)



r_data_dir <- "./"
save(confseq_SL, confseq_glm, confseq_unadj, times, ATE,
     file=paste(r_data_dir, 'observational_study_SL.RData', sep=""))

