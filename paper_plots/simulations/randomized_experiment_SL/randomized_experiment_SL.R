library(drconfseq)
library(parallel)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)

ATE <- 1
n = 10000
start_time <- 500
d = 3
X_mtx <- cbind(1, matrix(rnorm(n * d), nrow = n))
X_data <- X_mtx %>%
  as.data.frame() %>%
  mutate(V1 = NULL)

beta_mu <- c(1,-1,-2, 3)

reg_true <- function(x)
{
  beta_mu %*% c(x[1], x[2] ^ 2, sin(x[3]), abs(x[4]))
}

prop_score_true <- function(x) {
  1 / 2
}

reg_observed <- apply(X_mtx, MARGIN = 1, FUN = reg_true)
p <- apply(X_mtx, MARGIN = 1, FUN = prop_score_true)
treatment <- rbinom(n, 1, p)
y <- reg_observed + treatment * ATE + rt(n, df = 5)

# Get SuperLearner prediction function for $\mu^1$.
# Using default ML algorithm choices
sl_reg_1 <- get_SL_fn()

# Do the same for $\mu^0$.
sl_reg_0 <- sl_reg_1

# Get SuperLearner prediction function for $\pi$
pi_fn <- get_SL_fn(family = binomial)

# Get GLM superlearner
glm_reg_1 = get_SL_fn(SL.library = "SL.glm")
glm_reg_0 <- glm_reg_1


times <- unique(round(logseq(start_time, 10000, n = 30)))
alpha <- 0.1
n_cores <- detectCores()

# Split the sample (if we don't do this explicitly,
# confseq_ate or drate_variables_sequential can automatically)
train_idx <- rbinom(n, p = 0.5, size = 1) == 1

confseq_SL <-
  confseq_ate(
    y,
    X_data,
    treatment,
    regression_fn_1 = sl_reg_1,
    regression_fn_0 = sl_reg_0,
    propensity_score_fn = function(y, X, newX) {
      1 / 2
    },
    train_idx = train_idx,
    t_opt = 1000,
    alpha = alpha,
    times = times,
    n_cores = n_cores,
    cross_fit = TRUE
  )
confseq_glm <-
  confseq_ate(
    y,
    X_data,
    treatment,
    regression_fn_1 = glm_reg_1,
    regression_fn_0 = glm_reg_0,
    propensity_score_fn = function(y, X, newX) {
      1 / 2
    },
    train_idx = train_idx,
    t_opt = 1000,
    alpha = alpha,
    times = times,
    n_cores = n_cores,
    cross_fit = TRUE
  )
confseq_unadj <-
  confseq_ate_unadjusted(
    y = y,
    treatment = treatment,
    propensity_score = 1 / 2,
    t_opt = 1000,
    alpha = alpha,
    times = times
  )

r_data_dir <- "./"
save(
  confseq_SL,
  confseq_glm,
  confseq_unadj,
  times,
  ATE,
  file = paste(r_data_dir, 'randomized_experiment_SL.RData', sep = "")
)
