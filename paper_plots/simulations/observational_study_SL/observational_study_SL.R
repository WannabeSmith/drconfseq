library(drconfseq)
library(parallel)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)

t_opt <- 1000
ATE <- 1
n <- 10000
start_time <- 500
d <- 3
X_mtx <- cbind(1, matrix(rnorm(n * d), nrow = n))
X_data <- X_mtx %>%
  as.data.frame() %>%
  mutate(V1 = NULL)

beta_mu <- c(1, -1, -2, 3)

reg_true <- function(x) {
  beta_mu %*% c(x[1], x[2]^2, sin(x[3]), abs(x[4]))
}

prop_score_true <- function(x) {
  logodds <- reg_true(x)
  pi <- exp(logodds) / (1 + exp(logodds))
  # Ensure bounded away from 0 and 1
  pi <- pi * 0.6 + 0.2
}

reg_observed <- apply(X_mtx, MARGIN = 1, FUN = reg_true)
p <- apply(X_mtx, MARGIN = 1, FUN = prop_score_true)
treatment <- rbinom(n, 1, p)
y <- reg_observed + treatment * ATE + rt(n, df = 5)

drate_variables_oracle <-
  pseudo_outcome_abstract(
    y = y,
    reg_1 = reg_observed,
    reg_0 = reg_observed,
    propensity_score = p,
    treatment = treatment
  )



# Get SuperLearner prediction function for $\mu^1$.
# Using default ML algorithm choices
sl_reg_1 <- get_SL_fn()

# Do the same for $\mu^0$.
sl_reg_0 <- get_SL_fn()

# Get SuperLearner prediction function for $\pi$
pi_fn <- get_SL_fn(family = binomial)

# Get GLM superlearner
glm_reg_1 <- get_SL_fn(SL.library = "SL.glm")


times <- unique(round(logseq(start_time, n, n = 30)))
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
    regression_fn_0 = sl_reg_1,
    propensity_score_fn = pi_fn,
    train_idx = train_idx,
    t_opt = t_opt,
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
    regression_fn_0 = glm_reg_1,
    propensity_score_fn = get_SL_fn(
      SL.library =
        "SL.glm",
      family = binomial()
    ),
    train_idx = train_idx,
    t_opt = t_opt,
    alpha = alpha,
    times = times,
    n_cores = n_cores,
    cross_fit = TRUE
  )
## confseq_unadj <-
##   confseq_ate_unadjusted(
##     y = y,
##     treatment = treatment,
##     propensity_score = NULL,
##     t_opt = 1000,
##     alpha = alpha,
##     times = times
##   )

diff_in_means_influence_fns <- y * treatment / cumul_mean(treatment, regularizer_obs = 1) -
  y * (1 - treatment) / cumul_mean(1 - treatment, regularizer_obs = 1)
confseq_diff_in_means_all <-
  asymptotic_confseq(
    x = diff_in_means_influence_fns,
    t_opt = t_opt,
    alpha = alpha
  )
confseq_diff_in_means <- list(
  "l" = confseq_diff_in_means_all$l[times],
  "u" = confseq_diff_in_means_all$u[times]
)


r_data_dir <- "./"
save(
  confseq_SL,
  confseq_glm,
  confseq_diff_in_means,
  times,
  ATE,
  file = paste(r_data_dir, "observational_study_SL.RData", sep = "")
)
