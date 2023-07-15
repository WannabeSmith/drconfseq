library(drconfseq)
library(parallel)
library(reticulate)
## use_virtualenv("./venv_drconfseq")
library(here)
library(rje) # for expit
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
library(gsDesign)
library(rpact)

n <- 10000
N <- n * 2
start_time <- 500
d <- 3
t_opt <- 1000
alpha <- 0.1
prop_score_min_personalized <- 0.1
prop_score_max_personalized <- 1 - prop_score_min_personalized
num_repeats <- 1000
## num_repeats_width <- 5
beta_mu <- c(1, -1, -2, 0.5)

times <- unique(round(logseq(start_time, n, 100)))

ATE <- 0.2

prop_score_settings <-
  list(
    "bernoulli_expt" =
      list(
        "prop_score_true" =
          function(y, X, newX) {
            return(1 / 2)
          },
        "prop_score_min" = 1 / 2,
        "prop_score_max" = 1 / 2
      ),
    "personalized_expt" =
      list(
        "prop_score_true" =
          function(y, X, newX) {
            newX_matrix <- model.matrix(~., data = newX)
            regs <- apply(newX_matrix, MARGIN = 1, FUN = reg_true)
            untruncated_pscores <- expit(newX_matrix[, 1] / 2)
            return(pmin(prop_score_max, pmax(prop_score_min, untruncated_pscores)))
          },
        "prop_score_min" = prop_score_min_personalized,
        "prop_score_max" = prop_score_max_personalized
      )
  )

prop_score_setting_name <- "bernoulli_expt"
prop_score_setting <- prop_score_settings[[prop_score_setting_name]]
prop_score_true <- prop_score_setting[["prop_score_true"]]


reg_true <- function(x) {
  rbinom(1, 1,
    prob = 1 / 4 +
      1 / 2 * expit(beta_mu %*% c(
        x[1],
        sqrt(abs(x[2])),
        sin(x[3]),
        abs(x[4])
      ))
  )
}


data_generator_fn <- function() {
  X_mtx <- cbind(1, matrix(rnorm(N * d), nrow = N))
  X_data <- X_mtx %>%
    as.data.frame() %>%
    mutate(V1 = NULL)

  reg_observed <- apply(X_mtx, MARGIN = 1, FUN = reg_true)
  p <- prop_score_true(y = NULL, X = NULL, newX = X_data)
  treatment <- rbinom(N, 1, p)
  y <- rbinom(N, 1, prob = 0.6 * treatment + 0.4 * (1 - treatment))


  # Get GLM superlearner
  glm_reg_1 <- get_SL_fn(SL.library = "SL.mean")
  glm_reg_0 <- glm_reg_1

  stopifnot(n %% 2 == 0)
  train_idx <- sample(c(rep(1, N / 2), rep(0, N / 2)))

  pseudo_outcomes <- pseudo_outcome_sequential(
    y = y,
    X = X_data,
    treatment = treatment,
    regression_fn_1 = glm_reg_1,
    regression_fn_0 = glm_reg_0,
    propensity_score_fn = prop_score_true,
    train_idx = train_idx,
    times = NULL,
    cross_fit = FALSE
  )

  return(pseudo_outcomes[[as.character(N)]])
}

acs_fn <- function(y) {
  asymptotic_confseq(
    x = y,
    t_opt = t_opt,
    alpha = alpha,
    return_all_times = TRUE
  )
}

grpseq_general <- function(y, typeOfDesign) {
  kMax <- 10

  design <- getDesignGroupSequential(kMax = kMax, alpha = alpha, sided = 2, typeOfDesign = typeOfDesign)

  grpseq_times <- unique(round(logseq(from = start_time, to = n, n = kMax)))

  grpseq_means <- cumul_mean(y)[grpseq_times]
  grpseq_vars <- cumul_var(y)[grpseq_times]

  ## data <- getDataSet(
  ##   design = design,
  ##   n = grpseq_times,
  ##   means = grpseq_means,
  ##   stDevs = sqrt(grpseq_vars)
  ## )

  ## grpseq_cis <- getRepeatedConfidenceIntervals(design, dataInput = data)

  grpseq_margin <- design$criticalValues * sqrt(grpseq_vars) / sqrt(grpseq_times)

  grpseq_l <- grpseq_means - grpseq_margin
  grpseq_u <- grpseq_means + grpseq_margin

  grpseq_l_full <- rep(-Inf, n)
  grpseq_u_full <- rep(Inf, n)

  ## names(grpseq_l) <- 1:n
  ## names(grpseq_u_full) <- 1:n

  grpseq_l_full[grpseq_times] <- grpseq_l
  grpseq_u_full[grpseq_times] <- grpseq_u

  for(i in 1:length(grpseq_l_full)) {
    if (grpseq_l_full[i] == -Inf && i > 1) {
      grpseq_l_full[i] <- grpseq_l_full[i - 1]
    }
  }

  for(i in 1:length(grpseq_u_full)) {
    if (grpseq_u_full[i] == -Inf && i > 1) {
      grpseq_u_full[i] <- grpseq_u_full[i - 1]
    }
  }

  # Replace all zeros with the last non-zero value in lower confidence intervals
  noninf_indices <- which(grpseq_l_full != -Inf)
  noninf_values <- grpseq_l_full[noninf_indices]
  grpseq_l_full <- rep(noninf_values, diff(c(0, noninf_indices)))

  # Replace all zeros with the last non-zero value in upper confidence intervals
  noninf_indices <- which(grpseq_u_full != Inf)
  noninf_values <- grpseq_u_full[noninf_indices]
  grpseq_u_full <- rep(noninf_values, diff(c(0, noninf_indices)))

  return(list("l" = grpseq_l_full, "u" = grpseq_u_full))
}

grpseq_pocock_fn <- function(y) {
  grpseq_general(y, typeOfDesign = "P")
}


grpseq_obrien_fleming_fn <- function(y) {
  grpseq_general(y, typeOfDesign = "OF")
}

clt_fn <-
  function(y) {
    naive_confidence_intervals(x = y, alpha = alpha, return_all_times = TRUE)
  }


print("Simulating time-uniform confidence sequence miscoverage...")
acs_miscoverage <-
  get_cumul_miscoverage_rate(
    data_generator_fn = data_generator_fn,
    conf_set_fn = acs_fn,
    times = start_time:n,
    num_repeats = num_repeats,
    mu = ATE,
    n_cores = parallel::detectCores()
  )

print("Simulating fixed-time confidence interval miscoverage...")
clt_miscoverage <-
  get_cumul_miscoverage_rate(
    data_generator_fn = data_generator_fn,
    conf_set_fn = clt_fn,
    times = start_time:n,
    num_repeats = num_repeats,
    mu = ATE,
    n_cores = parallel::detectCores()
  )

print("Simulating pocock miscoverage...")
pocock_miscoverage <-
  get_cumul_miscoverage_rate(
    data_generator_fn = data_generator_fn,
    conf_set_fn = grpseq_pocock_fn,
    times = start_time:n,
    num_repeats = num_repeats,
    mu = ATE,
    n_cores = parallel::detectCores()
  )

print("Simulating O'Brien-Fleming miscoverage...")
obrien_fleming_miscoverage <-
  get_cumul_miscoverage_rate(
    data_generator_fn = data_generator_fn,
    conf_set_fn = grpseq_obrien_fleming_fn,
    times = start_time:n,
    num_repeats = num_repeats,
    mu = ATE,
    n_cores = parallel::detectCores()
  )

## print("Simulating sub-exponential miscoverage...")
## subexponential_miscoverage <-
##   get_cumul_miscoverage_rate(
##     data_generator_fn = data_generator_fn,
##     conf_set_fn = subexponential_fn,
##     times = start_time:n,
##     num_repeats = num_repeats,
##     mu = ATE,
##     n_cores = parallel::detectCores()
##   )

print("Simulating asympcs width...")
acs_width <- get_avg_width(
  data_generator_fn = data_generator_fn,
  conf_set_fn = acs_fn,
  times = start_time:n,
  num_repeats = num_repeats,
  n_cores = parallel::detectCores()
)

print("Simulating clt width...")
clt_width <- get_avg_width(
  data_generator_fn = data_generator_fn,
  conf_set_fn = clt_fn,
  times = start_time:n,
  num_repeats = num_repeats,
  n_cores = parallel::detectCores()
)

print("Simulating pocock width...")
pocock_width <-
  get_avg_width(
    data_generator_fn = data_generator_fn,
    conf_set_fn = grpseq_pocock_fn,
    times = start_time:n,
    num_repeats = num_repeats,
    n_cores = parallel::detectCores()
  )

print("Simulating O'Brien-Fleming width...")
obrien_fleming_width <-
  get_avg_width(
    data_generator_fn = data_generator_fn,
    conf_set_fn = grpseq_obrien_fleming_fn,
    times = start_time:n,
    num_repeats = num_repeats,
    n_cores = parallel::detectCores()
  )

## print("Simulating subGaussian width...")
## subGaussian_width <- get_avg_width(
##   data_generator_fn = data_generator_fn,
##   conf_set_fn = subGaussian_fn,
##   times = start_time:n,
##   num_repeats = num_repeats,
##   n_cores = parallel::detectCores()
## )

## print("Simulating subexponential width...")
## subexponential_width <- get_avg_width(
##   data_generator_fn = data_generator_fn,
##   conf_set_fn = subexponential_fn,
##   times = start_time:n,
##   num_repeats = num_repeats,
##   n_cores = parallel::detectCores()
## )





## y <- data_generator_fn()
## acs <- acs_fn(y)
## clt <- clt_fn(y)
## script.dir <- dirname(sys.frame(1)$ofile)

r_data_dir <- "paper_plots/simulations/group_sequential/"

save(
  acs_miscoverage,
  clt_miscoverage,
  pocock_miscoverage,
  obrien_fleming_miscoverage,
  times,
  acs_width,
  clt_width,
  pocock_width,
  obrien_fleming_width,
  file = paste(r_data_dir,
    "grpseq_widths",
    "_",
    prop_score_setting_name,
    ".RData",
    sep = ""
  )
)
