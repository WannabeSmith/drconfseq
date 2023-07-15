library(drconfseq)
library(parallel)
library(reticulate)
## use_virtualenv("./venv_drconfseq")
library(here)
library(rje) # for expit
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)

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

for (prop_score_setting_name in names(prop_score_settings)) {
  prop_score_setting <- prop_score_settings[[prop_score_setting_name]]

  prop_score_true <- prop_score_setting[["prop_score_true"]]
  prop_score_min <- prop_score_setting[["prop_score_min"]]
  prop_score_max <- prop_score_setting[["prop_score_max"]]

  lower_bd <- -1 / prop_score_min - 1
  upper_bd <- 1 / prop_score_min + 1


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

  clt_fn <-
    function(y) {
      naive_confidence_intervals(x = y, alpha = alpha, return_all_times = TRUE)
    }


  subGaussian_fn <-
    function(y) {
      robbins_subGaussian_cs(
        x = y,
        alpha = alpha,
        t_opt = t_opt,
        sigma2 = (upper_bd - lower_bd)^2 / 4
      )
    }

  subexponential_fn <- function(y) {
    howard_subexponential_cs(
      x = y,
      alpha = alpha,
      v_opt = 1 / 4 * t_opt,
      running_intersection = FALSE,
      lower_bd = lower_bd,
      upper_bd = upper_bd
    )
  }

  ## betting_fn <-
  ##   function(y) {
  ##     betting_cs(x = y, alpha = alpha)
  ##   }


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

  print("Simulating sub-Gaussian miscoverage...")
  subGaussian_miscoverage <-
    get_cumul_miscoverage_rate(
      data_generator_fn = data_generator_fn,
      conf_set_fn = subGaussian_fn,
      times = start_time:n,
      num_repeats = num_repeats,
      mu = ATE,
      n_cores = parallel::detectCores()
    )


  print("Simulating sub-exponential miscoverage...")
  subexponential_miscoverage <-
    get_cumul_miscoverage_rate(
      data_generator_fn = data_generator_fn,
      conf_set_fn = subexponential_fn,
      times = start_time:n,
      num_repeats = num_repeats,
      mu = ATE,
      n_cores = parallel::detectCores()
    )
  ## print("Simulating time-uniform betting CS miscoverage...")
  ## betting_miscoverage <-
  ##   get_cumul_miscoverage_rate(
  ##     data_generator_fn = data_generator_fn,
  ##     conf_set_fn = betting_fn,
  ##     times = start_time:n,
  ##     num_repeats = num_repeats,
  ##     mu = p,
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

  print("Simulating subGaussian width...")
  subGaussian_width <- get_avg_width(
    data_generator_fn = data_generator_fn,
    conf_set_fn = subGaussian_fn,
    times = start_time:n,
    num_repeats = num_repeats,
    n_cores = parallel::detectCores()
  )

  print("Simulating subexponential width...")
  subexponential_width <- get_avg_width(
    data_generator_fn = data_generator_fn,
    conf_set_fn = subexponential_fn,
    times = start_time:n,
    num_repeats = num_repeats,
    n_cores = parallel::detectCores()
  )

  ## y <- data_generator_fn()
  ## acs <- acs_fn(y)
  ## clt <- clt_fn(y)
  ## script.dir <- dirname(sys.frame(1)$ofile)
  r_data_dir <- "paper_plots/simulations/widths_coverage/"

  save(
    acs_width,
    clt_width,
    subGaussian_width,
    subexponential_width,
    acs_miscoverage,
    clt_miscoverage,
    subGaussian_miscoverage,
    subexponential_miscoverage,
    ## acs,
    ## clt,
    ## times,
    ## start_time,
    ## n,
    ## file = "widths_coverage_ate.RData"
    file = paste(r_data_dir,
      "widths_coverage_ate",
      "_",
      prop_score_setting_name,
      ".RData",
      sep = ""
    )
  )
}
