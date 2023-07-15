library(drconfseq)
library(parallel)
library(reticulate)
## use_virtualenv("./venv_drconfseq")
library(here)
n <- 10000
start_time <- 100
t_opt <- 1000
alpha <- 0.1
p <- 0.5
num_repeats <- 1000

data_generators <- list(
  "bernoulli" = function() {
    rbinom(n, 1, 1 / 2)
  },
  "beta" = function() {
    rbeta(n, 1, 1)
  }
)

for (data_gen_name in names(data_generators)) {
  data_generator_fn <- data_generators[[data_gen_name]]


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
      robbins_subGaussian_cs(x = y, alpha = alpha, t_opt = t_opt, sigma2 = 1 / 4)
    }

  subexponential_fn <- function(y) {
    howard_subexponential_cs(
      x = y,
      alpha = alpha,
      v_opt = 1 / 4 * t_opt,
      running_intersection = FALSE
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
      mu = p,
      n_cores = parallel::detectCores()
    )

  print("Simulating fixed-time confidence interval miscoverage...")
  clt_miscoverage <-
    get_cumul_miscoverage_rate(
      data_generator_fn = data_generator_fn,
      conf_set_fn = clt_fn,
      times = start_time:n,
      num_repeats = num_repeats,
      mu = p,
      n_cores = parallel::detectCores()
    )

  print("Simulating sub-Gaussian miscoverage...")
  subGaussian_miscoverage <-
    get_cumul_miscoverage_rate(
      data_generator_fn = data_generator_fn,
      conf_set_fn = subGaussian_fn,
      times = start_time:n,
      num_repeats = num_repeats,
      mu = p,
      n_cores = parallel::detectCores()
    )


  print("Simulating sub-exponential miscoverage...")
  subexponential_miscoverage <-
    get_cumul_miscoverage_rate(
      data_generator_fn = data_generator_fn,
      conf_set_fn = subexponential_fn,
      times = start_time:n,
      num_repeats = num_repeats,
      mu = p,
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

  print("Simulating sub-Gaussian width...")
  subGaussian_width <- get_avg_width(
    data_generator_fn = data_generator_fn,
    conf_set_fn = subGaussian_fn,
    times = start_time:n,
    num_repeats = num_repeats,
    n_cores = parallel::detectCores()
  )

  print("Simulating sub-exponential width...")
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
    start_time,
    n,
    file = paste(r_data_dir, "widths_coverage",
      "_", data_gen_name,
      ".RData",
      sep = ""
    )
  )
}
