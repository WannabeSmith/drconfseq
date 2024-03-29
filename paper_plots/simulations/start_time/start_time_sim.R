library(drconfseq)
library(parallel)
start_times <- c(5, 10, 25, 50)
n <- 10000
num_repeats <- 10000
n_cores <- parallel::detectCores()
df <- 4
alpha <- 0.1
stopifnot(all(n >= start_times))

miscoverage_rates <- list()

data_generator_fn = function() {
  rt(n, df) + rexp(n, 1) - 1
}

for (start_time in start_times)
{
  conf_set_fn <- function(y)
  {
    asymptotic_confseq(
      x = y,
      t_opt = 5 * start_time,
      alpha = alpha,
      return_all_times = TRUE
    )
  }

  print(paste("Computing miscoverage for start time:", start_time))
  miscoverage_rate <-
    get_cumul_miscoverage_rate(
      data_generator_fn = data_generator_fn,
      conf_set_fn = conf_set_fn,
      times = start_time:n,
      num_repeats = num_repeats,
      n_cores = n_cores
    )

  miscoverage_rates[[as.character(start_time)]] <- miscoverage_rate
}

save(miscoverage_rates, alpha,
     file = "./start_time.RData")
