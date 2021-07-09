library(sequential.causal)
library(parallel)
n <- 100000
start_time <- 500
t_opt <- 1000
alpha <- 0.1
p_1 <- 0.8
p_2 <- 0.4
prop <- 0.5

data_generator_fn <- function() {
  treatment <- rbinom(n, 1, prop)

  y <- rbinom(n, 1, ifelse(treatment, p_1, p_2))
  t <- 1:n
  t_1 <- cumsum(treatment[treatment == 1])
  t_2 <- cumsum(1 - treatment[treatment == 0])

  pseudo_outcomes <- (treatment / prop - (1 - treatment) / prop) * y
  return(pseudo_outcomes)
}

acs_fn <- function(y) {
  asymptotic_confseq(x = y,
                     t_opt = t_opt,
                     return_all_times = TRUE)
}
clt_fn <-
  function(y) {
    naive_confidence_intervals(x = y, return_all_times = TRUE)
  }

print("Simulating time-uniform confidence sequence miscoverage...")
acs_miscoverage <-
  get_cumul_miscoverage_rate(
    data_generator_fn = data_generator_fn,
    conf_set_fn = acs_fn,
    times = start_time:n,
    num_repeats = 1000,
    mu = p_1 - p_2,
    n_cores = parallel::detectCores()
  )
print("Simulating fixed-time confidence interval miscoverage...")
clt_miscoverage <-
  get_cumul_miscoverage_rate(
    data_generator_fn = data_generator_fn,
    conf_set_fn = clt_fn,
    times = start_time:n,
    num_repeats = 1000,
    mu = p_1 - p_2,
    n_cores = parallel::detectCores()
  )

y <- data_generator_fn()
acs <- acs_fn(y)
clt <- clt_fn(y)

r_data_dir <- "./"

save(
  acs_miscoverage,
  clt_miscoverage,
  acs,
  clt,
  start_time,
  p_1,
  p_2,
  n,
  file = paste(r_data_dir, 'cs_vs_ci.RData', sep = '')
)
