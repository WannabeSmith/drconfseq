library(drconfseq)

N <- 10000
alpha <- 0.1
t <- 1:N
mu_t <- (1 - sin(2 * log(exp(1) + 10*t)) / log(exp(1) + 0.001*t)) / 2
mu_t_tilde <- cumsum(mu_t) / t

x <- rbinom(n=N, size=1, prob=mu_t)

rho2 <- best_rho2_exact(t_opt = 500, alpha_opt=alpha)
asympcs <- lyapunov_asympcs(x=x, rho2=rho2, alpha=alpha)

lower_cs_lyapunov <- asympcs$l
upper_cs_lyapunov <- asympcs$u


r_data_dir <- "./paper_plots/simulations/time_varying/"
save(
  t,
  mu_t_tilde,
  lower_cs_lyapunov,
  upper_cs_lyapunov,
  N,
  file = paste(r_data_dir, "time_varying.RData", sep = "")
)
