# Extract the 2015 equilibrium snapshot used by calibration.
idx <- function(s, k, h, i)
  s * 4 * 5 * 6 + (k - 1) * 5 * 6 + h * 6 + i + 2L

j_pop_idx <- function(i) as.vector(sapply(1:4, function(k)
  sapply(0:4, function(h) idx(1, k, h, i - 1L))))
j_sero_idx <- function(i) as.vector(sapply(1:4, function(k)
  sapply(1:4, function(h) idx(1, k, h, i - 1L))))
j_viremic_idx <- function(i) as.vector(sapply(1:4, function(k)
  sapply(1:3, function(h) idx(1, k, h, i - 1L))))

J_summary_at <- function(out, row) {
  age     <- seq_along(targets$age_groups)
  N       <- vapply(age, function(i) sum(out[row, j_pop_idx(i)]), numeric(1))
  sero    <- vapply(age, function(i) sum(out[row, j_sero_idx(i)]), numeric(1))
  viremic <- vapply(age, function(i) sum(out[row, j_viremic_idx(i)]), numeric(1))
  data.frame(age_group = targets$age_groups,
  	     N_hat     = N,
             p_sero    = sero / N,
             p_viremic = viremic / N)
}

summary_row <- function(out, time = tail(out[, 1], 1))
  which.min(abs(out[, 1] - time))

J_summary_final <- function(out) J_summary_at(out, summary_row(out))

fit_metrics <- function(s) {
  p <- s$p_sero
  p_safe <- pmin(pmax(p, 1e-10), 1 - 1e-10)
  log_resid <- log(s$N_hat) - log(targets$population)
  data.frame(
    prevalence_rmse = sqrt(mean((p - targets$p_prev_binom)^2)),
    prevalence_max_abs_err = max(abs(p - targets$p_prev_binom)),
    population_mape = mean(abs(s$N_hat - targets$population) / targets$population),
    population_max_ape = max(abs(s$N_hat - targets$population) / targets$population),
    binomial_prevalence_deviance = 2 * sum(
      dbinom(targets$x_prev, targets$population, targets$p_prev_binom, log = TRUE) -
        dbinom(targets$x_prev, targets$population, p_safe, log = TRUE)),
    population_srss = sum((log_resid / targets$population_sigma_log)^2),
    nll_prev = -sum(dbinom(targets$x_prev, targets$population, p_safe, log = TRUE)),
    nll_pop = 0.5 * sum((log_resid / targets$population_sigma_log)^2)
  )
}
