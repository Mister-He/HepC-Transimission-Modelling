# =============================================================================
# model_metrics.R — J (prison) target extraction and fit metrics
#
# The model tracks 5 HCV states: u (susceptible), a (acute), c (chronic),
# t (treatment), s (seropositive cleared/post-SVR). The prison screening
# target is anti-HCV serology, so the primary fitted prevalence is
# seroprevalence (a+c+t+s)/N; viremic prevalence (a+c+t)/N is a sensitivity.
# =============================================================================

# Compartment column for (stratum s, stage k, state h, age i):
# 3-strata model (D=0, J=1, X=2): C++ idx = s*120 + (k-1)*30 + h*6 + i;
# output column = idx + 2 (col 1 = time).
idx <- function(s, k, h, i) s * 4 * 5 * 6 + (k - 1) * 5 * 6 + h * 6 + i + 2L

j_pop_idx <- function(i) {
  as.vector(sapply(1:4, function(k) sapply(0:4, function(h) idx(1, k, h, i - 1L))))
}
j_sero_idx <- function(i) {
  as.vector(sapply(1:4, function(k) sapply(c(1, 2, 3, 4), function(h) idx(1, k, h, i - 1L))))
}
j_viremic_idx <- function(i) {
  as.vector(sapply(1:4, function(k) sapply(1:3, function(h) idx(1, k, h, i - 1L))))
}

J_summary_at <- function(out, row) {
  ages <- seq_along(cal_targets$age_groups)
  N_hat <- sapply(ages, function(i) sum(out[row, j_pop_idx(i)]))
  I_sero <- sapply(ages, function(i) sum(out[row, j_sero_idx(i)]))
  I_vir  <- sapply(ages, function(i) sum(out[row, j_viremic_idx(i)]))
  data.frame(
    age_group = cal_targets$age_groups,
    i         = ages,
    N_hat     = N_hat,
    p_sero    = I_sero / N_hat,
    p_viremic = I_vir / N_hat
  )
}

J_summary_final <- function(out, target_time = NULL) {
  if (is.null(target_time) && exists("TARGET_TIME", inherits = TRUE)) {
    target_time <- get("TARGET_TIME", inherits = TRUE)
  }
  if (!is.null(target_time)) {
    times <- out[, 1]
    row <- which.min(abs(times - target_time))
  } else {
    row <- nrow(out)
  }
  J_summary_at(out, row)
}

fit_metrics <- function(s, target_mode = TARGET_MODE) {
  N_obs <- cal_targets$prison_total
  prev_obs <- cal_targets$prev_binom
  p_hat <- if (target_mode == "sero") s$p_sero else s$p_viremic
  prev_err <- p_hat - prev_obs
  pop_ape <- abs(s$N_hat - N_obs) / N_obs

  prevalence_rmse <- sqrt(mean(prev_err^2))
  prevalence_max_abs_err <- max(abs(prev_err))
  population_mape <- mean(pop_ape)
  population_max_ape <- max(pop_ape)

  p_safe <- pmin(pmax(p_hat, 1e-10), 1 - 1e-10)
  dev_prev <- pmax(
    2 * sum(
      dbinom(cal_targets$x_prev, size = cal_targets$n_prev,
             prob = prev_obs, log = TRUE) -
        dbinom(cal_targets$x_prev, size = cal_targets$n_prev,
               prob = p_safe, log = TRUE)
    ), 0
  )
  srss_pop <- sum(((log(s$N_hat) - log(N_obs)) / 0.10)^2)

  data.frame(
    prevalence_rmse        = prevalence_rmse,
    prevalence_max_abs_err = prevalence_max_abs_err,
    population_mape        = population_mape,
    population_max_ape     = population_max_ape,
    binomial_prevalence_deviance = dev_prev,
    population_srss              = srss_pop,
    nll_prev = -sum(dbinom(cal_targets$x_prev, size = cal_targets$n_prev,
                           prob = p_safe, log = TRUE)),
    nll_pop  = 0.5 * srss_pop
  )
}
