# =============================================================================
# equilibrium.R — stability checks for the J target summaries
#
# "Equilibrium" = the age-specific summaries are stable over the final 10
# years of the simulation horizon (T vs T-10).
# Population stability (overall & per age grp): |log(N(T) / N(T-10)| \leq 0.01
# Default convergence: max(|p_a(T) - p_a(T-10)|) \leq 0.005
# =============================================================================

# Run settings for the burn-in to the 2015 snapshot (model year = horizon): start at 0, one stored
# row per model year (callers read only yearly rows: the snapshot and T-10). Set dt before calling.
equilibrium_data <- function(data, t_end = targets$equilibrium_horizon_years) {
  data$t_start <- 0; data$t_end <- t_end; data$save_every <- round(1 / data$dt); data
}

EQ_COMP <- state_cols(0:3); EQ_HCV <- state_cols(0:3, 1:4)
EQ_DC <- state_cols(0:3, k = 3); EQ_HCC <- state_cols(0:3, k = 4)

check_equilibrium <- function(out, t_lag = 10,
                              crit_log_pop = 0.01,
                              crit_prev = 0.005,
                              crit_state = 0.01,
                              crit_comp = 0.02) {
  times <- out[, 1]
  row_T  <- nrow(out)
  T      <- times[row_T]
  row_T5 <- which.min(abs(times - (T - t_lag)))

  s_T  <- J_summary_at(out, row_T)
  s_T5 <- J_summary_at(out, row_T5)

  log_ratio   <- abs(log(s_T$N_hat / s_T5$N_hat))
  prev_change <- abs(s_T$p_sero - s_T5$p_sero)

  # All-state convergence: total HCV, DC, HCC, and every compartment (state_cols from model_metrics.R).
  comp_T  <- out[row_T, EQ_COMP]
  comp_T5 <- out[row_T5, EQ_COMP]
  max_comp_log_ratio <- max(abs(log(pmax(comp_T, 1e-6) / pmax(comp_T5, 1e-6))))

  sum_log_ratio <- function(cols) abs(log(sum(out[row_T, cols]) / sum(out[row_T5, cols])))
  state_log_ratio <- max(sum_log_ratio(EQ_HCV), sum_log_ratio(EQ_DC), sum_log_ratio(EQ_HCC))

  violation <- sum(pmax(log_ratio / crit_log_pop - 1, 0)) +
    pmax(max(prev_change) / crit_prev - 1, 0) +
    pmax(state_log_ratio / crit_state - 1, 0) +
    pmax(max_comp_log_ratio / crit_comp - 1, 0)
  pass <- isTRUE(max(log_ratio) <= crit_log_pop &&
    max(prev_change) <= crit_prev &&
    state_log_ratio <= crit_state &&
    max_comp_log_ratio <= crit_comp)

  list(
    pass              = pass,
    T                 = T,
    T_minus_lag       = times[row_T5],
    max_log_pop_ratio = max(log_ratio),
    max_prev_change   = max(prev_change),
    violation         = violation,
    state_log_ratio   = state_log_ratio,
    max_comp_log_ratio = max_comp_log_ratio,
    by_age = data.frame(
      age_group   = s_T$age_group,
      N_hat_T     = s_T$N_hat,
      N_hat_T5    = s_T5$N_hat,
      log_ratio   = log_ratio,
      p_hat_T     = s_T$p_sero,
      p_hat_T5    = s_T5$p_sero,
      prev_change = prev_change
    )
  )
}
