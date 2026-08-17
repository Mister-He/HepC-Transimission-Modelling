# =============================================================================
# equilibrium.R — stability checks for the J target summaries
#
# "Equilibrium" = the age-specific summaries are stable over the final 5
# years of the simulation horizon (T vs T-5).
# =============================================================================

check_equilibrium <- function(out, t_lag = 5,
                              target_mode = if (exists("TARGET_MODE")) TARGET_MODE else "sero",
                              crit_log_pop = 0.01,
                              crit_prev = 0.005,
                              crit_total_pop = 0.01,
                              crit_state = 0.01,
                              crit_comp = 0.02) {
  times <- out[, 1]
  row_T  <- nrow(out)
  T      <- times[row_T]
  row_T5 <- which.min(abs(times - (T - t_lag)))

  s_T  <- J_summary_at(out, row_T)
  s_T5 <- J_summary_at(out, row_T5)

  log_ratio   <- abs(log(s_T$N_hat / s_T5$N_hat))
  prev_col <- if (target_mode == "sero") "p_sero" else "p_viremic"
  prev_change <- abs(s_T[[prev_col]] - s_T5[[prev_col]])

  tot_T  <- sum(out[row_T, -1])
  tot_T5 <- sum(out[row_T5, -1])
  total_log_ratio <- abs(log(tot_T / tot_T5))

  # All-state convergence: total HCV, DC, HCC, and every compartment.
  idx <- function(s, k, h, i) s * 4 * 5 * 6 + (k - 1) * 5 * 6 + h * 6 + i + 2L
  comp_idx <- as.vector(sapply(0:2, function(s) sapply(1:4, function(k)
    sapply(0:4, function(h) sapply(0:5, function(i) idx(s, k, h, i))))))
  comp_T  <- out[row_T, comp_idx]
  comp_T5 <- out[row_T5, comp_idx]
  max_comp_log_ratio <- max(abs(log(pmax(comp_T, 1e-6) / pmax(comp_T5, 1e-6))))

  hcv_idx <- as.vector(sapply(0:2, function(s) sapply(1:4, function(k)
    sapply(1:3, function(h) sapply(0:5, function(i) idx(s, k, h, i))))))
  dc_idx <- as.vector(sapply(0:2, function(s) sapply(0:4, function(h)
    sapply(0:5, function(i) idx(s, 3, h, i)))))
  hcc_idx <- as.vector(sapply(0:2, function(s) sapply(0:4, function(h)
    sapply(0:5, function(i) idx(s, 4, h, i)))))
  state_log_ratio <- max(
    abs(log(sum(out[row_T, hcv_idx]) / sum(out[row_T5, hcv_idx]))),
    abs(log(sum(out[row_T, dc_idx]) / sum(out[row_T5, dc_idx]))),
    abs(log(sum(out[row_T, hcc_idx]) / sum(out[row_T5, hcc_idx])))
  )

  # The gate applies to the J target summaries; the total-PWID ratio
  # (including the accumulating X pool) is reported as an inspection metric.
  pass <- max(log_ratio) <= crit_log_pop &&
    max(prev_change) <= crit_prev &&
    state_log_ratio <= crit_state &&
    max_comp_log_ratio <= crit_comp

  list(
    pass              = pass,
    T                 = T,
    T_minus_lag       = times[row_T5],
    max_log_pop_ratio = max(log_ratio),
    max_prev_change   = max(prev_change),
    total_log_ratio   = total_log_ratio,
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
