# =============================================================================
# targets.R — calibration target data and Binomial count reconstruction
# =============================================================================

cal_targets <- list(
  age_groups    = c("<20", "20-29", "30-39", "40-49", "50-59", "60+"),
  prev_supplied = c(0.1118421, 0.1731044, 0.2684954, 0.4301165,
                    0.4821029, 0.3544304),
  prison_total  = c(99, 1244, 1467, 1841, 1628, 409)
)

cal_targets$n_prev <- cal_targets$prison_total
# Integer HCV-positive counts (raw integers unavailable; reconstructed and
# recorded per AGENTS.md). Rounding differences are tabulated.
cal_targets$x_prev <- round(cal_targets$prev_supplied * cal_targets$n_prev)
cal_targets$prev_binom    <- cal_targets$x_prev / cal_targets$n_prev
cal_targets$rounding_diff <- cal_targets$prev_binom - cal_targets$prev_supplied

targets_table <- data.frame(
  age_group            = cal_targets$age_groups,
  n_prev               = cal_targets$n_prev,
  x_prev               = cal_targets$x_prev,
  target_prev_supplied = cal_targets$prev_supplied,
  target_prev_binom    = cal_targets$prev_binom,
  rounding_difference  = cal_targets$rounding_diff
)
