# =============================================================================
# model_metrics.R — J-stratum summaries (model prison size, seroprevalence, viraemic
# prevalence) at the 2015 snapshot, the two count targets (first arrests by age, community
# PWID) and fit metrics against all calibration targets.
# =============================================================================

idx <- function(s, k, h, i)
  s * 4 * 5 * 6 + (k - 1) * 5 * 6 + h * 6 + i + 2L

j_pop_idx <- function(i) as.vector(sapply(1:4, function(k)
  sapply(0:4, function(h) idx(1, k, h, i - 1L))))
j_sero_idx <- function(i) as.vector(sapply(1:4, function(k)
  sapply(1:4, function(h) idx(1, k, h, i - 1L))))
stratum_idx <- function(s, i) as.vector(sapply(1:4, function(k)
  sapply(0:4, function(h) idx(s, k, h, i - 1L))))
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

# Viraemic (states a, c, t) share of all PWID (D, J, F, X; every stage and age) at a row.
# Reported only: mortality is set per state, so the share varies over time.
viraemic_share <- function(out, row = nrow(out)) {
  y <- out[row, 2:481]
  h <- outer((0:3) * 30, (1:3) * 6, `+`)
  sum(y[as.vector(outer(0:3 * 120, outer(as.vector(h), 1:6, `+`), `+`))]) / sum(y)
}

summary_row <- function(out, time = tail(out[, 1], 1))
  which.min(abs(out[, 1] - time))

J_summary_final <- function(out) J_summary_at(out, summary_row(out))

# Calibration summaries at a row (default: the final, 2015, row):
# - p_sero: J seroprevalence by age;
# - first_arrests: D -> J flow per year by age, c_first * lambda3_i * D_i (pm = the run's params);
# - community_pwid: D + F, all ages (injecting, not in prison). Entry is constant and treatment
#   starts in 2019, so the 2015 equilibrium value also stands for the 2017 survey year.
calib_summary <- function(out, pm, row = summary_row(out)) {
  age <- seq_along(targets$age_groups)
  D <- vapply(age, function(i) sum(out[row, stratum_idx(0, i)]), numeric(1))
  F <- vapply(age, function(i) sum(out[row, stratum_idx(2, i)]), numeric(1))
  list(p_sero = J_summary_at(out, row)$p_sero,
       first_arrests = pm$c_first * pm$lambda3 * D,
       community_pwid = sum(D) + sum(F))
}

fit_metrics <- function(cs) {
  p <- cs$p_sero; y <- targets$seroprevalence
  data.frame(
    prevalence_rmse = sqrt(mean((p - y)^2)),
    prevalence_max_abs_err = max(abs(p - y)),
    first_arrests_max_rel_err = max(abs(cs$first_arrests - targets$first_arrests) /
                                      pmax(targets$first_arrests, 1)),
    community_pwid = cs$community_pwid,
    nll_prev = nll_prev(p),
    nll_first_arrests = nll_first_arrests(cs$first_arrests),
    nll_community = nll_community(cs$community_pwid)
  )
}
