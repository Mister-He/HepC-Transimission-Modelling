# =============================================================================
# model_metrics.R — J-stratum summaries (model prison size, seroprevalence, viraemic
# prevalence) at the 2015 snapshot and fit metrics against the seroprevalence targets.
# =============================================================================

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

fit_metrics <- function(s) {
  p <- s$p_sero; y <- targets$seroprevalence
  data.frame(
    prevalence_rmse = sqrt(mean((p - y)^2)),
    prevalence_max_abs_err = max(abs(p - y)),
    nll_prev = nll_prev(p)
  )
}
