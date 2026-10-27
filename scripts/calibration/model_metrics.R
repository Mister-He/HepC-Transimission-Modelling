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

# Used by scenarios_2026.R; tests/integration/test_sim_and_pipeline.R checks the deaths against sim.cpp.
# Compartment table (position in the 480-state vector), built once for rate_components().
comp <- expand.grid(age = 0:5, h = 0:4, k = 1:4, s = 0:3)
comp$pos <- idx(comp$s, comp$k, comp$h, comp$age) - 1L
inj <- comp$s %in% c(0, 2)  # D and F inject

# Annual new infections (u, r -> a in D and F) and HCV deaths at a state, as sim.cpp computes them.
rate_components <- function(state, time, pm) {
  y <- as.numeric(state)[comp$pos]
  active <- tapply(y[inj], comp$age[inj], sum)
  infectious <- tapply(y[inj & comp$h %in% 1:2], comp$age[inj & comp$h %in% 1:2], sum)
  share <- ifelse(active > 0, infectious / active, 0)
  gam <- pm$q * as.vector(pm$C_contact %*% share)  # force of infection by age
  sus <- inj & comp$h %in% c(0, 4)
  inc <- gam[comp$age + 1] * y * sus
  # HCV deaths: liver deaths in DC/HCC plus the non-liver excess of viraemic people (a, c, t)
  # over the non-viraemic rate, mu * (smr_viraemic - smr_nonviraemic); cure removes both.
  treated <- comp$h %in% c(3, 4)
  liver <- ifelse(comp$k == 3, pm$mu_DC * ifelse(treated, pm$psi_DC, 1),
                  ifelse(comp$k == 4, pm$mu_HCC * ifelse(treated, pm$psi_HCC, 1), 0))
  excess <- ifelse(comp$h %in% 1:3, pm$mu[comp$age + 1] * (pm$smr_viraemic - pm$smr_nonviraemic), 0)
  deaths <- (liver + excess) * y
  c(incidence = sum(inc), mortality = sum(deaths),
    fx_incidence = sum(inc[comp$s == 2]), fx_mortality = sum(deaths[comp$s >= 2]))
}

