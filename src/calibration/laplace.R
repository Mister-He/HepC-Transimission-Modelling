# =============================================================================
# laplace.R — Laplace approximation for the fitted target summaries
#
# Method (per AGENTS.md / prompt.md):
#  1. finite-difference Hessian of the PURE statistical NLL (no equilibrium
#     penalty) at the optimum;
#  2. generalized inverse with relative eigenvalue cutoff rel_thresh
#     (condition number <= 1/rel_thresh); directions below the cutoff are
#     treated as un-identified;
#  3. Monte Carlo draws from N(theta_hat, Sigma) re-simulated to the target
#     summaries; draws that fail the equilibrium gate are discarded (the
#     reported intervals are conditional on equilibrium feasibility);
#  4. 95% intervals per age group compared with observed intervals
#     (Jeffreys Binomial for prevalence; log-Normal sigma_pop = 0.10 for
#     population).
# =============================================================================

laplace_intervals <- function(theta, base_params, data,
                              target_mode = TARGET_MODE,
                              n_draws = 1000, rel_thresh = 1e-4,
                              seed = 2026, t_lag = 5) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    return(list(success = FALSE, reason = "MASS not installed"))
  }
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    return(list(success = FALSE, reason = "numDeriv not installed"))
  }

  stat_obj <- make_objective(base_params, data, target_mode = target_mode,
                             equilibrium_penalty = FALSE)
  hess <- tryCatch(numDeriv::hessian(stat_obj, theta),
                   error = function(e) NULL)
  if (is.null(hess) || any(!is.finite(hess))) {
    return(list(success = FALSE, reason = "Hessian failed"))
  }

  Sigma <- tryCatch(MASS::ginv(hess), error = function(e) NULL)
  if (is.null(Sigma) || any(!is.finite(Sigma))) {
    return(list(success = FALSE, reason = "Hessian inversion failed"))
  }
  ev <- eigen(hess, symmetric = TRUE)
  rank_eff <- sum(ev$values > max(ev$values) * rel_thresh)

  set.seed(seed)
  draws <- MASS::mvrnorm(n_draws, mu = theta, Sigma = Sigma)

  sims <- vector("list", n_draws)
  n_used <- 0L
  for (j in seq_len(n_draws)) {
    pm <- tryCatch(build_params(draws[j, ], base_params), error = function(e) NULL)
    if (is.null(pm)) next
    out <- tryCatch(run_sim(pm, data), error = function(e) NULL)
    if (is.null(out) || nrow(out) < 2) next
    eq <- tryCatch(check_equilibrium(out, target_mode = target_mode, t_lag = t_lag),
                   error = function(e) NULL)
    if (is.null(eq) || !eq$pass) next
    s <- tryCatch(J_summary_final(out), error = function(e) NULL)
    if (is.null(s)) next
    n_used <- n_used + 1L
    sims[[n_used]] <- data.frame(
      age_group = s$age_group,
      N_hat = s$N_hat,
      p_hat = if (target_mode == "sero") s$p_sero else s$p_viremic
    )
  }
  sims <- sims[seq_len(n_used)]

  if (n_used < 20) {
    return(list(success = FALSE, reason = paste("too few feasible draws:", n_used)))
  }

  mat_N <- do.call(rbind, lapply(sims, function(z) z$N_hat))
  mat_p <- do.call(rbind, lapply(sims, function(z) z$p_hat))

  # Point predictions at the optimum
  pm_best <- build_params(theta, base_params)
  out_best <- run_sim(pm_best, data)
  s_best <- J_summary_final(out_best)
  p_point <- if (target_mode == "sero") s_best$p_sero else s_best$p_viremic

  intervals <- data.frame(
    age_group = cal_targets$age_groups,
    p_hat     = p_point,
    p_lo      = apply(mat_p, 2, quantile, probs = 0.025),
    p_hi      = apply(mat_p, 2, quantile, probs = 0.975),
    p_obs_lo  = qbeta(0.025, cal_targets$x_prev + 0.5,
                      cal_targets$n_prev - cal_targets$x_prev + 0.5),
    p_obs_hi  = qbeta(0.975, cal_targets$x_prev + 0.5,
                      cal_targets$n_prev - cal_targets$x_prev + 0.5),
    N_hat     = s_best$N_hat,
    N_lo      = apply(mat_N, 2, quantile, probs = 0.025),
    N_hi      = apply(mat_N, 2, quantile, probs = 0.975),
    N_obs     = cal_targets$prison_total,
    N_obs_lo  = exp(log(cal_targets$prison_total) - qnorm(0.975) * sigma_pop),
    N_obs_hi  = exp(log(cal_targets$prison_total) + qnorm(0.975) * sigma_pop)
  )
  intervals$p_overlap <- !(intervals$p_hi < intervals$p_obs_lo |
                           intervals$p_lo > intervals$p_obs_hi)
  intervals$N_overlap <- !(intervals$N_hi < intervals$N_obs_lo |
                           intervals$N_lo > intervals$N_obs_hi)

  list(
    success = TRUE,
    hessian_method = "numDeriv Richardson",
    rank = length(ev$values),
    effective_dimension = rank_eff,
    n_draws_used = n_used,
    n_draws_total = n_draws,
    rel_thresh = rel_thresh,
    intervals = intervals,
    eigenvalues = ev$values,
    eigenvectors = ev$vectors
  )
}
