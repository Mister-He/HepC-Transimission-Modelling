# =============================================================================
# mcmc.R — exact log posterior and adaptive random-walk Metropolis sampler used for
# the reference MCMC posterior, plus posterior summaries.
# =============================================================================

# Exact target shared with NPE: independent Normal prior on theta (log scales),
# Binomial prevalence + log-normal population likelihood, and the 2015
# equilibrium check as a hard support restriction (failure -> -Inf, no penalty).
# Returns list(lp, ode) so that ODE calls can be counted separately.
log_posterior_exact <- function(theta, base_params, data, prior_mean, prior_sd) {
  lprior <- sum(dnorm(theta, prior_mean, prior_sd, log = TRUE))
  if (!within_bounds(theta)) return(list(lp = -Inf, ode = 0L))
  pm <- tryCatch(build_params(theta, base_params), error = function(e) NULL)
  if (is.null(pm)) return(list(lp = -Inf, ode = 0L))
  out <- tryCatch(run_sim(pm, data), error = function(e) NULL)
  if (is.null(out) || nrow(out) < 2) return(list(lp = -Inf, ode = 1L))
  s <- J_summary_final(out)
  if (!all(is.finite(s$N_hat)) || !all(is.finite(s$p_sero)) || any(s$N_hat <= 0))
    return(list(lp = -Inf, ode = 1L))
  eq <- tryCatch(check_equilibrium(out), error = function(e) NULL)
  if (is.null(eq) || !eq$pass) return(list(lp = -Inf, ode = 1L))
  list(lp = lprior - nll_prev(s$p_sero) - nll_pop(s$N_hat), ode = 1L)
}

# Adaptive random-walk Metropolis with a fixed evaluation budget (Haario et al. 2001 style).
# - init: draw from the prior until the target is finite; every try uses budget.
# - warmup: n_warmup proposals; proposal covariance = Sigma_init for the first
#   adapt_start proposals, then a mixture: with prob safety_p Sigma_init, else
#   (2.38^2/d) * (empirical covariance of warmup states + reg), refreshed every
#   adapt_every proposals. Adaptation is frozen before the first retained proposal.
# - sampling: retained proposals until evaluations reach max_evals. Every retained row
#   stores cumulative evaluations, ODE calls and elapsed time.
arwm_budget <- function(logpost, draw_prior, Sigma_init, reg, n_warmup = 4999L,
                        max_evals = 50000L, adapt_start = 200L,
                        adapt_every = 50L, safety_p = 0.05, seed = 1L,
                        t_offset = 0, max_init_tries = 100L) {
  set.seed(seed)
  t0 <- proc.time()[["elapsed"]] - t_offset
  elapsed <- function() proc.time()[["elapsed"]] - t0
  d <- ncol(Sigma_init); c_d <- 2.38^2 / d
  L_init <- chol(Sigma_init)
  evals <- 0L; ode <- 0L
  ev <- function(th) { r <- logpost(th); evals <<- evals + 1L; ode <<- ode + r$ode; r$lp }

  # Initialisation from the prior (budgeted).
  lp <- -Inf; tries <- 0L
  while (!is.finite(lp) && tries < max_init_tries && evals < max_evals) {
    theta <- draw_prior(); lp <- ev(theta); tries <- tries + 1L
  }
  if (!is.finite(lp)) return(list(status = "init_failed", init_tries = tries, evals = evals, ode = ode))
  init_evals <- evals

  warm <- matrix(NA_real_, n_warmup + 1L, d); warm[1, ] <- theta
  L_adapt <- NULL; acc_w <- 0L; status <- "ok"
  for (t in seq_len(n_warmup)) {
    if (evals >= max_evals) { status <- "budget_exhausted_in_warmup"; break }
    if (t > adapt_start && ((t - adapt_start - 1L) %% adapt_every == 0L || is.null(L_adapt))) {
      S <- stats::cov(warm[seq_len(t), , drop = FALSE]) + reg
      L_adapt <- tryCatch(chol(c_d * S), error = function(e) L_init)
    }
    L <- if (t <= adapt_start || runif(1) < safety_p) L_init else L_adapt
    prop <- theta + as.numeric(crossprod(L, rnorm(d)))
    lp_prop <- ev(prop)
    if (is.finite(lp_prop) && log(runif(1)) < lp_prop - lp) { theta <- prop; lp <- lp_prop; acc_w <- acc_w + 1L }
    warm[t + 1L, ] <- theta
  }
  warmup_evals <- evals; warmup_sec <- elapsed()
  if (status != "ok") return(list(status = status, evals = evals, ode = ode, warmup = warm))
  if (is.null(L_adapt)) L_adapt <- L_init  # frozen from here on

  cap <- max(1000L, max_evals - evals) * 4L
  keep <- matrix(NA_real_, cap, d); info <- matrix(NA_real_, cap, 5)
  colnames(info) <- c("lp", "accepted", "evals", "ode", "elapsed_sec")
  n <- 0L
  while (evals < max_evals) {
    L <- if (runif(1) < safety_p) L_init else L_adapt
    prop <- theta + as.numeric(crossprod(L, rnorm(d)))
    lp_prop <- ev(prop); acc <- 0
    if (is.finite(lp_prop) && log(runif(1)) < lp_prop - lp) { theta <- prop; lp <- lp_prop; acc <- 1 }
    n <- n + 1L
    if (n > nrow(keep)) {
      keep <- rbind(keep, matrix(NA_real_, cap, d)); info <- rbind(info, matrix(NA_real_, cap, 5))
    }
    keep[n, ] <- theta; info[n, ] <- c(lp, acc, evals, ode, elapsed())
  }
  list(status = status, samples = keep[seq_len(n), , drop = FALSE],
       info = info[seq_len(n), , drop = FALSE], warmup = warm,
       init_tries = tries, init_evals = init_evals, warmup_evals = warmup_evals,
       warmup_sec = warmup_sec, warmup_acceptance = acc_w / n_warmup,
       proposal_chol_adapt = L_adapt, total_evals = evals, total_ode = ode,
       total_sec = elapsed())
}

posterior_summaries <- function(samples, theta_names = paste0("theta", 1:12)) {
  d <- ncol(samples)
  data.frame(
    parameter = theta_names,
    median = apply(samples, 2, median),
    mean   = apply(samples, 2, mean),
    sd     = apply(samples, 2, sd),
    q025   = apply(samples, 2, quantile, probs = 0.025),
    q975   = apply(samples, 2, quantile, probs = 0.975)
  )
}
