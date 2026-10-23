# =============================================================================
# mcmc.R — exact log posterior and adaptive random-walk Metropolis sampler used for
# the reference MCMC posterior, plus posterior summaries.
# =============================================================================

# Exact target shared with NPE and Nelder-Mead: the truncated Normal prior of likelihood.R
# (-nlp, up to a constant), the Beta seroprevalence likelihood, and the 2015 equilibrium check
# as a hard support restriction (failure -> -Inf, no penalty).
# Returns list(lp, ode) so that ODE calls can be counted separately.
log_posterior_exact <- function(theta, base_params, data) {
  lprior <- -nlp(theta)
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
  list(lp = lprior - nll_prev(s$p_sero), ode = 1L)
}

# Adaptive random-walk Metropolis with a fixed evaluation budget (Haario et al. 2001 style).
# - init: draw from the prior until the target is finite; every try uses budget.
# - warmup: n_warmup proposals; proposal covariance = Sigma_init for the first
#   adapt_start proposals, then a mixture: with prob safety_p Sigma_init, else
#   (2.38^2/d) * (empirical covariance of warmup states + reg), refreshed every
#   adapt_every proposals. Adaptation is frozen before the first retained proposal.
# - sampling: retained proposals until evaluations reach max_evals. Every retained row
#   stores cumulative evaluations, ODE calls, elapsed time and proposal jump size.
# Diagnostics returned: adaptation trace (relative Frobenius change of the adapted
# covariance between updates and its condition number), warmup jump sizes, CPU time.
arwm_budget <- function(logpost, draw_prior, Sigma_init, reg, n_warmup = 4999L,
                        max_evals = 50000L, adapt_start = 200L,
                        adapt_every = 50L, safety_p = 0.05, seed = 1L,
                        t_offset = 0, max_init_tries = 100L) {
  set.seed(seed)
  cpu0 <- sum(proc.time()[c("user.self", "sys.self")])
  t0 <- proc.time()[["elapsed"]] - t_offset
  elapsed <- function() proc.time()[["elapsed"]] - t0
  d <- ncol(Sigma_init); c_d <- 2.38^2 / d
  L_init <- chol(Sigma_init)
  evals <- 0L; ode <- 0L
  ev <- function(th) { r <- logpost(th); evals <<- evals + 1L; ode <<- ode + r$ode; r$lp }

  # Initialisation: prior draws until the target is finite (all budgeted).
  lp <- -Inf; tries <- 0L; start_source <- "prior"
  while (!is.finite(lp) && tries < max_init_tries && evals < max_evals) {
    theta <- draw_prior(); lp <- ev(theta); tries <- tries + 1L
  }
  if (!is.finite(lp)) return(list(status = "init_failed", init_tries = tries, evals = evals, ode = ode))
  init_evals <- evals

  warm <- matrix(NA_real_, n_warmup + 1L, d); warm[1, ] <- theta
  warm_jump <- matrix(NA_real_, n_warmup, 2, dimnames = list(NULL, c("jump", "accepted")))
  L_adapt <- NULL; acc_w <- 0L; status <- "ok"; S_prev <- NULL; S_adapt <- NULL
  adapt_trace <- list()
  for (t in seq_len(n_warmup)) {
    if (evals >= max_evals) { status <- "budget_exhausted_in_warmup"; break }
    if (t > adapt_start && ((t - adapt_start - 1L) %% adapt_every == 0L || is.null(L_adapt))) {
      # Recent half of warmup only: the walk-in from the prior start would inflate the proposal.
      S <- stats::cov(warm[max(1L, t %/% 2L):t, , drop = FALSE]) + reg
      L_adapt <- tryCatch(chol(c_d * S), error = function(e) L_init)
      ev_S <- eigen(c_d * S, symmetric = TRUE, only.values = TRUE)$values
      drift <- if (is.null(S_prev)) NA_real_ else norm(S - S_prev, "F") / max(norm(S_prev, "F"), 1e-12)
      adapt_trace[[length(adapt_trace) + 1L]] <- c(iteration = t, drift = drift,
                                                   condition_number = max(ev_S) / min(ev_S))
      S_prev <- S; S_adapt <- c_d * S
    }
    L <- if (t <= adapt_start || runif(1) < safety_p) L_init else L_adapt
    step <- as.numeric(crossprod(L, rnorm(d))); prop <- theta + step
    lp_prop <- ev(prop); acc <- 0
    if (is.finite(lp_prop) && log(runif(1)) < lp_prop - lp) { theta <- prop; lp <- lp_prop; acc_w <- acc_w + 1L; acc <- 1 }
    warm[t + 1L, ] <- theta; warm_jump[t, ] <- c(sqrt(sum(step^2)), acc)
  }
  warmup_evals <- evals; warmup_sec <- elapsed()
  if (status != "ok") return(list(status = status, evals = evals, ode = ode, warmup = warm))
  if (is.null(L_adapt)) L_adapt <- L_init  # frozen from here on

  cap <- max_evals - evals  # each retained iteration costs exactly one evaluation
  keep <- matrix(NA_real_, cap, d); info <- matrix(NA_real_, cap, 6)
  colnames(info) <- c("lp", "accepted", "evals", "ode", "elapsed_sec", "jump")
  n <- 0L
  while (evals < max_evals) {
    L <- if (runif(1) < safety_p) L_init else L_adapt
    step <- as.numeric(crossprod(L, rnorm(d))); prop <- theta + step
    lp_prop <- ev(prop); acc <- 0
    if (is.finite(lp_prop) && log(runif(1)) < lp_prop - lp) { theta <- prop; lp <- lp_prop; acc <- 1 }
    n <- n + 1L
    keep[n, ] <- theta; info[n, ] <- c(lp, acc, evals, ode, elapsed(), sqrt(sum(step^2)))
  }
  list(status = status, samples = keep[seq_len(n), , drop = FALSE],
       info = info[seq_len(n), , drop = FALSE], warmup = warm,
       init_tries = tries, init_evals = init_evals, warmup_evals = warmup_evals,
       warmup_sec = warmup_sec, warmup_acceptance = acc_w / n_warmup,
       proposal_chol_adapt = L_adapt, proposal_cov_adapt = S_adapt,
       adapt_trace = if (length(adapt_trace)) do.call(rbind, adapt_trace) else NULL,
       warmup_jump = warm_jump, start_source = start_source,
       total_evals = evals, total_ode = ode, total_sec = elapsed(),
       cpu_sec = sum(proc.time()[c("user.self", "sys.self")]) - cpu0)
}

posterior_summaries <- function(samples, theta_names = paste0("theta", seq_len(ncol(samples)))) {
  data.frame(
    parameter = theta_names,
    median = apply(samples, 2, median),
    mean   = apply(samples, 2, mean),
    sd     = apply(samples, 2, sd),
    q025   = apply(samples, 2, quantile, probs = 0.025),
    q975   = apply(samples, 2, quantile, probs = 0.975)
  )
}
