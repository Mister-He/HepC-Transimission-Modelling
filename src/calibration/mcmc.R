# =============================================================================
# mcmc.R — Bayesian posterior sampling for the 12 calibration parameters
#
# Components:
#   make_priors()      documented weakly informative priors (see below)
#   log_prior()        log-density of the priors
#   log_posterior()    log-likelihood + log-prior - equilibrium penalty
#   am_mcmc()          adaptive Metropolis sampler (Haario et al. 2001)
#   posterior_summaries()  medians/quantiles for parameters
#   predictive_quantiles() equal-tailed 95% CrI for target summaries
#
# Priors (documented in docs/calibration/bayes_methodology.md):
#   log(contact_scale) ~ Normal(0, 2^2)         theta[1:6]
#   log(beta_scale)    ~ Student-t(3, 0, 2)     theta[7:12]
# Rationale: the base contact matrix is a model-specific calibrated guess
# (weak prior); base beta is anchored to CNB official new-drug-abuser data
# but the per-age-group values are placeholders, so a heavy-tailed
# Student-t prior centres the inflow scale near the official anchor while
# tolerating the large deviations the data require (e.g. 60+ inflow).
# =============================================================================

make_priors <- function(contact_mean = 0, contact_sd = 2,
                        beta_mean = 0, beta_sd = 2, beta_df = 3) {
  list(
    contact_mean = contact_mean, contact_sd = contact_sd,
    beta_mean = beta_mean, beta_sd = beta_sd, beta_df = beta_df
  )
}

log_prior <- function(theta, priors) {
  lp <- sum(dnorm(theta[1:6], priors$contact_mean, priors$contact_sd, log = TRUE))
  z <- (theta[7:12] - priors$beta_mean) / priors$beta_sd
  lp <- lp + sum(dt(z, df = priors$beta_df, log = TRUE)) -
    6 * log(priors$beta_sd)
  lp
}

log_posterior <- function(theta, base_params, data,
                          priors, target_mode = TARGET_MODE,
                          equilibrium_penalty = TRUE) {
  pm <- tryCatch(build_params(theta, base_params), error = function(e) NULL)
  if (is.null(pm)) return(-Inf)
  out <- tryCatch(run_sim(pm, data), error = function(e) NULL)
  if (is.null(out) || nrow(out) < 2) return(-Inf)

  s <- J_summary_final(out)
  p_hat <- if (target_mode == "sero") s$p_sero else s$p_viremic
  if (!all(is.finite(s$N_hat)) || !all(is.finite(p_hat)) ||
      any(s$N_hat <= 0)) return(-Inf)

  lp <- -(nll_prev(p_hat) + nll_pop(s$N_hat)) + log_prior(theta, priors)
  if (equilibrium_penalty) {
    eq <- tryCatch(check_equilibrium(out, target_mode = target_mode),
                   error = function(e) NULL)
    if (!is.null(eq) && !eq$pass) {
      lp <- lp - EQ_PENALTY_FACTOR * (
        pmax(eq$max_log_pop_ratio / 0.01, 0) +
        pmax(eq$max_prev_change   / 0.005, 0) +
        pmax(eq$total_log_ratio   / 0.01, 0))
    }
  }
  lp
}

# Adaptive Metropolis (Haario, Saksman & Tamminen 2001): Gaussian proposal
# centred at the current state with covariance c_d * Sigma_t + eps*I, where
# Sigma_t is the empirical covariance of the chain so far and
# c_d = (2.38^2)/d. The proposal covariance is initialised from Sigma0
# (Laplace covariance or NPE posterior covariance). Supports continuation:
# pass `init_state` (a previous block's $state) to warm-restart from the
# last state instead of discarding earlier samples.
am_mcmc <- function(theta0, logpost, n_iter, burnin = 0, Sigma0,
                    seed = 1, adapt_from = 200, thin = 1,
                    init_state = NULL) {
  set.seed(seed)
  d <- length(theta0)
  c_d <- (2.38^2) / d
  eps <- 1e-8 * diag(d)

  if (is.null(init_state)) {
    theta <- theta0
    lp <- logpost(theta)
    mean_t <- theta
    M2 <- matrix(0, d, d)
    n <- 1L
    n_accept <- 0L
    total0 <- 0L
    burnin_eff <- burnin
  } else {
    theta <- init_state$theta
    lp <- init_state$logpost
    mean_t <- init_state$mean
    M2 <- init_state$M2
    n <- init_state$n
    n_accept <- init_state$n_accept
    total0 <- init_state$n_iter
    burnin_eff <- 0L
  }
  if (!is.finite(lp)) stop("log-posterior not finite at starting point")

  n_keep <- floor((n_iter - burnin_eff) / thin)
  out <- matrix(NA_real_, n_keep, d)
  lp_out <- numeric(n_keep)
  it_keep <- integer(n_keep)

  keep_idx <- 0L

  for (t in seq_len(n_iter)) {
    Sigma_prop <- if (t <= adapt_from) Sigma0 else {
      if (n > 1) M2 / (n - 1) else Sigma0
    }
    prop <- as.numeric(MASS::mvrnorm(1, theta, c_d * Sigma_prop + eps))
    lp_prop <- logpost(prop)
    if (is.finite(lp_prop) && log(runif(1)) < (lp_prop - lp)) {
      theta <- prop
      lp <- lp_prop
      n_accept <- n_accept + 1L
    }

    # Welford online mean/covariance update
    n <- n + 1L
    delta <- theta - mean_t
    mean_t <- mean_t + delta / n
    M2 <- M2 + outer(delta, theta - mean_t)

    if (t > burnin_eff && (t - burnin_eff) %% thin == 0) {
      keep_idx <- keep_idx + 1L
      out[keep_idx, ] <- theta
      lp_out[keep_idx] <- lp
      it_keep[keep_idx] <- total0 + t
    }
  }

  list(
    theta0 = theta0,
    samples = out[seq_len(keep_idx), , drop = FALSE],
    logpost = lp_out[seq_len(keep_idx)],
    iteration = it_keep[seq_len(keep_idx)],
    acceptance = n_accept / (total0 + n_iter),
    n_iter = n_iter, burnin = burnin, thin = thin, seed = seed,
    state = list(theta = theta, logpost = lp, mean = mean_t, M2 = M2,
                 n = n, n_accept = n_accept,
                 n_iter = total0 + n_iter)
  )
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

# Posterior predictive credible intervals for the 12 target summaries.
# Simulates at each retained draw, keeps equilibrium-feasible draws, and
# returns equal-tailed 95% CrI + medians per age group.
predictive_intervals <- function(samples, base_params, data,
                                 target_mode = TARGET_MODE,
                                 n_cores = 1, seed = 1) {
  set.seed(seed)
  run_one <- function(j) {
    tryCatch({
      th <- samples[j, ]
      pm <- tryCatch(build_params(th, base_params), error = function(e) NULL)
      if (is.null(pm)) return(NULL)
      out <- tryCatch(run_sim(pm, data), error = function(e) NULL)
      if (is.null(out)) return(NULL)
      eq <- tryCatch(check_equilibrium(out, target_mode = target_mode),
                     error = function(e) NULL)
      if (is.null(eq) || !eq$pass) return(NULL)
      s <- J_summary_final(out)
      c(s$p_sero, s$N_hat)
    }, error = function(e) NULL)
  }

  if (n_cores > 1 && requireNamespace("parallel", quietly = TRUE)) {
    res <- parallel::mclapply(seq_len(nrow(samples)), run_one,
                              mc.cores = n_cores)
  } else {
    res <- lapply(seq_len(nrow(samples)), run_one)
  }
  keep <- sapply(res, function(z) is.numeric(z) && length(z) == 12L)
  res <- res[keep]
  if (length(res) == 0) stop("no equilibrium-feasible posterior draws")

  mat <- do.call(rbind, res)
  p_mat <- mat[, 1:6, drop = FALSE]
  N_mat <- mat[, 7:12, drop = FALSE]

  out <- data.frame(
    age_group = cal_targets$age_groups,
    p_median = apply(p_mat, 2, median),
    p_lo = apply(p_mat, 2, quantile, probs = 0.025),
    p_hi = apply(p_mat, 2, quantile, probs = 0.975),
    N_median = apply(N_mat, 2, median),
    N_lo = apply(N_mat, 2, quantile, probs = 0.025),
    N_hi = apply(N_mat, 2, quantile, probs = 0.975),
    n_draws_used = length(res)
  )
  attr(out, "predictive_matrix") <- list(p = p_mat, N = N_mat)
  out
}
