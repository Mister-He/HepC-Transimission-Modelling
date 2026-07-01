# =============================================================================
# hmc_calibration.R  —  HMC Calibration for the HCV PWID Compartmental Model
#
# Append this file after setup.R (which sources sim.cpp and defines params,
# data, idx(), y0, etc.).
#
# Parameters estimated (20 total with 10 age groups):
#   theta[1:10]    = log_C_contact_scale  row-specific contact scalings
#   theta[11:20]   = log_tot_in_scaling_fct  row-specific total-in scalings
#
# Scaling factors:
#   c_true[k]          = c_composite[k] / tot_in_scaling_fct[k]
#   C_contact_scale[j] = exp(theta[j])   j = 1..10
#   tot_in_scaling_fct[k] = 1 + exp(theta[10 + k])   k = 1..10
#
# Likelihood (joint observation model, J-stratum only, per age group i = 1..10):
#   ALR(q_obs)[a]  ~ Normal(ALR(pi_model)[a], sigma_pop[a])   population composition
#   d2_ell[t, a]   ~ StudentT(nu_shape, 0, sigma_shape)        shape regularization
#   logit(obs_prev[a]) ~ Normal(logit(q_a(theta)), sigma_prev[a]) prevalence
#   sigma_prev[a] = sqrt(se_logit[a]^2 + tau_prev^2), where
#   se_logit[a] = obs_prev_se[a] / (obs_prev[a] * (1 - obs_prev[a]))
#
# Parameter relationships (new c_composite design):
#   c_composite[k]    = tot_in_scaling_fct[k] * c_true[k]   (k = 1..10)
#   C_contact_scale[k] = exp(theta[k])   k = 1..10
#   tot_in_scaling_fct[k] = 1 + exp(theta[10 + k])   k = 1..10
#   lambda1[k]        = lambda3[k] * c_true[k]   (pre-computed before each sim call)
#
# Prior:
#   C_contact_scale  ~ LogNormal(0, 1)      [theta[1:10] ~ N(0, 0.25)]
#   (tot_in_scaling_fct - 1) ~ LogNormal(0, 1)     [theta[11:20] ~ N(0, 0.25)]
#   Note: tot_in_scaling_fct[k] > 1 ensures that c_true[k] < c_composite[k], which is consistent with the model design.
#
# Gradient: central finite differences through the C++ ODE solver
#   Cost per HMC step: (L+1) gradient evals × 2 × N_PARAMS ODE runs
#                    = (L+1) × 40 ODE runs  [plus 1 for current lp]
#   Tip: keep L small (5-10) and use parallel chains.
#
# Outputs:
#   hmc_chains        — raw chain list (all iterations, all chains)
#   post_warmup_list  — post-warmup samples per chain
#   all_samples       — pooled post-warmup matrix (n_samp × 20)
#   diag_table        — R-hat and ESS summary data.frame
#   ppc_out           — posterior predictive replicates
#   plots             — trace, PPC histogram, PPC interval
# =============================================================================

# =============================================================================
# 1.  PARAMETER TRANSFORMS
# =============================================================================

#' Unconstrained theta -> constrained named list.
constrain_theta <- function(theta) {
  if (length(theta) != 2L * N_AGE || any(!is.finite(theta))) {
    stop(sprintf("theta must contain %d finite values", 2L * N_AGE))
  }
  c_scales   <- exp(theta[1:N_AGE])
  tot_in_scaling_fct <- 1.0 + exp(theta[(N_AGE + 1):(2 * N_AGE)])

  list(
    C_contact_scale = c_scales,
    tot_in_scaling_fct = tot_in_scaling_fct,
    log_C_contact_scale = theta[1:N_AGE],
    log_tot_in_scaling_fct = theta[(N_AGE + 1):(2 * N_AGE)]
  )
}

#' Build a full parameter list from unconstrained theta.
#'
#' C_contact_scale[k] scales the contact matrix row k and
#' tot_in_scaling_fct[k] scales the total inflow beta
#' where determines c_true[k] = c_composite[k] / tot_in_scaling_fct[k].
#' lambda1 is then updated to lambda3 * c_true before passing to sim.cpp.
build_params_from_theta <- function(theta, base_params) {
  p  <- constrain_theta(theta)
  pm <- base_params

  if (is.null(base_params$c_composite) ||
      any(!is.finite(base_params$c_composite)) ||
      any(base_params$c_composite <= 0)) {
    stop("base_params$c_composite must be a positive finite vector")
  }

  # Scaling contact matrix
  for (j in seq_along(p$C_contact_scale)) {
    pm$C_contact[j, ] <- base_params$C_contact[j, ] * p$C_contact_scale[j]
  }

  # Scaling total inflow beta
  pm$beta <- base_params$beta * p$tot_in_scaling_fct

  c_true    <- base_params$c_composite / p$tot_in_scaling_fct
  pm$lambda1 <- base_params$lambda3 * c_true
  pm
}

#' Vectorised back-transform: theta matrix -> interpretable parameter matrix.
#'
#' theta[,1:N_AGE] log(C_contact_scale) -> exp() = C_contact_scale
#' theta[,N_AGE + (1:N_AGE)] log(tot_in_scaling_fct - 1)
#'   -> 1 + exp() = tot_in_scaling_fct
theta_to_orig <- function(samps) {
  if (is.null(dim(samps))) {
    samps <- matrix(samps, nrow = 1L)
  }
  contact_cols <- do.call(cbind, lapply(seq_len(N_AGE), function(j) {
    exp(samps[, j])
  }))
  colnames(contact_cols) <- paste0("C_contact_scale_", seq_len(N_AGE))
  tot_in_cols <- do.call(cbind, lapply(seq_len(N_AGE), function(j) {
    1.0 + exp(samps[, N_AGE + j])
  }))
  colnames(tot_in_cols) <- paste0("tot_in_scaling_fct_", seq_len(N_AGE))

  cbind(
    contact_cols,
    tot_in_cols
  )
}


# =============================================================================
# 2.  AGE-STRUCTURED MODEL SUMMARIES FROM FINAL ODE STATE
# =============================================================================

logit <- function(p) {
  p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
  log(p / (1 - p))
}

inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

prevalence_logit_sd <- function(obs_prev, obs_tot = NULL, extra_sd = 0.25, obs_se = NULL) {
  p <- pmin(pmax(obs_prev, 1e-6), 1 - 1e-6)
  if (!is.null(obs_se)) {
    se_logit <- pmax(obs_se, 1e-9) / (p * (1 - p))
  } else {
    se_logit <- sqrt(1 / pmax(obs_tot * p * (1 - p), 1e-12))
  }
  sqrt(se_logit^2 + extra_sd^2)
}


#' Compute model-implied age-group totals and HCV prevalence from the
#' J-stratum (s = 1) at the final ODE state.
compute_age_quantities <- function(y_final) {
  n_age <- length(obs_tot)
  age_total <- numeric(n_age)
  age_pos <- numeric(n_age)

  for (i in seq_len(n_age) - 1L) {
    total_i <- 0.0
    pos_i <- 0.0
    for (k in 1:4) {
      total_i <- total_i + y_final[idx(s = 1L, k = k, h = 0L, i = i)]
      total_i <- total_i + y_final[idx(s = 1L, k = k, h = 1L, i = i)]
      total_i <- total_i + y_final[idx(s = 1L, k = k, h = 2L, i = i)]
      total_i <- total_i + y_final[idx(s = 1L, k = k, h = 3L, i = i)]

      # Preserve the current positive-state definition used in the model.
      pos_i <- pos_i + y_final[idx(s = 1L, k = k, h = 1L, i = i)]
      pos_i <- pos_i + y_final[idx(s = 1L, k = k, h = 2L, i = i)]
      pos_i <- pos_i + y_final[idx(s = 1L, k = k, h = 3L, i = i)]
    }
    age_total[i + 1L] <- total_i
    age_pos[i + 1L] <- pos_i
  }

  n_model_total <- sum(age_total)
  if (!is.finite(n_model_total) || n_model_total <= 0) return(NULL)
  if (any(!is.finite(age_total)) || any(!is.finite(age_pos))) return(NULL)
  if (any(age_total <= 0)) return(NULL)
  if (any(age_pos < 0) || any(age_pos > age_total)) return(NULL)

  list(
    total_by_age = age_total,
    pos_by_age   = age_pos,
    n_model_total = n_model_total,
    p_age = age_total / n_model_total,
    q_age = pmin(pmax(age_pos / age_total, 1e-12), 1 - 1e-12)
  )
}


# =============================================================================
# 3.  LOG-POSTERIOR AND GRADIENT
# =============================================================================

# ── 3a. Log-prior ────────────────────────────────────────────────────────────
#
# theta[1:10]  = log(C_contact_scale[k]) ~ Normal(0, 1)
# theta[11:20] = log(tot_in_scaling_fct[k] - 1) ~ Normal(0, 1)
CONTACT_PRIOR_MEANS <- c(log(0.5), rep(0.05, 4), rep(0.75, 5))
TOT_IN_PRIOR_MEANS  <- rep(0.0, N_AGE)

log_prior <- function(theta) {
    sum(dnorm(theta[1:N_AGE], mean = CONTACT_PRIOR_MEANS, sd = 0.25, log = TRUE)) +
    sum(dnorm(theta[(N_AGE + 1):(2 * N_AGE)], mean = TOT_IN_PRIOR_MEANS, sd = 0.25, log = TRUE))
}

# ── 3b. Population composition log-likelihood (ALR-normal) ───────────────────
#
# Compares observed age composition q_obs[a] = Y_pop_obs[a] / sum(Y_pop_obs)
# with model-implied composition pi_model[a] = N_model[a] / sum(N_model)
# via an additive log-ratio (ALR) normal likelihood, reference group = A.
#
# For a = 1, ..., A-1:
#   z_obs[a]   = log(q_obs[a]  / q_obs[A])
#   z_model[a] = log(pi[a]     / pi[A])
#   z_obs[a]  ~ Normal(z_model[a], sigma_pop[a])
#
compute_population_composition_loglik <- function(pi_model, q_obs, sigma_pop) {
  # Choose last age group as reference
  n_age <- length(pi_model)
  eps   <- 1e-12

  pi_safe <- pmax(pi_model, eps); pi_safe <- pi_safe / sum(pi_safe)
  q_safe  <- pmax(q_obs,    eps); q_safe  <- q_safe  / sum(q_safe)

  if (length(sigma_pop) == 1L) sigma_pop <- rep(sigma_pop, n_age)

  ref <- n_age
  ll  <- 0.0
  for (a in seq_len(n_age - 1L)) {
    z_obs   <- log(q_safe[a]  / q_safe[ref])
    z_model <- log(pi_safe[a] / pi_safe[ref])
    ll <- ll + dnorm(z_obs, mean = z_model, sd = sigma_pop[a], log = TRUE)
  }
  ll
}

# ── 3c. Age-shape log-prior (Student-t on second differences of centered log π)
#
# Regularizes curvature of the log age-composition profile without forcing
# smoothness. For each internal age group a = 2, ..., A-1:
#   ell[a] = log(pi[a]) - mean(log(pi))
#   d2[a]  = ell[a+1] - 2*ell[a] + ell[a-1]
#   d2[a] ~ StudentT(nu_shape, 0, sigma_shape)
#
compute_age_shape_logprior <- function(pi_model, nu_shape = 4L, sigma_shape = 0.3) {
  n_age <- length(pi_model)
  if (n_age < 3L) return(0.0)

  eps     <- 1e-12
  pi_safe <- pmax(pi_model, eps); pi_safe <- pi_safe / sum(pi_safe)

  log_pi <- log(pi_safe)
  ell    <- log_pi - mean(log_pi)

  lp <- 0.0
  for (a in 2L:(n_age - 1L)) {
    d2 <- ell[a + 1L] - 2.0 * ell[a] + ell[a - 1L]
    lp <- lp + dt(d2 / sigma_shape, df = nu_shape, log = TRUE) - log(sigma_shape)
  }
  lp
}

# ── 3d. Prevalence log-likelihood (logit-normal with delta-method SE) ─────────
#
# logit(prev_obs[a]) ~ Normal(logit(p_model[a]), sqrt(se_logit[a]^2 + tau^2))
# se_logit computed from obs_prev_se via delta method: se / (p * (1-p))
#
compute_prevalence_loglik <- function(q_age, obs_prev, obs_prev_se, prev_extra_sd = 0.25) {
  prev_sd <- prevalence_logit_sd(obs_prev, obs_se = obs_prev_se, extra_sd = prev_extra_sd)
  sum(dnorm(logit(obs_prev), mean = logit(q_age), sd = prev_sd, log = TRUE))
}

# ── 3e. Analytical gradient of log-prior ─────────────────────────────────────
#
# Gradients of direct log-scale priors.
#
# Note: not called during sampling (numerical gradients are used); kept for
# reference and unit-testing against finite differences.
#
grad_log_prior_analytical <- function(theta) {
  grad <- numeric(length(theta))
  grad[1:N_AGE] <- -(theta[1:N_AGE] - CONTACT_PRIOR_MEANS)
  grad[(N_AGE + 1):(2 * N_AGE)] <- -(theta[(N_AGE + 1):(2 * N_AGE)] - TOT_IN_PRIOR_MEANS)
  grad
}

# ── 3f. Log-likelihood ───────────────────────────────────────────────────────
log_likelihood <- function(theta, base_params, data) {
  pm  <- build_params_from_theta(theta, base_params)
  out <- run_sim(pm, data)

  if (!is.matrix(out) || nrow(out) == 0L) return(-Inf)
  if (any(!is.finite(out)))               return(-Inf)

  y_final <- as.numeric(out[nrow(out), -1L])
  obs     <- compute_age_quantities(y_final)
  if (is.null(obs)) return(-Inf)

  q_obs    <- obs_tot / sum(obs_tot)
  pi_model <- obs$total_by_age / sum(obs$total_by_age)

  sigma_pop   <- if (!is.null(data$sigma_pop))   data$sigma_pop   else rep(0.05, length(obs_tot))
  sigma_shape <- if (!is.null(data$sigma_shape)) data$sigma_shape else 0.3
  nu_shape    <- if (!is.null(data$nu_shape))    data$nu_shape    else 4L
  prev_extra_sd <- if (!is.null(data$prev_logit_sd)) data$prev_logit_sd else 0.25

  ll_pop   <- compute_population_composition_loglik(pi_model, q_obs, sigma_pop)
  lp_shape <- compute_age_shape_logprior(pi_model, nu_shape, sigma_shape)
  ll_prev  <- compute_prevalence_loglik(obs$q_age, obs_prev, obs_prev_se, prev_extra_sd)

  if (!is.finite(ll_pop) || !is.finite(lp_shape) || !is.finite(ll_prev)) return(-Inf)

  ll_pop + ll_prev + lp_shape
}

# ── 3d. Log-posterior ────────────────────────────────────────────────────────
log_posterior <- function(theta, base_params, data) {
  lp <- log_prior(theta)
  if (!is.finite(lp)) return(-Inf)

  ll <- tryCatch(
    log_likelihood(theta, base_params, data),
    error = function(e) -Inf
  )
  if (!is.finite(ll)) return(-Inf)
  lp + ll
}

# ── 3e. Numerical gradient (central finite differences) ───────────────────────
#
# Analytical gradients are not available because the likelihood is computed by
# a black-box C++ ODE solver. We therefore use central finite differences:
#
#   ∂ log π / ∂ theta_j  ≈
#       [log π(theta + eps·e_j) - log π(theta - eps·e_j)] / (2·eps)
#
# This requires 2 × N_PARAMS ODE runs per gradient evaluation.
# With L leapfrog steps, each HMC proposal costs ≈ (L + 1) × 2 × N_PARAMS ODE runs.
#
# Fallback: if either perturbed evaluation is -Inf, the component is set to 0
# (safe — the chain will reject or stay put, not diverge).
#
grad_log_posterior <- function(theta, base_params, data, eps = 1e-4) {
  n    <- length(theta)
  grad <- numeric(n)
  lp0  <- NULL   # lazily evaluated if needed for fallback

  for (j in seq_len(n)) {
    theta_p    <- theta; theta_p[j] <- theta[j] + eps
    theta_m    <- theta; theta_m[j] <- theta[j] - eps

    lp_p <- log_posterior(theta_p, base_params, data)
    lp_m <- log_posterior(theta_m, base_params, data)

    if (is.finite(lp_p) && is.finite(lp_m)) {
      grad[j] <- (lp_p - lp_m) / (2 * eps)
    } else if (is.finite(lp_p)) {
      if (is.null(lp0)) lp0 <- log_posterior(theta, base_params, data)
      grad[j] <- (lp_p - lp0) / eps
    } else if (is.finite(lp_m)) {
      if (is.null(lp0)) lp0 <- log_posterior(theta, base_params, data)
      grad[j] <- (lp0 - lp_m) / eps
    } else {
      grad[j] <- 0.0   # both sides infeasible — zero gradient (safe)
    }

  }
  grad
}

# =============================================================================
# 4.  LEAPFROG INTEGRATOR
# =============================================================================
#
# Störmer–Verlet / leapfrog scheme for Hamiltonian H(θ, r) = U(θ) + K(r):
#   U(θ) = -log π(θ)   (potential energy)
#   K(r) = ½ rᵀr       (kinetic energy, unit mass matrix)
#
# Algorithm:
#   r_{½}     ← r₀ + (ε/2) ∇log π(θ₀)          half-step momentum
#   θ_l       ← θ_{l-1} + ε · r_{l-½}            full position step  } ×L
#   r_{l+½}   ← r_{l-½} + ε · ∇log π(θ_l)       full momentum step  }
#   r_L       ← r_{L-½} + (ε/2) ∇log π(θ_L)     final half-step
#
# Total gradient evaluations: L + 1
#
leapfrog <- function(theta, r, grad_fn, eps, L) {
  g <- grad_fn(theta)             # ∇log π(θ₀)
  r <- r + 0.5 * eps * g         # half-step momentum

  for (l in seq_len(L)) {
    theta <- theta + eps * r      # full position update
    g     <- grad_fn(theta)       # ∇log π(θ_l)
    if (l < L) {
      r <- r + eps * g            # full momentum update (not last step)
    }
  }
  r <- r + 0.5 * eps * g         # final half-step (g = ∇log π(θ_L))

  list(theta = theta, r = r)
}


# =============================================================================
# 5.  HMC CHAIN WITH DUAL-AVERAGING STEP-SIZE ADAPTATION
# =============================================================================
#
# Step-size adaptation follows Hoffman & Gelman (2014) Algorithm 5:
#   H̄_m  ← (1 - 1/(m+t₀)) H̄_{m-1} + (δ - α_m)/(m+t₀)
#   ε_m   ← exp(μ - √m/γ · H̄_m)
#   ε̄_m  ← exp(m^{-κ} log ε_m + (1-m^{-κ}) log ε̄_{m-1})
# After warmup: fix ε = ε̄_{n_warmup}.
#
run_hmc_chain <- function(theta_init,
                           n_iter,
                           n_warmup,
                           eps_init    = 0.01,
                           L           = 10L,
                           adapt_delta = 0.65,
                           base_params,
                           data,
                           seed        = NULL,
                           chain_id    = 1L) {
  if (!is.null(seed)) set.seed(seed)

  n_p      <- length(theta_init)
  samples  <- matrix(NA_real_, n_iter, n_p,
                     dimnames = list(NULL, param_names_log))
  accepted  <- logical(n_iter)
  lp_trace  <- numeric(n_iter)
  eps_trace <- numeric(n_iter)

  theta <- theta_init
  lp    <- log_posterior(theta, base_params, data)

  # Dual-averaging state
  eps     <- eps_init
  eps_bar <- eps_init
  H_bar   <- 0.0
  mu      <- log(10 * eps_init)   # μ = log(10ε₀)
  gamma   <- 0.05
  t0      <- 10
  kappa   <- 0.75

  cat(sprintf(
    "[Chain %d] Starting  n_iter=%d | n_warmup=%d | eps0=%.4f | L=%d\n",
    chain_id, n_iter, n_warmup, eps_init, L))

  for (iter in seq_len(n_iter)) {
    # Sample momenta from N(0, I)
    r0 <- rnorm(n_p)
    K0 <- 0.5 * sum(r0^2)

    # Leapfrog proposal
    grad_fn <- function(t) grad_log_posterior(t, base_params, data)
    prop    <- leapfrog(theta, r0, grad_fn, eps, L)

    lp_prop <- log_posterior(prop$theta, base_params, data)
    K_prop  <- 0.5 * sum(prop$r^2)

    # Metropolis-Hastings acceptance
    log_alpha <- (lp_prop - K_prop) - (lp - K0)
    alpha     <- if (is.finite(log_alpha)) min(1.0, exp(log_alpha)) else 0.0

    if (runif(1L) < alpha) {
      theta <- prop$theta
      lp    <- lp_prop
      accepted[iter] <- TRUE
    }

    samples[iter, ]  <- theta
    lp_trace[iter]   <- lp
    eps_trace[iter]  <- eps

    # ── Dual-averaging adaptation (warmup only) ────────────────────────────
    if (iter <= n_warmup) {
      H_bar   <- (1 - 1 / (iter + t0)) * H_bar +
                 (adapt_delta - alpha) / (iter + t0)
      log_eps <- mu - sqrt(iter) / gamma * H_bar
      eps     <- exp(log_eps)
      eps_bar <- exp(iter^(-kappa) * log_eps +
                     (1 - iter^(-kappa)) * log(eps_bar))

      if (iter == n_warmup) {
        eps <- eps_bar
        cat(sprintf("[Chain %d] Warmup done. Adapted eps = %.5f\n",
                    chain_id, eps))
      }
    }

    # Progress report every 100 iterations
    if (iter %% 100L == 0L) {
      phase   <- if (iter <= n_warmup) "warmup" else "sample"
      win_acc <- mean(accepted[max(1L, iter - 99L):iter])
      cat(sprintf("[Chain %d] iter=%4d (%s) | lp=%8.2f | acc=%.2f | eps=%.5f\n",
                  chain_id, iter, phase, lp, win_acc, eps))
    }
  }
  
  idx_post <- seq(n_warmup + 1L, n_iter)
  list(
    samples     = samples[idx_post, , drop = FALSE],  # post-warmup only
    samples_all = samples,                             # full trace (for plots)
    lp_trace    = lp_trace,
    eps_trace   = eps_trace,
    accepted    = accepted,
    acc_rate    = mean(accepted[idx_post]),
    eps_final   = eps,
    n_warmup    = n_warmup
  )
}
