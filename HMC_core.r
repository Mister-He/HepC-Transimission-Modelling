# =============================================================================
# hmc_calibration.R  —  HMC Calibration for the HCV PWID Compartmental Model
#
# Append this file after setup.R (which sources sim.cpp and defines params,
# data, idx(), y0, etc.).
#
# Parameters estimated (11 total with 10 age groups):
#   theta[1]       = mu_hier           log-scale mean of contact row scalings
#   theta[2]       = log(sigma_hier)   log-scale std dev of contact row scalings
#   theta[3:11]    = eta[1:9]          standardised deviates (non-centred)
#
# Scaling factors generated (non-centred log-normal hierarchy):
#   C_contact_scale[j] = exp(mu_hier + sigma_hier * eta[j])   j = 1..9
#   C_contact_scale[10] = 1  (reference row, fixed for identifiability)
#
# Likelihood (two-step observation model, J-stratum only, per age group i = 1..10):
#   obs_tot[1:10]          ~ Multinomial(N_age_obs, p_1:10(theta))
#   obs_pos[a] | obs_tot[a] ~ BetaBinomial(obs_tot[a], q_a(theta), phi)
#
# Prior:
#   mu_hier     ~ Normal(0, 1)              [log-scale mean; exp(0)=1 reference]
#   sigma_hier  ~ LogNormal(log(0.5), 0.5)  [theta[2] ~ N(log(0.5), 0.5)]
#   eta[j]      ~ Normal(0, 1)  j = 1..9   [non-centred deviates]
#
# Non-centred reparameterisation removes the delta+excess ridge and decouples
# the hierarchy hyperparameters from the unit-level contact row scalings.
#
# Gradient: central finite differences through the C++ ODE solver
#   Cost per HMC step: (L+1) gradient evals × 2 × N_PARAMS ODE runs
#                    = (L+1) × 22 ODE runs  [plus 1 for current lp]
#   Tip: keep L small (5-10) and use parallel chains.
#
# Outputs:
#   hmc_chains        — raw chain list (all iterations, all chains)
#   post_warmup_list  — post-warmup samples per chain
#   all_samples       — pooled post-warmup matrix (n_samp × 11)
#   diag_table        — R-hat and ESS summary data.frame
#   ppc_out           — posterior predictive replicates
#   plots             — trace, PPC histogram, PPC interval
# =============================================================================

# =============================================================================
# 1.  PARAMETER TRANSFORMS
# =============================================================================

#' Unconstrained theta -> constrained named list (non-centred log-normal hierarchy)
#'
#' C_contact_scale[j] = exp(mu_hier + sigma_hier * eta[j]) for all free rows.
#' The final row is held at 1 as the reference for identifiability.
constrain_theta <- function(theta) {
  mu_hier    <- theta[1L]
  sigma_hier <- exp(theta[2L])
  eta        <- theta[-c(1L, 2L)]

  c_scales <- c(exp(mu_hier + sigma_hier * eta), 1.0)

  list(
    C_contact_scale = c_scales,
    mu_hier         = mu_hier,
    sigma_hier      = sigma_hier,
    eta             = eta
  )
}

#' Build a full parameter list from unconstrained theta
#'
#' The free row scaling factors are hierarchical draws on the log scale,
#' while the final row remains fixed at 1.
build_params_from_theta <- function(theta, base_params) {
  p  <- constrain_theta(theta)
  pm <- base_params

  # Scale all but the final row; the final row is the reference row.
  for (j in seq_along(p$eta)) {
    pm$C_contact[j, ] <- base_params$C_contact[j, ] * p$C_contact_scale[j]
  }
  pm
}

#' Vectorised back-transform: theta matrix -> interpretable parameter matrix
#'
#' theta[,1]   mu_hier           -> identity         = mu_hier (log-scale mean)
#' theta[,2]   log(sigma_hier)   -> exp()           = sigma_hier
#' theta[,3:n] eta               -> exp(mu + sig*eta) = C_contact_scale
theta_to_orig <- function(samps) {
  mu    <- samps[, 1L]
  sigma <- exp(samps[, 2L])
  n_free <- ncol(samps) - 2L
  c_cols <- do.call(cbind, lapply(seq_len(n_free), function(j) {
    exp(mu + sigma * samps[, 2L + j])
  }))
  colnames(c_cols) <- paste0("C_contact_scale_", seq_len(n_free))
  cbind(mu_hier    = samps[, 1L],
        sigma_hier = sigma,
        c_cols)
}


# =============================================================================
# 2.  POISSON MEANS FROM FINAL ODE STATE
# =============================================================================

log_beta_binomial <- function(x, size, prob, phi) {
  if (length(x) != length(size) || length(x) != length(prob)) {
    stop("x, size, and prob must have the same length")
  }
  if (!is.finite(phi) || phi <= 0) return(-Inf)
  if (any(!is.finite(x)) || any(!is.finite(size)) || any(!is.finite(prob))) return(-Inf)
  if (any(x < 0) || any(size < 0) || any(x > size)) return(-Inf)

  prob <- pmin(pmax(prob, 1e-12), 1 - 1e-12)
  a <- prob * phi
  b <- (1 - prob) * phi

  sum(lchoose(size, x) + lbeta(x + a, size - x + b) - lbeta(a, b))
}

#' Compute model-implied age-group totals and HCV-positive counts from the
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
# Non-centred log-normal hierarchy (all theta are unconstrained or log-transformed):
#   mu_hier     ~ Normal(0, 1)              → theta[1] ~ N(0, 1)
#   sigma_hier  ~ LogNormal(log(0.5), 0.5)  → theta[2] ~ N(log(0.5), 0.5)
#   eta[j]      ~ Normal(0, 1)              → theta[3:n] ~ N(0, 1)
#
log_prior <- function(theta) {
  lp_mu_hier    <- dnorm(theta[1L],    mean = 0.0,      sd = 1.0, log = TRUE)
  lp_log_sigma  <- dnorm(theta[2L],    mean = log(0.5), sd = 0.5, log = TRUE)
  lp_eta        <- sum(dnorm(theta[-c(1L, 2L)], mean = 0.0, sd = 1.0, log = TRUE))

  lp_mu_hier + lp_log_sigma + lp_eta
}

# ── 3b. Analytical gradient of log-prior ─────────────────────────────────────
#
# Gradients of the non-centred log-normal hierarchy:
#   d/d theta[1]    = -(theta[1] - 0) / 1^2
#   d/d theta[2]    = -(theta[2] - log(0.5)) / 0.5^2
#   d/d theta[j]    = -theta[j]   for j = 3..n  (eta ~ N(0,1))
#
# Note: not called during sampling (numerical gradients are used); kept for
# reference and unit-testing against finite differences.
#
grad_log_prior_analytical <- function(theta) {
  grad <- numeric(length(theta))
  grad[1L]   <- -(theta[1L] - 0.0) / 1.0^2
  grad[2L]   <- -(theta[2L] - log(0.5)) / 0.5^2
  grad[-c(1L, 2L)] <- -theta[-c(1L, 2L)]
  grad
}

# ── 3c. Log-likelihood ───────────────────────────────────────────────────────
log_likelihood <- function(theta, base_params, data) {
  pm  <- build_params_from_theta(theta, base_params)
  out <- run_sim(pm, data)

  if (!is.matrix(out) || nrow(out) == 0L) return(-Inf)
  if (any(!is.finite(out)))               return(-Inf)

  y_final <- as.numeric(out[nrow(out), -1L])   # drop time column
  obs     <- compute_age_quantities(y_final)
  if (is.null(obs)) return(-Inf)

  phi     <- if (!is.null(data$phi_overdisp)) data$phi_overdisp else 50.0

  ll_age   <- dmultinom(obs_tot, prob = obs$p_age, log = TRUE)
  ll_case  <- log_beta_binomial(obs_pos, size = obs_tot, prob = obs$q_age, phi = phi)

  ll_age + ll_case
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
