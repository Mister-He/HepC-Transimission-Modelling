# =============================================================================
# hmc_calibration.R  —  HMC Calibration for the HCV PWID Compartmental Model
#
# Append this file after setup.R (which sources sim.cpp and defines params,
# data, idx(), y0, etc.).
#
# Parameters estimated (11 total with 10 age groups):
#   theta[1]       = log(beta_scaling_fct - 1), so beta_scaling_fct > 1
#   theta[2:11]    = log_C_contact_scale  row-specific contact scalings
#
# Scaling factors:
#   beta_scaling_fct   = 1 + exp(theta[1])
#   c_true             = beta_scaling_fct / c_composite
#                      where c_composite = base_params$c
#   C_contact_scale[j] = exp(theta[1 + j])   j = 1..10
#
# Likelihood (two-part observation model, J-stratum only, per age group i = 1..10):
#   log(obs_tot[1:10])     ~ Normal(log(total_by_age(theta)), count_log_sd)
#   logit(obs_prev[a])      ~ Normal(logit(q_a(theta)), se_logit[a])
#   se_logit[a]^2           = 1 / (obs_tot[a] * obs_prev[a] * (1 - obs_prev[a]))
#                              + prev_logit_sd^2
#
# Prior:
#   beta_scaling_fct - 1 ~ LogNormal(0, 2)  [theta[1] ~ N(0, 2)]
#   C_contact_scale  ~ LogNormal(0, 1)      [theta[2:11] ~ N(0, 1)]
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

#' Unconstrained theta -> constrained named list.
constrain_theta <- function(theta) {
  n_contact  <- if (exists("N_CONTACT", inherits = TRUE)) {
    get("N_CONTACT", inherits = TRUE)
  } else {
    length(theta) - 1L
  }
  beta_scaling_fct <- 1.0 + exp(theta[1L])
  c_scales   <- exp(theta[1L + seq_len(n_contact)])

  list(
    C_contact_scale = c_scales,
    beta_scaling_fct = beta_scaling_fct,
    log_C_contact_scale = theta[1L + seq_len(n_contact)]
  )
}

#' Build a full parameter list from unconstrained theta
#'
#' All contact rows are scaled directly on the log scale.
build_params_from_theta <- function(theta, base_params) {
  p  <- constrain_theta(theta)
  pm <- base_params

  if (is.null(base_params$c) || length(base_params$c) != 1L ||
      !is.finite(base_params$c) || base_params$c <= 0) {
    stop("base_params$c (c_composite) must be one positive finite value")
  }

  for (j in seq_along(p$C_contact_scale)) {
    pm$C_contact[j, ] <- base_params$C_contact[j, ] * p$C_contact_scale[j]
  }
  pm$beta <- base_params$beta * p$beta_scaling_fct
  pm$c <- p$beta_scaling_fct / base_params$c
  pm
}

#' Vectorised back-transform: theta matrix -> interpretable parameter matrix
#'
#' theta[,1]   log(beta_scaling_fct - 1) -> 1 + exp() = beta_scaling_fct
#' c_true      beta_scaling_fct / c_composite, where c_composite = params$c
#' theta[,2:n] log_C_contact_scale  -> exp() = C_contact_scale
theta_to_orig <- function(samps, c_composite = params$c) {
  if (is.null(dim(samps))) {
    samps <- matrix(samps, nrow = 1L)
  }
  n_contact <- if (exists("N_CONTACT", inherits = TRUE)) {
    get("N_CONTACT", inherits = TRUE)
  } else {
    ncol(samps) - 1L
  }
  c_cols <- do.call(cbind, lapply(seq_len(n_contact), function(j) {
    exp(samps[, 1L + j])
  }))
  colnames(c_cols) <- paste0("C_contact_scale_", seq_len(n_contact))
  if (length(c_composite) != 1L || !is.finite(c_composite) || c_composite <= 0) {
    stop("c_composite must be one positive finite value")
  }
  beta_scaling_fct <- 1.0 + exp(samps[, 1L])
  cbind(
    beta_scaling_fct = beta_scaling_fct,
    c_true = beta_scaling_fct / c_composite,
    c_cols
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

prevalence_logit_sd <- function(obs_prev, obs_tot, extra_sd = 0.25) {
  p <- pmin(pmax(obs_prev, 1e-6), 1 - 1e-6)
  sqrt(1 / pmax(obs_tot * p * (1 - p), 1e-12) + extra_sd^2)
}

get_count_log_sd <- function(data) {
  if (!is.null(data$count_log_sd)) data$count_log_sd else 0.35
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
# Direct log-scale priors:
#   log(beta_scaling_fct - 1) ~ Normal(0, 2)
#   log_C_contact_scale_j ~ Normal(0, 1)
#
log_prior <- function(theta) {
  dnorm(theta[1L], mean = 0.0, sd = 2.0, log = TRUE) +
    sum(dnorm(theta[-1L], mean = 0.0, sd = 1.0, log = TRUE))
}

# ── 3b. Analytical gradient of log-prior ─────────────────────────────────────
#
# Gradients of direct log-scale priors.
#
# Note: not called during sampling (numerical gradients are used); kept for
# reference and unit-testing against finite differences.
#
grad_log_prior_analytical <- function(theta) {
  grad <- numeric(length(theta))
  grad[1L] <- -theta[1L] / 2.0^2
  grad[-1L] <- -theta[-1L]
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

  prev_extra_sd <- if (!is.null(data$prev_logit_sd)) data$prev_logit_sd else 0.25
  prev_sd       <- prevalence_logit_sd(obs_prev, obs_tot, prev_extra_sd)

  ll_count <- sum(dnorm(log(obs_tot), mean = log(obs$total_by_age),
                        sd = get_count_log_sd(data), log = TRUE))
  ll_prev <- sum(dnorm(logit(obs_prev), mean = logit(obs$q_age), sd = prev_sd, log = TRUE))

  ll_count + ll_prev
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
