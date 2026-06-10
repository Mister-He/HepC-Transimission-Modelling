# =============================================================================
# hmc_calibration.R  —  HMC Calibration for the HCV PWID Compartmental Model
#
# Append this file after setup.R (which sources sim.cpp and defines params,
# data, idx(), y0, etc.).
#
# Parameters estimated (12 total, all in log-space for positivity):
#   theta[1]       = log(beta_scale)        scalar multiplier on beta[1:9]
#   theta[2]       = log(delta)             shift parameter for shifted Gamma prior
#   theta[3]       = log(alpha)             shape parameter for Gamma hyperprior
#   theta[4]       = log(beta_rate)         rate parameter for Gamma hyperprior
#   theta[5:12]    = log(C_contact_scale[1:8])  8 age-specific row scalings
#
# Scaling factors generated:
#   C_contact_scale[1:8] = delta + Gamma(alpha, beta_rate)
#   C_contact_scale[9]    = 1  (reference row, fixed for identifiability)
#
# Likelihood (three-step observation model, J-stratum only, per age group i = 1..9):
#   N_total_obs            ~ LogNormal(log(N_total_model), sigma_N)
#   obs_tot[1:9] | N_total ~ Multinomial(N_total_obs, p_1:9(theta))
#   obs_pos[a] | obs_tot[a] ~ BetaBinomial(obs_tot[a], q_a(theta), phi)
#
# Prior: LogNormal hyperpriors on (beta_scale, delta, alpha, beta_rate)
#   beta_scale ~ LogNormal(0, 1)
#   delta      ~ LogNormal(log(0.25), 0.5)   [shift/floor for row scalings]
#   alpha      ~ LogNormal(log(2.5), 0.3)    [Gamma shape hyperparameter]
#   beta_rate  ~ LogNormal(log(1.0), 0.3)    [Gamma rate hyperparameter]
#
# The 8 free row scalings follow a shifted Gamma hierarchy, with the 9th row
# fixed to 1 as the reference row for identifiability.
#
# Gradient: central finite differences through the C++ ODE solver
#   Cost per HMC step: (L+1) gradient evals × 2 × N_PARAMS ODE runs
#                    = (L+1) × 24 ODE runs  [plus 1 for current lp]
#   Tip: keep L small (5-10) and use parallel chains.
#
# Outputs:
#   hmc_chains        — raw chain list (all iterations, all chains)
#   post_warmup_list  — post-warmup samples per chain
#   all_samples       — pooled post-warmup matrix (n_samp × 12)
#   diag_table        — R-hat and ESS summary data.frame
#   ppc_out           — posterior predictive replicates
#   plots             — trace, PPC histogram, PPC interval
# =============================================================================

 =============================================================================
# 1.  PARAMETER TRANSFORMS
# =============================================================================

#' Unconstrained theta (log-space) -> constrained named list
#' 
#' Generates 8 free contact scaling factors from a shifted Gamma hierarchy.
#' The 9th row is held at 1 as a reference for identifiability.
constrain_theta <- function(theta) {
  beta_scale <- exp(theta[1L])
  delta      <- exp(theta[2L])
  alpha      <- exp(theta[3L])
  beta_rate  <- exp(theta[4L])

  # 8 estimated row scalings; row 9 remains the reference row (scale = 1)
  c_scales <- c(delta + exp(theta[5:12]), 1.0)
  
  list(
    beta_scale      = beta_scale,
    C_contact_scale = c_scales,
    delta           = delta,
    alpha           = alpha,
    beta_rate       = beta_rate
  )
}

#' Build a full parameter list from unconstrained theta
#'
#' The 8 free row scaling factors are hierarchical draws on the log scale,
#' while the 9th row remains fixed at 1.
build_params_from_theta <- function(theta, base_params) {
  p  <- constrain_theta(theta)
  pm <- base_params

  # Scale beta uniformly across all age groups
  pm$beta <- base_params$beta * p$beta_scale

  # Scale rows 1:8; row 9 is the reference row and stays unchanged.
  for (j in 1:8) {
    pm$C_contact[j, ] <- base_params$C_contact[j, ] * p$C_contact_scale[j]
  }
  pm
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
  age_total <- numeric(9L)
  age_pos <- numeric(9L)

  for (i in 0:8) {
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
# Prior specification on original scale:
#   beta_scale ~ LogNormal(0, 1)                 [median 1]
#   delta      ~ LogNormal(log(0.25), 0.5)       [positive shift/floor]
#   alpha      ~ LogNormal(log(2.5), 0.3)        [Gamma shape hyperparameter]
#   beta_rate  ~ LogNormal(log(1.0), 0.3)        [Gamma rate hyperparameter]
#
# The 8 free row scalings are modeled as delta + Gamma(alpha, beta_rate),
# which places the prior mass for the resulting scales around 1–2.
#
# Change-of-variables: theta_j = log(X_j), X_j = exp(theta_j)
#   log p(theta_j) = log p_X(exp(theta_j)) + theta_j
#
log_prior <- function(theta) {
  beta_scale <- exp(theta[1L])
  delta      <- exp(theta[2L])
  alpha      <- exp(theta[3L])
  beta_rate  <- exp(theta[4L])
  excess     <- exp(theta[5:12])
  
  # LogNormal priors on each parameter (using dnorm on log scale + Jacobian)
  lp_beta_scale <- dnorm(theta[1L], mean =  0.0, sd = 1.0, log = TRUE)
  lp_delta      <- dnorm(theta[2L], mean = log(0.25), sd = 0.5, log = TRUE)
  lp_alpha      <- dnorm(theta[3L], mean = log(2.5), sd = 0.3, log = TRUE)
  lp_beta_rate  <- dnorm(theta[4L], mean = log(1.0), sd = 0.3, log = TRUE)
  lp_contact    <- sum(dgamma(excess, shape = alpha, rate = beta_rate, log = TRUE) + theta[5:12])
  
  # Sum across all parameters
  lp_beta_scale + lp_delta + lp_alpha + lp_beta_rate + lp_contact
}

# ── 3b. Analytical gradient of log-prior ─────────────────────────────────────
#
# With shifted Gamma priors on the row scalings:
#   d/d theta_j [dgamma(exp(theta_j), shape, rate, log=TRUE) + theta_j]
#     = alpha - beta_rate * exp(theta_j)
#
grad_log_prior_analytical <- function(theta) {
  grad <- numeric(12L)
  grad[1L] <- -(theta[1L] - 0.0) / 1.0^2
  grad[2L] <- -(theta[2L] - log(0.25)) / 0.5^2
  grad[3L] <- -(theta[3L] - log(2.5)) / 0.3^2
  grad[4L] <- -(theta[4L] - log(1.0)) / 0.3^2
  alpha     <- exp(theta[3L])
  beta_rate <- exp(theta[4L])
  grad[5:12] <- alpha - beta_rate * exp(theta[5:12])
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

  sigma_N <- if (!is.null(data$sigma_N)) data$sigma_N else 0.10
  phi     <- if (!is.null(data$phi_overdisp)) data$phi_overdisp else 50.0

  ll_total <- dlnorm(N_total_obs, meanlog = log(obs$n_model_total), sdlog = sigma_N, log = TRUE)
  ll_age   <- dmultinom(obs_tot, prob = obs$p_age, log = TRUE)
  ll_case  <- log_beta_binomial(obs_pos, size = obs_tot, prob = obs$q_age, phi = phi)

  ll_total + ll_age + ll_case
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
# This requires 2 × N_PARAMS = 20 ODE runs per gradient evaluation.
# With L leapfrog steps, each HMC proposal costs ≈ (L + 1) × 20 ODE runs.
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
