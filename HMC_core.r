# =============================================================================
# hmc_calibration.R  —  HMC Calibration for the HCV PWID Compartmental Model
#
# Append this file after setup.R (which sources sim.cpp and defines params,
# data, idx(), y0, etc.).
#
# Model architecture:
#   1. Two age-varying positive scaling functions modify the C++ ODE model:
#        C_contact_scale[k]    = exp((SPLINE_B %*% alpha)[k])
#        tot_in_scaling_fct[k] = exp((SPLINE_B %*% gamma)[k])
#   2. A cubic B-spline basis with 3 internal knots maps each 7-dimensional
#      coefficient vector to N_AGE = 10 age-specific values.
#   3. The final ODE state is compared with observed age composition and
#      age-specific HCV prevalence. No additional age-shape penalty is used.
#   4. HMC samples the spline coefficients using finite-difference gradients
#      through the black-box C++ ODE solver.
#
# Parameters estimated (2*SPLINE_K = 14 total with the current basis):
#   theta[1:K]       = alpha  coefficients for log_C_contact_scale
#   theta[K+1:2K]    = gamma  coefficients for log_tot_in_scaling_fct
#   SPLINE_K = 3 internal knots + cubic degree 3 + intercept 1 = 7.
#
# Scaling relationships:
#   c_true[k]          = c_composite[k] / tot_in_scaling_fct[k]
#   c_composite[k]     = tot_in_scaling_fct[k] * c_true[k]
#   lambda1[k]         = lambda3[k] * c_true[k]
#   Each relationship is updated before every simulation call.
#
# Joint observation model (J-stratum only, per age group i = 1..10):
#   ALR(q_obs)[a] ~ Normal(ALR(pi_model)[a], sigma_pop[a])  population composition
#   obs_pos[a]    ~ Binomial(obs_tot[a], q_a(theta))         HCV prevalence
#   obs_pos[a] = round(obs_prev[a] * obs_tot[a]) HCV positives per age group.
#
# Prior (centered P-spline / RW2 on deviations from a reference curve):
#   mean0 is obtained by projecting the original per-age log prior medians
#   (CONTACT_PRIOR_MEANS or TOT_IN_PRIOR_MEANS) onto SPLINE_B.
#   delta = coef - mean0
#   delta[1] ~ Normal(0, SPLINE_ANCHOR_SD)
#   delta[2] ~ Normal(0, SPLINE_ANCHOR_SD)
#   delta[k] ~ Normal(2*delta[k-1] - delta[k-2], SPLINE_RW_SD), k >= 3
#   Thus the complete projected prior curve is the prior center, while smooth
#   departures from that curve are allowed.
#
# Gradient: central finite differences through the C++ ODE solver
#   N_PARAMS = 14, so one gradient evaluation normally requires
#   2 × N_PARAMS = 28 ODE runs.
#   One HMC proposal with L leapfrog steps costs approximately
#   (L + 1) × 28 ODE runs, plus one ODE run for the proposed log posterior.
#   Tip: keep L small (5-10) and use parallel chains.
#
# Outputs:
#   hmc_chains        — raw chain list (all iterations, all chains)
#   post_warmup_list  — post-warmup samples per chain
#   all_samples       — pooled post-warmup matrix (n_samp × 14)
#   diag_table        — R-hat and ESS summary data.frame
#   ppc_out           — posterior predictive replicates
#   plots             — trace, PPC histogram, PPC interval
# =============================================================================

# =============================================================================
# 1.  PARAMETER TRANSFORMS
# =============================================================================

#' Unconstrained theta (spline coefficients) -> constrained named list.
#'
#' theta = c(alpha, gamma) with SPLINE_K coefficients each. The N_AGE age values
#' are reconstructed via the fixed B-spline basis: log_vec = SPLINE_B %*% coef.
constrain_theta <- function(theta) {
  if (length(theta) != 2L * SPLINE_K || any(!is.finite(theta))) {
    stop(sprintf("theta must contain %d finite values", 2L * SPLINE_K))
  }
  alpha <- theta[1:SPLINE_K]
  gamma <- theta[(SPLINE_K + 1):(2 * SPLINE_K)]

  log_C_contact_scale    <- as.numeric(SPLINE_B %*% alpha)
  log_tot_in_scaling_fct <- as.numeric(SPLINE_B %*% gamma)

  list(
    C_contact_scale = exp(log_C_contact_scale),
    tot_in_scaling_fct = exp(log_tot_in_scaling_fct),
    log_C_contact_scale = log_C_contact_scale,
    log_tot_in_scaling_fct = log_tot_in_scaling_fct
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

#' Vectorised back-transform: spline-coefficient matrix -> age-level parameters.
#'
#' Each row of samps is c(alpha, gamma) (2*SPLINE_K coefficients). The N_AGE
#' age values are reconstructed via the basis and exponentiated:
#'   C_contact_scale    = exp(alpha %*% t(SPLINE_B))
#'   tot_in_scaling_fct = exp(gamma %*% t(SPLINE_B))
theta_to_orig <- function(samps) {
  if (is.null(dim(samps))) {
    samps <- matrix(samps, nrow = 1L)
  }
  alpha <- samps[, 1:SPLINE_K, drop = FALSE]
  gamma <- samps[, (SPLINE_K + 1):(2 * SPLINE_K), drop = FALSE]

  contact_cols <- exp(alpha %*% t(SPLINE_B))
  tot_in_cols  <- exp(gamma %*% t(SPLINE_B))
  colnames(contact_cols) <- paste0("C_contact_scale_", seq_len(N_AGE))
  colnames(tot_in_cols)  <- paste0("tot_in_scaling_fct_", seq_len(N_AGE))

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

# ── 3a. Spline basis (dimension reduction) and log-prior ─────────────────────
#
# Instead of estimating N_AGE free values for each age-varying log vector, each
# vector is represented by SPLINE_K B-spline coefficients:
#   log_C_contact_scale    = SPLINE_B %*% alpha   (alpha = theta[1:K])
#   log_tot_in_scaling_fct = SPLINE_B %*% gamma   (gamma = theta[K + 1:K])
# so theta has length 2*SPLINE_K instead of 2*N_AGE.
#
# Per-age prior medians (unchanged targets, used to centre the coefficients).
CONTACT_PRIOR_MEANS <- log(c(0.2506, 0.8799, 0.6041, 2.2305, 6.2280, 102.0037))
TOT_IN_PRIOR_MEANS <- log(c(0.96, 0.23, 0.06, 0.33, 1.4, 2.0))

SPLINE_N_KNOTS   <- 1L     # internal B-spline knots per input vector
SPLINE_RW_SD     <- 0.11   # tau: RW2 scale for deviations from the prior curve
SPLINE_ANCHOR_SD <- 0.375   # sd on the first two deviation coefficients

# Cubic B-spline basis over the age index, using exactly three internal knots.
SPLINE_AGE_GRID <- seq_len(N_AGE)
SPLINE_INTERNAL_KNOTS <- as.numeric(stats::quantile(
  SPLINE_AGE_GRID,
  probs = seq(0, 1, length.out = SPLINE_N_KNOTS + 2L)[-c(1L, SPLINE_N_KNOTS + 2L)],
  names = FALSE,
  type = 7
))
SPLINE_B <- splines::bs(
  SPLINE_AGE_GRID,
  knots = SPLINE_INTERNAL_KNOTS,
  degree = 3L,
  Boundary.knots = range(SPLINE_AGE_GRID),
  intercept = TRUE
)
SPLINE_K <- ncol(SPLINE_B)

# Parameter names are rebuilt from SPLINE_K so downstream matrices remain
# consistent when the number of internal knots changes.
param_names_log <- c(
  paste0("alpha_", seq_len(SPLINE_K)),
  paste0("gamma_", seq_len(SPLINE_K))
)

# Least-squares projection of the per-age prior means onto the basis. The full
# projected coefficient vector is the reference curve for the centered RW2 prior
# and can also be used for initialization.
project_to_spline <- function(y) {
  as.numeric(solve(crossprod(SPLINE_B), crossprod(SPLINE_B, y)))
}
CONTACT_COEF_PRIOR_MEANS <- project_to_spline(CONTACT_PRIOR_MEANS)
TOT_IN_COEF_PRIOR_MEANS  <- project_to_spline(TOT_IN_PRIOR_MEANS)
CONTACT_COEF_PRIOR_SDS   <- rep(SPLINE_ANCHOR_SD, SPLINE_K)
TOT_IN_COEF_PRIOR_SDS    <- rep(SPLINE_ANCHOR_SD, SPLINE_K)

# Update spline-coefficient prior centers and scales. HMC uses this to turn the
# Nelder-Mead MAP coefficients into valid priors without changing constraints:
# theta remains unconstrained and exp(SPLINE_B %*% theta) enforces positivity.
configure_spline_priors <- function(contact_mean = CONTACT_COEF_PRIOR_MEANS,
                                    tot_in_mean = TOT_IN_COEF_PRIOR_MEANS,
                                    contact_sd = SPLINE_ANCHOR_SD,
                                    tot_in_sd = SPLINE_ANCHOR_SD,
                                    rw_sd = SPLINE_RW_SD) {
  stopifnot(length(contact_mean) == SPLINE_K, length(tot_in_mean) == SPLINE_K)
  if (length(contact_sd) == 1L) contact_sd <- rep(contact_sd, SPLINE_K)
  if (length(tot_in_sd) == 1L) tot_in_sd <- rep(tot_in_sd, SPLINE_K)
  stopifnot(length(contact_sd) == SPLINE_K, length(tot_in_sd) == SPLINE_K)
  if (any(!is.finite(contact_mean)) || any(!is.finite(tot_in_mean)) ||
      any(!is.finite(contact_sd)) || any(!is.finite(tot_in_sd)) ||
      any(contact_sd <= 0) || any(tot_in_sd <= 0) ||
      !is.finite(rw_sd) || rw_sd <= 0) {
    stop("Spline prior means must be finite; standard deviations must be positive and finite")
  }

  CONTACT_COEF_PRIOR_MEANS <<- as.numeric(contact_mean)
  TOT_IN_COEF_PRIOR_MEANS  <<- as.numeric(tot_in_mean)
  CONTACT_COEF_PRIOR_SDS   <<- as.numeric(contact_sd)
  TOT_IN_COEF_PRIOR_SDS    <<- as.numeric(tot_in_sd)
  SPLINE_RW_SD             <<- as.numeric(rw_sd)
  invisible(TRUE)
}

# Centered P-spline log-prior for one coefficient vector.
#
# Let mean0 denote the complete projected reference coefficient curve and define
# delta = coef - mean0. The first two deviations anchor the overall level and
# slope relative to that curve. Subsequent deviations follow an RW2 prior:
#   delta[k] - 2*delta[k-1] + delta[k-2] ~ Normal(0, SPLINE_RW_SD).
# This preserves the curvature of the full reference curve in the prior mean,
# while allowing data-driven departures that vary smoothly across coefficients.
rw2_log_prior <- function(coef, mean0, sd0) {
  K <- length(coef)
  if (length(mean0) != K || length(sd0) != K) {
    stop("coef, mean0 and sd0 must have the same length")
  }

  delta <- coef - mean0

  lp <- dnorm(delta[1], mean = 0, sd = sd0[1], log = TRUE) +
        dnorm(delta[2], mean = 0, sd = sd0[2], log = TRUE)

  if (K >= 3L) {
    for (k in 3:K) {
      lp <- lp + dnorm(
        delta[k],
        mean = 2 * delta[k - 1L] - delta[k - 2L],
        sd = SPLINE_RW_SD,
        log = TRUE
      )
    }
  }
  lp
}

log_prior <- function(theta) {
  alpha <- theta[1:SPLINE_K]
  gamma <- theta[(SPLINE_K + 1):(2 * SPLINE_K)]
  rw2_log_prior(alpha, CONTACT_COEF_PRIOR_MEANS, CONTACT_COEF_PRIOR_SDS) +
  rw2_log_prior(gamma, TOT_IN_COEF_PRIOR_MEANS, TOT_IN_COEF_PRIOR_SDS)
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

# ── 3c. Prevalence log-likelihood (binomial) ─────────────────────────────────
#
# case_a ~ Binomial(N_a, prev_a), a = age group index
#   case_a  = observed HCV positives in age group a  (obs_pos)
#   N_a     = sample size / number tested in age group a  (obs_tot)
#   prev_a  = model-implied prevalence q_age[a]
#
compute_prevalence_loglik <- function(q_age, obs_pos, obs_tot) {
  p <- pmin(pmax(q_age, 1e-12), 1 - 1e-12)
  sum(dbinom(obs_pos, size = obs_tot, prob = p, log = TRUE))
}

# ── 3d. Analytical gradient of log-prior ─────────────────────────────────────
#
# Gradients of direct log-scale priors.
#
# Note: not called during sampling (numerical gradients are used); kept for
# reference and unit-testing against finite differences.
#
grad_log_prior_analytical <- function(theta) {
  grad_coef <- function(coef, mean0, sd0) {
    K     <- length(coef)
    delta <- coef - mean0
    g     <- numeric(K)

    g[1] <- -delta[1] / sd0[1]^2
    g[2] <- -delta[2] / sd0[2]^2

    if (K >= 3L) {
      for (k in 3:K) {
        d <- (delta[k] - 2 * delta[k - 1L] + delta[k - 2L]) /
             SPLINE_RW_SD^2
        g[k]      <- g[k]      - d
        g[k - 1L] <- g[k - 1L] + 2 * d
        g[k - 2L] <- g[k - 2L] - d
      }
    }
    g
  }
  c(
    grad_coef(theta[1:SPLINE_K], CONTACT_COEF_PRIOR_MEANS, CONTACT_COEF_PRIOR_SDS),
    grad_coef(theta[(SPLINE_K + 1):(2 * SPLINE_K)], TOT_IN_COEF_PRIOR_MEANS, TOT_IN_COEF_PRIOR_SDS)
  )
}

# ── 3e. Log-likelihood ───────────────────────────────────────────────────────
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

  sigma_pop <- if (!is.null(data$sigma_pop)) {
    data$sigma_pop
  } else {
    rep(0.05, length(obs_tot))
  }

  ll_pop  <- compute_population_composition_loglik(pi_model, q_obs, sigma_pop)
  ll_prev <- compute_prevalence_loglik(obs$q_age, obs_pos, obs_tot)

  if (!is.finite(ll_pop) || !is.finite(ll_prev)) return(-Inf)

  ll_pop + ll_prev
}

# ── 3f. Log-posterior ────────────────────────────────────────────────────────
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

# ── 3g. Numerical gradient (central finite differences) ───────────────────────
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

  if (length(theta_init) != 2L * SPLINE_K || any(!is.finite(theta_init))) {
    stop(sprintf("theta_init must contain %d finite values", 2L * SPLINE_K))
  }

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
