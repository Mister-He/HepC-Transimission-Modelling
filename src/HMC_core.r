# =============================================================================
# HMC_core.r  —  HMC Calibration for the HCV PWID Compartmental Model
#
# Append this file after setup.R (which sources sim.cpp and defines params,
# data, idx(), y0, etc.).
#
# Model architecture:
#   1. Twelve positive age-specific parameters modify the C++ ODE model:
#        C_contact_scale[k] = exp(theta[k])
#        inflow_scale[k]    = exp(theta[N_AGE + k])
#   2. Each contact parameter scales one contact-matrix row and each inflow
#      parameter scales the baseline total inflow for that age group.
#   3. The final ODE state is compared with observed age composition and
#      age-specific HCV prevalence. No additional age-shape penalty is used.
#   4. HMC samples the direct log parameters using finite-difference gradients
#      through the black-box C++ ODE solver.
#
# Parameters estimated (2*N_AGE = 12 for six age groups):
#   theta[1:N_AGE]                 = log contact-row scaling factors
#   theta[(N_AGE+1):(2*N_AGE)]     = log age-specific inflow scaling factors
#
# Scaling relationships:
#   c_true[k]          = c_composite[k] / inflow_scale[k]
#   lambda1[k]         = lambda3[k] * c_true[k]
#   Each relationship is updated before every simulation call.
#
# Joint observation model (J-stratum only, per age group i = 1..10):
#   ALR(q_obs)[a] ~ Normal(ALR(pi_model)[a], sigma_pop[a])  population composition
#   obs_pos[a]    ~ Binomial(obs_tot[a], q_a(theta))         HCV prevalence
#   obs_pos[a] = round(obs_prev[a] * obs_tot[a]) HCV positives per age group.
#
# Prior: independent normal priors on each direct log-scale parameter.
#
# Gradient: central finite differences through the C++ ODE solver
#   N_PARAMS = 12, so one gradient evaluation normally requires
#   2 × N_PARAMS = 24 ODE runs.
#   One HMC proposal with L leapfrog steps costs approximately
#   (L + 1) × 24 ODE runs, plus one ODE run for the proposed log posterior.
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

#' Unconstrained direct log parameters -> constrained named list.
constrain_theta <- function(theta) {
  if (length(theta) != 2L * N_AGE || any(!is.finite(theta))) {
    stop(sprintf("theta must contain %d finite values", 2L * N_AGE))
  }
  log_C_contact_scale <- theta[seq_len(N_AGE)]
  log_inflow_scale <- theta[N_AGE + seq_len(N_AGE)]

  list(
    C_contact_scale = exp(log_C_contact_scale),
    inflow_scale = exp(log_inflow_scale),
    log_C_contact_scale = log_C_contact_scale,
    log_inflow_scale = log_inflow_scale
  )
}

#' Build a full parameter list from unconstrained theta.
#'
#' C_contact_scale[k] scales the contact matrix row k and
#' inflow_scale[k] multiplies baseline beta[k] and determines
#' c_true[k] = c_composite[k] / inflow_scale[k].
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

  if (length(base_params$beta) != N_AGE || any(!is.finite(base_params$beta)) ||
      any(base_params$beta <= 0)) {
    stop("base_params$beta must contain one positive finite baseline inflow per age group")
  }

  # Scale each baseline age-specific total inflow.
  pm$beta <- base_params$beta * p$inflow_scale

  c_true    <- base_params$c_composite / p$inflow_scale
  pm$lambda1 <- base_params$lambda3 * c_true
  pm
}

#' Vectorised back-transform from direct log parameters to age-level values.
theta_to_orig <- function(samps) {
  if (is.null(dim(samps))) {
    samps <- matrix(samps, nrow = 1L)
  }
  if (ncol(samps) != 2L * N_AGE) stop(sprintf("samps must have %d columns", 2L * N_AGE))
  contact_cols <- exp(samps[, seq_len(N_AGE), drop = FALSE])
  tot_in_cols  <- exp(samps[, N_AGE + seq_len(N_AGE), drop = FALSE])
  colnames(contact_cols) <- paste0("C_contact_scale_", seq_len(N_AGE))
  colnames(tot_in_cols)  <- paste0("inflow_scale_", seq_len(N_AGE))

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

# ── 3a. Direct age-specific log priors ───────────────────────────────────────
CONTACT_PRIOR_MEANS <- log(c(3.68983, 0.06629, 0.08773, 0.61861, 1.74367,1.24803))
INFLOW_SCALE_PRIOR_MEANS <- log(c(0.63544, 4.33889, 0.05326, 0.33554, 1.16247, 1.55939))
CONTACT_PRIOR_SDS <- rep(0.375, N_AGE)
INFLOW_SCALE_PRIOR_SDS <- rep(0.375, N_AGE)

param_names_log <- c(
  paste0("log_C_contact_scale_", seq_len(N_AGE)),
  paste0("log_inflow_scale_", seq_len(N_AGE))
)

configure_age_priors <- function(contact_mean = CONTACT_PRIOR_MEANS,
                                 inflow_scale_mean = INFLOW_SCALE_PRIOR_MEANS,
                                 contact_sd = CONTACT_PRIOR_SDS,
                                 inflow_scale_sd = INFLOW_SCALE_PRIOR_SDS) {
  if (length(contact_sd) == 1L) contact_sd <- rep(contact_sd, N_AGE)
  if (length(inflow_scale_sd) == 1L) inflow_scale_sd <- rep(inflow_scale_sd, N_AGE)
  stopifnot(length(contact_mean) == N_AGE, length(inflow_scale_mean) == N_AGE,
            length(contact_sd) == N_AGE, length(inflow_scale_sd) == N_AGE)
  if (any(!is.finite(contact_mean)) || any(!is.finite(inflow_scale_mean)) ||
      any(!is.finite(contact_sd)) || any(!is.finite(inflow_scale_sd)) ||
      any(contact_sd <= 0) || any(inflow_scale_sd <= 0)) {
    stop("Prior means must be finite; standard deviations must be positive and finite")
  }
  CONTACT_PRIOR_MEANS <<- as.numeric(contact_mean)
  INFLOW_SCALE_PRIOR_MEANS <<- as.numeric(inflow_scale_mean)
  CONTACT_PRIOR_SDS <<- as.numeric(contact_sd)
  INFLOW_SCALE_PRIOR_SDS <<- as.numeric(inflow_scale_sd)
  invisible(TRUE)
}

log_prior <- function(theta) {
  if (length(theta) != 2L * N_AGE || any(!is.finite(theta))) return(-Inf)
  sum(dnorm(theta[seq_len(N_AGE)], CONTACT_PRIOR_MEANS, CONTACT_PRIOR_SDS, log = TRUE)) +
    sum(dnorm(theta[N_AGE + seq_len(N_AGE)], INFLOW_SCALE_PRIOR_MEANS,
              INFLOW_SCALE_PRIOR_SDS, log = TRUE))
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
  c(
    -(theta[seq_len(N_AGE)] - CONTACT_PRIOR_MEANS) / CONTACT_PRIOR_SDS^2,
    -(theta[N_AGE + seq_len(N_AGE)] - INFLOW_SCALE_PRIOR_MEANS) /
      INFLOW_SCALE_PRIOR_SDS^2
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

  if (length(theta_init) != 2L * N_AGE || any(!is.finite(theta_init))) {
    stop(sprintf("theta_init must contain %d finite values", 2L * N_AGE))
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
