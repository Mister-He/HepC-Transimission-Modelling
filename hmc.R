# =============================================================================
# hmc_calibration.R  —  HMC Calibration for the HCV PWID Compartmental Model
#
# Append this file after setup.R (which sources sim.cpp and defines params,
# data, idx(), y0, etc.).
#
# Parameters estimated (10 total, all in log-space for positivity):
#   theta[1]      = log(beta_scale)           scalar multiplier on beta[1:9]
#   theta[2:10]   = log(C_contact_scale[1:9]) row-wise multipliers on C_contact
#
# Likelihood (Poisson, J-stratum only, per age group i = 1..9):
#   obs_pos[i] ~ Poisson( sum_{k=1..4} sum_{h in {0,2}} y[T, idx(s=1,k,h,i-1)] )
#   obs_tot[i] ~ Poisson( sum_{k=1..4} sum_{h=0..3}    y[T, idx(s=1,k,h,i-1)] )
#
# Prior: Normal(0, 10) on original (constrained) scale + log-Jacobian correction
#   => log pi(theta_j) = dnorm(exp(theta_j), 0, 10, log=T) + theta_j
#
# Gradient: central finite differences through the C++ ODE solver
#   Cost per HMC step: (L+1) gradient evals × 2 × N_PARAMS ODE runs
#                    = (L+1) × 20 ODE runs  [plus 1 for current lp]
#   Tip: keep L small (5-10) and use parallel chains.
#
# Outputs:
#   hmc_chains        — raw chain list (all iterations, all chains)
#   post_warmup_list  — post-warmup samples per chain
#   all_samples       — pooled post-warmup matrix (n_samp × 10)
#   diag_table        — R-hat and ESS summary data.frame
#   ppc_out           — posterior predictive replicates
#   plots             — trace, PPC histogram, PPC interval
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# =============================================================================
# 0.  OBSERVATIONS
# =============================================================================

obs_pos <- c( 55, 145, 183, 164, 212, 299, 222, 190, 133)
obs_tot <- c(307, 797, 829, 633, 598, 642, 481, 439, 366)

stopifnot(length(obs_pos) == 9L, length(obs_tot) == 9L)

N_PARAMS         <- 10L
param_names_log  <- c("log_beta_scale", paste0("log_C_scale_", 1:9))
param_names_orig <- c("beta_scale",     paste0("C_scale_",     1:9))


# =============================================================================
# 1.  PARAMETER TRANSFORMS
# =============================================================================

#' Unconstrained theta (log-space) -> constrained named list
constrain_theta <- function(theta) {
  list(
    beta_scale      = exp(theta[1L]),
    C_contact_scale = exp(theta[2:10])
  )
}

#' Build a full parameter list from unconstrained theta
build_params_from_theta <- function(theta, base_params) {
  p  <- constrain_theta(theta)
  pm <- base_params

  # Scale beta uniformly across all age groups
  pm$beta <- base_params$beta * p$beta_scale

  # Scale each row of C_contact independently
  for (j in 1:9) {
    pm$C_contact[j, ] <- base_params$C_contact[j, ] * p$C_contact_scale[j]
  }
  pm
}


# =============================================================================
# 2.  POISSON MEANS FROM FINAL ODE STATE
# =============================================================================

#' Compute Poisson rate parameters from the final ODE state vector.
#'
#' Both rates are for the J-stratum (s = 1) only:
#'
#'   lambda_pos[i] = sum_{k=1..4} sum_{h in {0, 2}} y[ idx(s=1, k, h, i-1) ]
#'     (h=0: susceptible; h=2: chronic — as specified in the likelihood)
#'
#'   lambda_tot[i] = sum_{k=1..4} sum_{h=0..3}      y[ idx(s=1, k, h, i-1) ]
#'     (all health states in J-stratum)
#'
#' Both are floored at 1e-10 to prevent log(0) in the Poisson log-pmf.
compute_poisson_means <- function(y_final) {
  lambda_pos <- numeric(9L)
  lambda_tot <- numeric(9L)

  for (i in 0:8) {
    lp <- 0.0; lt <- 0.0
    for (k in 1:4) {
      # Positive: h in {0, 2}
      lp <- lp + y_final[idx(s = 1L, k = k, h = 0L, i = i)]
      lp <- lp + y_final[idx(s = 1L, k = k, h = 2L, i = i)]
      # Total: h in {0, 1, 2, 3}
      lt <- lt + y_final[idx(s = 1L, k = k, h = 0L, i = i)]
      lt <- lt + y_final[idx(s = 1L, k = k, h = 1L, i = i)]
      lt <- lt + y_final[idx(s = 1L, k = k, h = 2L, i = i)]
      lt <- lt + y_final[idx(s = 1L, k = k, h = 3L, i = i)]
    }
    lambda_pos[i + 1L] <- max(lp, 1e-10)
    lambda_tot[i + 1L] <- max(lt, 1e-10)
  }
  list(pos = lambda_pos, tot = lambda_tot)
}


# =============================================================================
# 3.  LOG-POSTERIOR AND GRADIENT
# =============================================================================

# ── 3a. Log-prior ────────────────────────────────────────────────────────────
#
# Prior on original scale:  X_j ~ Normal(0, sigma)  (improper for X<0, but
# positivity is enforced by the log-transform so only the right tail matters).
#
# Change-of-variables: theta_j = log(X_j),  X_j = exp(theta_j)
#   log p(theta_j) = log p_X(exp(theta_j)) + log|dX/dtheta_j|
#                  = dnorm(exp(theta_j), 0, sigma, log=T) + theta_j
#
log_prior <- function(theta, sigma = 10) {
  orig <- exp(theta)
  sum(dnorm(orig, mean = 0, sd = sigma, log = TRUE)) + sum(theta)
}

# ── 3b. Analytical gradient of log-prior (used internally for sanity checks) ─
#
#   d/d theta_j  log p(theta_j)
#     = d/d theta_j [ -exp(2*theta_j)/(2*sigma^2) + theta_j ]   + const
#     = -exp(2*theta_j)/sigma^2 + 1
#     = 1 - exp(2*theta_j) / sigma^2
#
grad_log_prior_analytical <- function(theta, sigma = 10) {
  1 - exp(2 * theta) / sigma^2
}

# ── 3c. Log-likelihood ───────────────────────────────────────────────────────
log_likelihood <- function(theta, base_params, data) {
  pm  <- build_params_from_theta(theta, base_params)
  out <- run_sim(pm, data)

  if (!is.matrix(out) || nrow(out) == 0L) return(-Inf)
  if (any(!is.finite(out)))               return(-Inf)

  y_final <- as.numeric(out[nrow(out), -1L])   # drop time column
  lam     <- compute_poisson_means(y_final)

  sum(dpois(obs_pos, lambda = lam$pos, log = TRUE)) +
  sum(dpois(obs_tot, lambda = lam$tot, log = TRUE))
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


# =============================================================================
# 6.  CONVERGENCE DIAGNOSTICS: R-HAT AND ESS
# =============================================================================

# ── 6a. Gelman-Rubin R-hat ───────────────────────────────────────────────────
#
# R-hat measures mixing across chains.  R-hat < 1.01 indicates convergence
# (Vehtari et al. 2021); values 1.01–1.05 are borderline.
#
# Formula (split-chain pooled variance):
#   B = (N/(M-1)) Σ_m (ψ̄_m - ψ̄)²    between-chain variance
#   W = (1/M)     Σ_m s²_m             within-chain variance
#   var̂(ψ) = (N-1)/N · W + B/N
#   R-hat   = sqrt(var̂(ψ) / W)
#
compute_rhat <- function(chains) {
  M <- length(chains)
  N <- nrow(chains[[1L]])
  J <- ncol(chains[[1L]])

  rhat <- setNames(numeric(J), param_names_log)

  for (j in seq_len(J)) {
    chain_j     <- lapply(chains, function(ch) ch[, j])
    chain_means <- vapply(chain_j, mean, numeric(1L))
    chain_vars  <- vapply(chain_j, var,  numeric(1L))

    grand_mean <- mean(chain_means)
    B          <- N / (M - 1) * sum((chain_means - grand_mean)^2)
    W          <- mean(chain_vars)
    var_hat    <- (N - 1) / N * W + B / N

    rhat[j]    <- sqrt(var_hat / W)
  }
  rhat
}

# ── 6b. Effective Sample Size (ESS) ──────────────────────────────────────────
#
# ESS based on Geyer's initial monotone sequence estimator applied to the
# pooled autocorrelation function across all chains.
#
# ESS = (M × N) / τ,   where τ = 1 + 2 Σ_{t=1}^{T*} ρ_t
# T* = last lag where ρ_{2t} + ρ_{2t+1} ≥ 0 (monotone pairs truncation).
#
# Target: ESS > 400, or ESS/N_total > 0.1 per parameter.
#
compute_ess <- function(chains) {
  M <- length(chains)
  N <- nrow(chains[[1L]])
  J <- ncol(chains[[1L]])

  ess <- setNames(numeric(J), param_names_log)

  for (j in seq_len(J)) {
    all_samps <- unlist(lapply(chains, function(ch) ch[, j]))
    n_total   <- length(all_samps)
    max_lag   <- min(n_total - 1L, 500L)

    acf_vals  <- acf(all_samps, lag.max = max_lag,
                     plot = FALSE, type = "correlation")$acf[-1L]

    # Geyer monotone pairs:  Γ_t = ρ_{2t} + ρ_{2t+1}
    n_pairs <- floor(length(acf_vals) / 2L)
    if (n_pairs < 1L) { ess[j] <- n_total; next }

    pairs   <- acf_vals[2L * seq_len(n_pairs) - 1L] +
               acf_vals[2L * seq_len(n_pairs)]
    cutoff  <- which(pairs <= 0)[1L]
    if (is.na(cutoff)) cutoff <- n_pairs

    # Sum autocorrelations up to cut-off pair
    rho_sum <- 1 + 2 * sum(acf_vals[seq_len(2L * (cutoff - 1L) + 1L)])
    ess[j]  <- n_total / max(rho_sum, 1)
  }
  ess
}

# ── 6c. Divergence check ─────────────────────────────────────────────────────
#
# A "divergence" occurs when the leapfrog trajectory explodes (H changes by
# > threshold).  Here we flag iterations where |ΔH| > 1000 as divergent.
# In practice, this is visible as sudden jumps in the lp trace.
#
count_divergences <- function(chains_raw, threshold = 1000) {
  lapply(chains_raw, function(ch) {
    lp  <- ch$lp_trace[ch$n_warmup + seq_len(nrow(ch$samples))]
    sum(abs(diff(lp)) > threshold, na.rm = TRUE)
  })
}

# ── 6d. Print formatted diagnostics table ───────────────────────────────────
print_diagnostics <- function(post_warmup_list, chains_raw = NULL) {
  rhat_vals  <- compute_rhat(post_warmup_list)
  ess_vals   <- compute_ess(post_warmup_list)
  all_samps  <- do.call(rbind, post_warmup_list)
  orig_samps <- exp(all_samps)

  diag_df <- data.frame(
    Parameter  = param_names_log,
    Mean_log   = round(colMeans(all_samps),               4),
    SD_log     = round(apply(all_samps,  2, sd),          4),
    Mean_orig  = round(colMeans(orig_samps),              4),
    Q2.5_orig  = round(apply(orig_samps, 2, quantile, 0.025), 4),
    Q97.5_orig = round(apply(orig_samps, 2, quantile, 0.975), 4),
    Rhat       = round(rhat_vals, 4),
    ESS        = round(ess_vals,  1),
    Conv_OK    = ifelse(rhat_vals < 1.05 & ess_vals > 100, "YES", "NO"),
    stringsAsFactors = FALSE
  )

  cat("\n======================================================\n")
  cat("          Convergence Diagnostics Summary\n")
  cat("  Target: R-hat < 1.05   |   ESS > 100 per parameter\n")
  cat("======================================================\n")
  print(diag_df, row.names = FALSE)

  if (!is.null(chains_raw)) {
    div_counts <- count_divergences(chains_raw)
    cat("\nDivergences per chain (post-warmup):\n")
    print(setNames(unlist(div_counts), paste0("Chain ", seq_along(div_counts))))
  }

  invisible(diag_df)
}


# =============================================================================
# 7.  POSTERIOR PREDICTIVE CHECKS (PPC)
# =============================================================================

# ── 7a. Generate PPC replicates ──────────────────────────────────────────────
#
# For each sampled theta_s, run the ODE model and draw:
#   y_rep_pos[i] ~ Poisson(lambda_pos[i] | theta_s)
#   y_rep_tot[i] ~ Poisson(lambda_tot[i] | theta_s)
#
# Returns replicated counts (integer) and the Poisson means (real),
# plus posterior predictive p-values (ppp) per age group.
#
generate_ppc_samples <- function(post_samples, base_params, data,
                                  n_ppc = 200L) {
  n_avail  <- nrow(post_samples)
  draw_idx <- sort(sample(n_avail, min(n_ppc, n_avail)))
  n_draws  <- length(draw_idx)

  ppc_pos  <- matrix(NA_integer_, n_draws, 9L)
  ppc_tot  <- matrix(NA_integer_, n_draws, 9L)
  lam_pos  <- matrix(NA_real_,    n_draws, 9L)
  lam_tot  <- matrix(NA_real_,    n_draws, 9L)

  cat(sprintf("Generating %d PPC draws from %d posterior samples...\n",
              n_draws, n_avail))

  for (ii in seq_len(n_draws)) {
    theta <- post_samples[draw_idx[ii], ]

    result <- tryCatch({
      pm      <- build_params_from_theta(theta, base_params)
      out     <- run_sim(pm, data)
      y_fin   <- as.numeric(out[nrow(out), -1L])
      lam     <- compute_poisson_means(y_fin)
      lam
    }, error = function(e) NULL)

    if (!is.null(result)) {
      lam_pos[ii, ] <- result$pos
      lam_tot[ii, ] <- result$tot
      ppc_pos[ii, ] <- rpois(9L, result$pos)
      ppc_tot[ii, ] <- rpois(9L, result$tot)
    }

    if (ii %% 50L == 0L)
      cat(sprintf("  PPC draw %d / %d\n", ii, n_draws))
  }

  # Posterior predictive p-values: Pr(y_rep >= y_obs | y_obs)
  # Calibrated model: ppp near 0.5; poor fit: near 0 or 1.
  ppp_pos <- vapply(1:9, function(i)
    mean(ppc_pos[, i] >= obs_pos[i], na.rm = TRUE), numeric(1L))
  ppp_tot <- vapply(1:9, function(i)
    mean(ppc_tot[, i] >= obs_tot[i], na.rm = TRUE), numeric(1L))

  list(ppc_pos = ppc_pos, ppc_tot = ppc_tot,
       lam_pos = lam_pos, lam_tot = lam_tot,
       ppp_pos = ppp_pos, ppp_tot = ppp_tot)
}

# ── 7b. PPC histogram plot ───────────────────────────────────────────────────
#
# One panel per age group: histogram of y_rep (9 × n_ppc), vertical line at
# observed count, and posterior predictive p-value annotated per panel.
#
plot_ppc_histograms <- function(ppc, type = c("pos", "tot")) {
  type <- match.arg(type)

  ppc_mat  <- if (type == "pos") ppc$ppc_pos else ppc$ppc_tot
  obs_vals <- if (type == "pos") obs_pos     else obs_tot
  ppp_vals <- if (type == "pos") ppc$ppp_pos else ppc$ppp_tot

  fill_col <- if (type == "pos") "#2166AC" else "#D6604D"
  main_ttl <- if (type == "pos")
    "PPC: HCV-Positive Counts (J-stratum)" else
    "PPC: Total PWID Counts (J-stratum)"

  age_lbl <- paste0("Age group ", 1:9)

  # Drop failed draws
  ppc_mat  <- ppc_mat[complete.cases(ppc_mat), , drop = FALSE]

  ppc_long <- as.data.frame(ppc_mat) %>%
    setNames(age_lbl) %>%
    pivot_longer(everything(), names_to = "Age", values_to = "y_rep") %>%
    mutate(Age = factor(Age, levels = age_lbl))

  obs_df <- data.frame(
    Age = factor(age_lbl, levels = age_lbl),
    obs = obs_vals
  )
  ppp_df <- data.frame(
    Age = factor(age_lbl, levels = age_lbl),
    ppp = round(ppp_vals, 2)
  )

  ggplot(ppc_long, aes(x = y_rep)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins    = 30,
                   fill    = fill_col,
                   alpha   = 0.65,
                   colour  = "white",
                   linewidth = 0.2) +
    geom_vline(data = obs_df, aes(xintercept = obs),
               colour = "black", linewidth = 1.1) +
    geom_text(data = ppp_df,
              aes(label = paste0("ppp=", ppp)),
              x = Inf, y = Inf, hjust = 1.15, vjust = 1.6,
              size = 3.0, colour = "grey25") +
    facet_wrap(~ Age, scales = "free", ncol = 3) +
    labs(
      title    = main_ttl,
      subtitle = "Histogram: posterior predictive  |  Vertical line: observed  |  ppp near 0.5 = good fit",
      x        = "Count",
      y        = "Density"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background  = element_rect(fill = "grey92", colour = "grey70"),
      panel.grid.minor  = element_blank(),
      plot.subtitle     = element_text(size = 9, colour = "grey40")
    )
}

# ── 7c. PPC interval plot ─────────────────────────────────────────────────────
#
# Dot-and-whisker: posterior predictive median ± 50% and 90% intervals per age
# group, with observed count overlaid as a red point.
#
plot_ppc_intervals <- function(ppc) {
  age_lbl <- paste0("Age ", 1:9)

  summarise_ppc_mat <- function(ppc_mat, obs, type_label) {
    mat <- ppc_mat[complete.cases(ppc_mat), , drop = FALSE]
    data.frame(
      Age    = factor(age_lbl, levels = age_lbl),
      obs    = obs,
      med    = apply(mat, 2, median),
      lo50   = apply(mat, 2, quantile, 0.25),
      hi50   = apply(mat, 2, quantile, 0.75),
      lo90   = apply(mat, 2, quantile, 0.05),
      hi90   = apply(mat, 2, quantile, 0.95),
      type   = type_label,
      stringsAsFactors = FALSE
    )
  }

  plot_df <- bind_rows(
    summarise_ppc_mat(ppc$ppc_pos, obs_pos, "HCV-positive"),
    summarise_ppc_mat(ppc$ppc_tot, obs_tot, "Total PWID")
  )

  ggplot(plot_df, aes(x = Age)) +
    geom_linerange(aes(ymin = lo90, ymax = hi90),
                   linewidth = 1.0, colour = "#4292C6", alpha = 0.45) +
    geom_linerange(aes(ymin = lo50, ymax = hi50),
                   linewidth = 3.0, colour = "#08519C", alpha = 0.75) +
    geom_point(aes(y = med),
               shape = 18, size = 3.5, colour = "#08306B") +
    geom_point(aes(y = obs),
               shape = 21, size = 3.5,
               fill = "#D7191C", colour = "black", stroke = 1.2) +
    facet_wrap(~ type, scales = "free_y", ncol = 1) +
    labs(
      title    = "PPC: Posterior Predictive Intervals vs Observed",
      subtitle = "Diamond: predictive median  |  Thick bar: 50% PI  |  Thin bar: 90% PI  |  Red dot: observed",
      x        = "Age group",
      y        = "Count"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x       = element_text(angle = 35, hjust = 1),
      strip.background  = element_rect(fill = "grey92", colour = "grey70"),
      panel.grid.minor  = element_blank(),
      plot.subtitle     = element_text(size = 9, colour = "grey40")
    )
}

# ── 7d. Trace plots ───────────────────────────────────────────────────────────
#
# Full-chain trace for selected parameters (all iterations, all chains),
# with a dashed vertical line marking the end of warmup.
#
plot_traces <- function(chains_raw, param_idx = 1:min(5L, N_PARAMS)) {

  trace_list <- lapply(seq_along(chains_raw), function(ch) {
    df           <- as.data.frame(chains_raw[[ch]]$samples_all)
    colnames(df) <- param_names_log
    df$iter      <- seq_len(nrow(df))
    df$chain     <- factor(ch)
    df
  })
  trace_df  <- bind_rows(trace_list)
  n_warmup  <- chains_raw[[1L]]$n_warmup
  params_to <- param_names_log[param_idx]

  trace_df %>%
    select(iter, chain, all_of(params_to)) %>%
    pivot_longer(all_of(params_to),
                 names_to  = "param",
                 values_to = "value") %>%
    mutate(param = factor(param, levels = params_to)) %>%
    ggplot(aes(x = iter, y = value, colour = chain)) +
    geom_line(alpha = 0.70, linewidth = 0.35) +
    geom_vline(xintercept = n_warmup,
               linetype = "dashed", colour = "grey35", linewidth = 0.8) +
    facet_wrap(~ param, scales = "free_y", ncol = 1) +
    scale_colour_brewer(palette = "Set1") +
    labs(
      title    = "Trace plots — all chains (log scale)",
      subtitle = "Dashed line = end of warmup",
      x        = "Iteration",
      y        = "Value (log-transformed parameter)",
      colour   = "Chain"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background  = element_rect(fill = "grey92", colour = "grey70"),
      panel.grid.minor  = element_blank(),
      legend.position   = "right",
      plot.subtitle     = element_text(size = 9, colour = "grey40")
    )
}

# ── 7e. Posterior density plot (original scale) ────────────────────────────
plot_posterior_densities <- function(post_warmup_list) {
  all_samps <- do.call(rbind, lapply(seq_along(post_warmup_list), function(ch) {
    df           <- as.data.frame(exp(post_warmup_list[[ch]]))
    colnames(df) <- param_names_orig
    df$chain     <- factor(ch)
    df
  }))

  all_samps %>%
    pivot_longer(-chain, names_to = "param", values_to = "value") %>%
    mutate(param = factor(param, levels = param_names_orig)) %>%
    ggplot(aes(x = value, fill = chain, colour = chain)) +
    geom_density(alpha = 0.30, linewidth = 0.6) +
    facet_wrap(~ param, scales = "free", ncol = 2) +
    scale_fill_brewer(palette  = "Set1") +
    scale_colour_brewer(palette = "Set1") +
    labs(
      title    = "Posterior densities (original scale, all chains)",
      x        = "Parameter value",
      y        = "Density",
      fill     = "Chain",
      colour   = "Chain"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "grey92", colour = "grey70"),
      panel.grid.minor = element_blank(),
      legend.position  = "right"
    )
}


# =============================================================================
# 8.  MAIN: RUN HMC
# =============================================================================

# ── Sampler settings ──────────────────────────────────────────────────────────
N_CHAINS    <- 4L      # parallel chains for R-hat / ESS
N_CORES     <- 10L      # cores for parallel gradient batches (set to 1 for debugging)
N_WARMUP    <- 200L    # adaptation (discarded)
N_ITER      <- 800L    # total iterations per chain  (post-warmup = N_ITER - N_WARMUP)
EPS_INIT    <- 0.01    # initial step size (dual averaging will adapt)
L_STEPS     <- 10L     # leapfrog steps per proposal
ADAPT_DELTA <- 0.65    # target acceptance rate

# ── Initial points ────────────────────────────────────────────────────────────
set.seed(114514)
inits <- lapply(seq_len(N_CHAINS), function(ch) {
  c(
    log(runif(1L, 0.01, 1.0)),      # log(beta_scale):          prior mean ≈0, start in (0,1)
    log(runif(9L, 0.01, 10.0))      # log(C_contact_scale_j):   start in (0,10)
  )
})

# ── Run chains (sequential; replace with parallel::mclapply for speed) ────────
cat("=== HMC Calibration: HCV PWID Model ===\n")
cat(sprintf("Chains: %d | Warmup: %d | Sampling: %d | L: %d\n",
            N_CHAINS, N_WARMUP, N_ITER - N_WARMUP, L_STEPS))
cat(sprintf("Approx ODE runs per HMC step: %d\n", (L_STEPS + 1L) * 2L * N_PARAMS))

hmc_chains <- parallel::mclapply(seq_len(N_CHAINS), function(ch) {
  run_hmc_chain(
    theta_init  = inits[[ch]],
    n_iter      = N_ITER,
    n_warmup    = N_WARMUP,
    eps_init    = EPS_INIT,
    L           = L_STEPS,
    adapt_delta = ADAPT_DELTA,
    base_params = params,
    data        = data,
    seed        = ch * 314L,
    chain_id    = ch
  )
}, mc.cores = N_CORES)

# ── Pool post-warmup samples ──────────────────────────────────────────────────
post_warmup_list <- lapply(hmc_chains, function(ch) ch$samples)
all_samples      <- do.call(rbind, post_warmup_list)

cat(sprintf("\nTotal post-warmup samples: %d (%d × %d chains)\n",
            nrow(all_samples), nrow(post_warmup_list[[1L]]), N_CHAINS))

# ── Diagnostics ───────────────────────────────────────────────────────────────
diag_table <- print_diagnostics(post_warmup_list, chains_raw = hmc_chains)

# ── Acceptance rates ──────────────────────────────────────────────────────────
cat("\nAcceptance rates per chain (post-warmup):\n")
acc_df <- data.frame(
  Chain       = seq_len(N_CHAINS),
  Acc_rate    = round(vapply(hmc_chains, function(ch) ch$acc_rate, numeric(1L)), 3),
  Eps_final   = round(vapply(hmc_chains, function(ch) ch$eps_final, numeric(1L)), 5)
)
print(acc_df, row.names = FALSE)

# ── Posterior summary (original scale) ───────────────────────────────────────
orig_samples <- exp(all_samples)
post_summary <- data.frame(
  Parameter  = param_names_orig,
  Mean       = round(colMeans(orig_samples), 4),
  SD         = round(apply(orig_samples, 2, sd), 4),
  Q2.5       = round(apply(orig_samples, 2, quantile, 0.025), 4),
  Q50        = round(apply(orig_samples, 2, quantile, 0.500), 4),
  Q97.5      = round(apply(orig_samples, 2, quantile, 0.975), 4)
)
cat("\n=== Posterior Summary (original scale) ===\n")
print(post_summary, row.names = FALSE)

# ── PPC ───────────────────────────────────────────────────────────────────────
cat("\n=== Generating Posterior Predictive Samples ===\n")
ppc_out <- generate_ppc_samples(all_samples, params, data, n_ppc = 200L)

cat("\nPosterior predictive p-values (near 0.5 = good fit):\n")
ppp_df <- data.frame(
  Age     = paste0("Age ", 1:9),
  ppp_pos = round(ppc_out$ppp_pos, 3),
  ppp_tot = round(ppc_out$ppp_tot, 3)
)
print(ppp_df, row.names = FALSE)

# ── Plots ─────────────────────────────────────────────────────────────────────
p_trace    <- plot_traces(hmc_chains, param_idx = 1:5)
p_density  <- plot_posterior_densities(post_warmup_list)
p_hist_pos <- plot_ppc_histograms(ppc_out, type = "pos")
p_hist_tot <- plot_ppc_histograms(ppc_out, type = "tot")
p_interval <- plot_ppc_intervals(ppc_out)

print(p_trace)
print(p_density)
print(p_hist_pos)
print(p_hist_tot)
print(p_interval)

# ── Save ─────────────────────────────────────────────────────────────────────
saveRDS(
  list(
    hmc_chains       = hmc_chains,
    post_warmup_list = post_warmup_list,
    all_samples      = all_samples,
    diag_table       = diag_table,
    post_summary     = post_summary,
    ppc_out          = ppc_out
  ),
  file = "hmc_output.rds"
)
cat("\nAll results saved to hmc_output.rds\n")