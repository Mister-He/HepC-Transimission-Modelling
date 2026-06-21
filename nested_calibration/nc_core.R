# =============================================================================
# nc_core.R  —  Core functions for nested two-stage calibration
#
# Parameter space (19 total):
#   Stage 1 — contact rates (10):
#     theta[1]     = mu_hier          log-mean of contact row scalings
#     theta[2]     = log_sigma_hier   log-SD of contact row scalings
#     theta[3:10]  = eta[1:8]         non-centred deviates; row 9 fixed = 1
#
#   Stage 2 — inflow deviations (9):
#     theta[1:9]   = log_beta_delta[1:9]
#                    log(beta[i] / baseline_beta[i]) for each age group
#
# Calibration targets:
#   obs_prop[i] = obs_tot[i] / sum(obs_tot)  — age proportion of PWID survey
#   obs_prev[i] = obs_pos[i] / obs_tot[i]    — HCV prevalence by age group
#
# Observation model:
#   Stage 1 (contact, nested): BetaBinomial(obs_pos[i] | obs_tot[i], q_age[i], phi)
#                              Age proportions enforced by construction (beta rescale)
#   Stage 2 (beta):            Multinomial(obs_tot | p_age) + BetaBinomial(prev)
#
# Must source setup.R and HMC_core.r before this file (for run_sim, compute_age_quantities).
# =============================================================================

# ── Observed targets ──────────────────────────────────────────────────────────
obs_pos  <- c(11,  51,  99, 141, 209, 339, 437, 367, 351)
obs_tot  <- c(223, 572, 790, 765, 747, 810, 770, 658, 803)
obs_prev <- obs_pos / obs_tot
obs_prop <- obs_tot / sum(obs_tot)        # target age distribution

# ── Constants ─────────────────────────────────────────────────────────────────
PHI_OVERDISP     <- 50.0   # BetaBinomial overdispersion (fixed; tunable)
SIGMA_BETA_PRIOR <- 0.75   # prior SD on log_beta_delta (covers ±2.25 log-units around Stage-1 beta)

N_CONTACT <- 10L           # contact-rate parameters (Stage 1)
N_BETA    <-  9L           # beta log-deviation parameters (Stage 2)
N_JOINT   <- 19L           # total across both stages

CONTACT_PARAM_NAMES <- c("mu_hier", "log_sigma_hier", paste0("eta_", 1:8))
BETA_PARAM_NAMES    <- paste0("log_beta_delta_", 1:9)
JOINT_PARAM_NAMES   <- c(CONTACT_PARAM_NAMES, BETA_PARAM_NAMES)

# =============================================================================
# 1.  PARAMETER TRANSFORMS
# =============================================================================

# Non-centred log-normal hierarchy for contact row scalings.
# Reference: row 1 fixed at scale 1.0 (base row-sum ≈ 60 contacts — smallest row).
# Fixing the row with the smallest base allows rows 7-9 to achieve scales > 1,
# which is needed to reproduce the high observed HCV prevalences for those groups.
# Previously row 9 was the reference, but rows 8 and 9 have identical base matrices,
# and fixing row 9 at 1.0 forced rows 7-8 to very low scales (0.17-0.23), preventing
# the model from reaching the 57% prevalence observed in age groups 7-8.
# C_contact_scale[1] = 1  (reference — identifiability anchor)
# C_contact_scale[j] = exp(mu_hier + sigma_hier * eta[j-1])  for j = 2..9
constrain_contact <- function(theta_c) {
  mu    <- theta_c[1L]
  sigma <- exp(theta_c[2L])
  eta   <- theta_c[3:10]
  list(
    C_contact_scale = c(1.0, exp(mu + sigma * eta)),
    mu_hier         = mu,
    sigma_hier      = sigma,
    eta             = eta
  )
}

# Build a full params list with scaled contact matrix rows 2:9 (row 1 is the reference).
build_contact_params <- function(theta_c, base_params) {
  p  <- constrain_contact(theta_c)
  pm <- base_params
  for (j in 2:9) {
    pm$C_contact[j, ] <- base_params$C_contact[j, ] * p$C_contact_scale[j]
  }
  pm
}

# Apply log-beta deviations: beta[i] = baseline * exp(log_delta[i])
apply_beta_delta <- function(theta_b, base_params) {
  base_params$beta * exp(theta_b)
}

# Back-transform Stage 1 posterior samples to interpretable contact scales.
contact_to_orig <- function(samps) {
  mu    <- samps[, 1L]
  sigma <- exp(samps[, 2L])
  c_cols <- sapply(1:8, function(j) exp(mu + sigma * samps[, 2L + j]))
  colnames(c_cols) <- paste0("C_contact_scale_", 2:9)
  cbind(mu_hier = mu, sigma_hier = sigma, C_contact_scale_1 = 1.0, c_cols)
}

# =============================================================================
# 2.  NESTED BETA BACK-CALCULATION
# =============================================================================
# Given pm with contact rates already set, run ODE, compute model age
# proportions, proportionally rescale beta so model proportions match
# obs_prop, then run ODE a second time.  Returns list(pm_adj, model_q)
# from the second run, or NULL on any failure.
#
# Why one rescale + one re-run suffices:
#   At equilibrium, N_i ∝ beta[i] / effective_outflow[i].  The proportional
#   rescale beta_adj[i] = beta[i] * (obs_prop[i] / model_prop[i]) directly
#   targets the mismatch.  A single pass is accurate when deviations are
#   moderate; add n_iter > 1 for tighter enforcement if needed.

nested_beta_adjust <- function(pm, data, target_prop, n_iter = 1L) {
  for (it in seq_len(n_iter)) {
    out <- tryCatch(run_sim(pm, data), error = function(e) NULL)
    if (is.null(out) || !is.matrix(out) || nrow(out) == 0L || any(!is.finite(out))) return(NULL)

    q <- compute_age_quantities(as.numeric(out[nrow(out), -1L]))
    if (is.null(q)) return(NULL)

    ratio <- target_prop / q$p_age
    if (any(!is.finite(ratio)) || any(ratio <= 0)) return(NULL)
    pm$beta <- pm$beta * ratio
  }

  # Final run with fully adjusted beta
  out2 <- tryCatch(run_sim(pm, data), error = function(e) NULL)
  if (is.null(out2) || !is.matrix(out2) || nrow(out2) == 0L || any(!is.finite(out2))) return(NULL)
  q2 <- compute_age_quantities(as.numeric(out2[nrow(out2), -1L]))
  if (is.null(q2)) return(NULL)

  list(pm_adj = pm, model_q = q2)
}

# =============================================================================
# 3.  LOG-PRIORS
# =============================================================================

# Contact rates: non-centred log-normal hierarchy (same as HMC_core.r)
log_prior_contact <- function(theta_c) {
  dnorm(theta_c[1L],    mean = 0.0,      sd = 1.0, log = TRUE) +   # mu_hier
  dnorm(theta_c[2L],    mean = log(0.5), sd = 0.5, log = TRUE) +   # log_sigma_hier
  sum(dnorm(theta_c[3:10], mean = 0.0,  sd = 1.0, log = TRUE))     # eta[1:8]
}

# Beta deviations: independent Normal centred at 0 (= baseline beta).
# SIGMA_BETA_PRIOR = 1.5 covers the manual-fit range:
#   log(c(1,1,0.8,0.1,0.1,0.1,1,2,40)*0.7) ≈ (-0.36, ..., 3.33)
log_prior_beta_delta <- function(theta_b) {
  sum(dnorm(theta_b, mean = 0.0, sd = SIGMA_BETA_PRIOR, log = TRUE))
}

# =============================================================================
# 4.  LOG-LIKELIHOODS
# =============================================================================

log_beta_binomial_nc <- function(x, size, prob, phi) {
  if (any(!is.finite(c(x, size, prob))) || !is.finite(phi) || phi <= 0) return(-Inf)
  if (any(x < 0) || any(x > size)) return(-Inf)
  prob <- pmin(pmax(prob, 1e-12), 1 - 1e-12)
  a <- prob * phi
  b <- (1 - prob) * phi
  sum(lchoose(size, x) + lbeta(x + a, size - x + b) - lbeta(a, b))
}

# =============================================================================
# 5.  LOG-POSTERIORS
# =============================================================================

# ── Stage 1: contact params only, beta derived via nested back-calculation ────
# Likelihood: BetaBinomial on prevalences.
# Age-proportion target is enforced structurally (not in the likelihood),
# so the sampler explores only the prevalence surface.
log_posterior_nested <- function(theta_c, base_params, data) {
  lp <- log_prior_contact(theta_c)
  if (!is.finite(lp)) return(-Inf)

  pm     <- build_contact_params(theta_c, base_params)
  result <- nested_beta_adjust(pm, data, obs_prop)
  if (is.null(result)) return(-Inf)

  ll <- log_beta_binomial_nc(obs_pos, obs_tot, result$model_q$q_age, PHI_OVERDISP)
  lp + ll
}

# ── Stage 2: beta params, contact rates fixed externally ─────────────────────
# Likelihood: BetaBinomial on prevalences ONLY.
# Age proportions are enforced structurally via nested_beta_adjust (same design
# as Stage 1), so no Multinomial term is needed.
#
# Why we dropped the Multinomial + free-ODE design:
#   At equilibrium, prevalence ≈ force-of-infection / outflow, which depends on
#   contact rates, NOT on beta.  So the BetaBinomial barely constrains beta, and
#   the Multinomial alone (with 9 free betas vs 8 independent proportion contrasts)
#   leaves a near-null space that produces catastrophic multimodality (R-hat > 10
#   for age groups 7-9, chains swapping modes).  Running the nested adjustment
#   inside Stage 2 breaks this degeneracy: each theta_b produces a UNIQUE final
#   beta (the one that exactly satisfies obs_prop), and the BetaBinomial then
#   identifies which direction in starting-beta space gives better prevalence fit.
#
# pm_contact: contact matrix fixed at Stage 1 posterior mean (C_contact scaled).
# beta_nested: the nested-adjusted beta from Stage 1 (NOT params$beta).
#   theta_b = 0 → no perturbation; nested adjustment gives back beta_nested.
log_posterior_stage2 <- function(theta_b, pm_contact, beta_nested, data) {
  lp <- log_prior_beta_delta(theta_b)
  if (!is.finite(lp)) return(-Inf)

  pm      <- pm_contact
  pm$beta <- beta_nested * exp(theta_b)   # perturbed starting beta

  # Enforce proportions exactly — same nested adjustment as Stage 1
  result <- nested_beta_adjust(pm, data, obs_prop)
  if (is.null(result)) return(-Inf)

  # Only BetaBinomial: proportions are satisfied structurally
  ll_prev <- log_beta_binomial_nc(obs_pos, obs_tot, result$model_q$q_age, PHI_OVERDISP)
  lp + ll_prev
}

# ── Joint reference: all 19 params simultaneously (for comparison / sanity) ───
# Likelihood: Multinomial + BetaBinomial.
# Useful for verifying the two-stage results agree with joint estimation.
log_posterior_joint <- function(theta, base_params, data) {
  lp <- log_prior_contact(theta[1:N_CONTACT]) +
        log_prior_beta_delta(theta[(N_CONTACT + 1L):N_JOINT])
  if (!is.finite(lp)) return(-Inf)

  pm      <- build_contact_params(theta[1:N_CONTACT], base_params)
  pm$beta <- apply_beta_delta(theta[(N_CONTACT + 1L):N_JOINT], base_params)

  out <- tryCatch(run_sim(pm, data), error = function(e) NULL)
  if (is.null(out) || !is.matrix(out) || nrow(out) == 0L || any(!is.finite(out))) return(-Inf)
  q <- compute_age_quantities(as.numeric(out[nrow(out), -1L]))
  if (is.null(q)) return(-Inf)

  ll_prop <- dmultinom(obs_tot, prob = q$p_age, log = TRUE)
  ll_prev <- log_beta_binomial_nc(obs_pos, obs_tot, q$q_age, PHI_OVERDISP)
  lp + ll_prop + ll_prev
}

# =============================================================================
# 6.  ADAPTIVE METROPOLIS-HASTINGS SAMPLER
# =============================================================================
# Uses an empirical covariance proposal (Haario et al. 2001) with a
# Robbins-Monro scale adjustment to target the desired acceptance rate.
# Works well for up to ~20 parameters; no gradient evaluations needed.
#
# Arguments:
#   theta_init     — starting point (named numeric vector)
#   log_post_fn    — function(theta) -> scalar log-posterior
#   n_iter         — total iterations (warmup + sampling)
#   n_warmup       — number of warmup / adaptation iterations (discarded)
#   target_accept  — target acceptance rate (0.234 is optimal for Gaussian)
#   adapt_every    — how often to update proposal covariance
#   init_scale     — initial isotropic proposal SD; defaults to 2.38/sqrt(d)
#   param_names    — column names for the sample matrix
#   chain_id       — integer label for progress messages
#   verbose_every  — print progress every this many iterations

run_adaptive_mh <- function(
  theta_init,
  log_post_fn,
  n_iter,
  n_warmup,
  target_accept = 0.234,
  adapt_every   = 50L,
  init_scale    = NULL,
  param_names   = NULL,
  chain_id      = 1L,
  verbose_every = 500L
) {
  n_p <- length(theta_init)
  if (is.null(init_scale)) init_scale <- 2.38 / sqrt(n_p)
  if (is.null(param_names)) param_names <- paste0("p_", seq_len(n_p))

  samples  <- matrix(NA_real_, n_iter, n_p, dimnames = list(NULL, param_names))
  accepted <- logical(n_iter)
  lp_trace <- numeric(n_iter)

  theta   <- theta_init
  lp_curr <- log_post_fn(theta)
  if (!is.finite(lp_curr)) {
    stop(sprintf("[Chain %d] Initial log-posterior is non-finite (%.3f). Fix starting values.",
                 chain_id, lp_curr))
  }

  # Proposal: start with scaled identity; switch to adapted covariance
  Sigma    <- (init_scale^2) * diag(n_p)
  L_chol   <- chol(Sigma)    # upper-triangular: t(L) %*% L = Sigma
  log_scale <- 0.0            # log of the global scale multiplier

  samp_buf <- matrix(theta_init, nrow = 1, ncol = n_p)   # accumulates draws for cov estimation

  cat(sprintf(
    "[Chain %d] AMH start | n_iter=%d n_warmup=%d n_params=%d init_lp=%.2f\n",
    chain_id, n_iter, n_warmup, n_p, lp_curr
  ))

  for (iter in seq_len(n_iter)) {

    # ── Propose ──────────────────────────────────────────────────────────────
    z          <- rnorm(n_p)
    theta_prop <- theta + as.numeric(crossprod(L_chol, z))   # t(L_chol) %*% z
    lp_prop    <- log_post_fn(theta_prop)

    # ── Accept / reject ───────────────────────────────────────────────────────
    log_alpha  <- lp_prop - lp_curr
    accept     <- is.finite(log_alpha) && (log_alpha >= 0 || log(runif(1L)) < log_alpha)
    if (accept) {
      theta   <- theta_prop
      lp_curr <- lp_prop
    }
    accepted[iter]  <- accept
    samples[iter, ] <- theta
    lp_trace[iter]  <- lp_curr

    # ── Adaptation (warmup only) ──────────────────────────────────────────────
    if (iter <= n_warmup) {
      samp_buf <- rbind(samp_buf, theta)

      if (iter %% adapt_every == 0L && nrow(samp_buf) > 2L * n_p) {
        # Robbins-Monro scale update toward target_accept
        acc_win    <- mean(accepted[max(1L, iter - adapt_every + 1L):iter])
        log_scale  <- log_scale + 0.5 * (acc_win - target_accept)

        S_emp    <- cov(samp_buf) + 1e-8 * diag(n_p)
        Sigma    <- exp(2 * log_scale) * (2.38^2 / n_p) * S_emp
        L_chol   <- tryCatch(chol(Sigma), error = function(e) {
          # Fallback: reset to scaled identity if covariance is degenerate
          chol(exp(2 * log_scale) * (init_scale^2) * diag(n_p))
        })
      }
    }

    # ── Progress ──────────────────────────────────────────────────────────────
    if (iter %% verbose_every == 0L) {
      phase      <- if (iter <= n_warmup) "warmup" else "sample"
      acc_recent <- mean(accepted[max(1L, iter - verbose_every + 1L):iter])
      cat(sprintf("[Chain %d] iter=%5d (%s) | lp=%8.2f | acc=%.3f\n",
                  chain_id, iter, phase, lp_curr, acc_recent))
    }
  }

  idx_post <- seq(n_warmup + 1L, n_iter)
  list(
    samples     = samples[idx_post, , drop = FALSE],
    samples_all = samples,
    lp_trace    = lp_trace,
    accepted    = accepted,
    acc_rate    = mean(accepted[idx_post]),
    Sigma_final = Sigma,
    n_warmup    = n_warmup,
    chain_id    = chain_id
  )
}
