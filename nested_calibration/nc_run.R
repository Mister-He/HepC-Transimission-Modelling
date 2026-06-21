# =============================================================================
# nc_run.R  —  Two-stage nested calibration: main entry point
#
# Run from the REPO ROOT:
#   source("nested_calibration/nc_run.R")
# or open this file in RStudio and set working dir to repo root first.
#
# What this does:
#   Stage 1 — samples 10 contact-rate parameters via AMH.  Inside every
#              likelihood call, beta is back-calculated proportionally so the
#              model age distribution matches the observed survey proportions.
#              Only the BetaBinomial prevalence likelihood enters the posterior;
#              the age-proportion target is enforced by construction.
#
#   Stage 2 — fixes contact rates at the Stage 1 posterior mean and samples
#              9 log(beta) deviation parameters.  The full Multinomial +
#              BetaBinomial likelihood is used so that residual age-proportion
#              variation is accounted for.
#
# Outputs written to nested_calibration/output/:
#   nc_stage1_results.rds    — Stage 1 chain objects + pooled samples
#   nc_stage2_results.rds    — Stage 2 chain objects + pooled samples
#   nc_combined_results.rds  — combined summary + PPC at posterior means
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# ── Dependencies ──────────────────────────────────────────────────────────────
# All three must be sourced from the repo root.
source("setup.R")                       # params, data, idx(), run_sim()
source("HMC_core.r")                    # compute_age_quantities(), log_beta_binomial()
source("nested_calibration/nc_core.R")  # all NC functions + constants

OUT_DIR <- "nested_calibration/output"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# MCMC settings
# Increase N_ITER_* for a production run; keep small for a smoke-test.
# =============================================================================
N_CHAINS      <- 2L
N_ITER_S1     <- 15000L    # Stage 1 total iterations per chain
N_WARMUP_S1   <-  5000L    # Stage 1 warmup (discarded)
N_ITER_S2     <- 10000L    # Stage 2 total iterations per chain
N_WARMUP_S2   <-  3000L    # Stage 2 warmup (discarded)

# =============================================================================
# Starting values: derived from manual fitting described in solutions.md
#
#   Contact scales (manual): c(0.06, 0.04, 0.04, 0.07, 0.2, 0.6, 0.1, 0.5)
#                             [rows 2-9; row 1 is now the reference = 1]
#     These are approximate starting guesses relative to row 1.
#     Rows 7-8 are initialised with moderate-to-high scales (0.6, 0.1→0.5) to
#     allow Stage 1 to explore the high-prevalence region for those age groups.
#
#   Beta scale factors (manual): c(1, 1, 0.8, 0.1, 0.1, 0.1, 1, 2, 40) * 0.7
# =============================================================================

# Contact starting values for rows 2-9 (row 1 is the reference, fixed at 1.0)
contact_scales_manual <- c(0.06, 0.04, 0.04, 0.07, 0.2, 0.6, 0.1, 0.5)
log_s   <- log(contact_scales_manual)
mu_0    <- mean(log_s)
sigma_0 <- max(sd(log_s), 0.1)
eta_0   <- (log_s - mu_0) / sigma_0

THETA_CONTACT_INIT <- c(mu_0, log(sigma_0), eta_0)
names(THETA_CONTACT_INIT) <- CONTACT_PARAM_NAMES

# =============================================================================
# Sanity checks on initial log-posteriors
# =============================================================================
cat("=== Nested Calibration: HCV PWID Model ===\n")
cat(sprintf("Stage 1: %d contact-rate params | Stage 2: %d beta-deviation params\n",
            N_CONTACT, N_BETA))

lp0_s1 <- log_posterior_nested(THETA_CONTACT_INIT, params, data)
cat(sprintf("Initial log-posterior  Stage 1: %.3f\n", lp0_s1))

if (!is.finite(lp0_s1)) {
  stop("Stage 1 initial log-posterior is non-finite. Adjust THETA_CONTACT_INIT.")
}

# =============================================================================
# STAGE 1: Contact-rate calibration (nested beta back-calculation)
# =============================================================================
cat(sprintf("\n--- Stage 1: %d chains × %d iter (%d warmup) ---\n",
            N_CHAINS, N_ITER_S1, N_WARMUP_S1))

set.seed(114514L)
chains_s1 <- vector("list", N_CHAINS)
for (k in seq_len(N_CHAINS)) {
  theta_k <- THETA_CONTACT_INIT + rnorm(N_CONTACT, sd = 0.05)  # per-chain jitter
  chains_s1[[k]] <- run_adaptive_mh(
    theta_init    = theta_k,
    log_post_fn   = function(th) log_posterior_nested(th, params, data),
    n_iter        = N_ITER_S1,
    n_warmup      = N_WARMUP_S1,
    param_names   = CONTACT_PARAM_NAMES,
    chain_id      = k,
    verbose_every = 2000L
  )
}

samples_s1          <- do.call(rbind, lapply(chains_s1, `[[`, "samples"))
theta_contact_mean  <- colMeans(samples_s1)
theta_contact_sd    <- apply(samples_s1, 2, sd)

cat(sprintf("\nStage 1: %d post-warmup samples | acceptance rates: %s\n",
            nrow(samples_s1),
            paste(sprintf("%.3f", sapply(chains_s1, `[[`, "acc_rate")), collapse = ", ")))
cat("Posterior means (contact):\n"); print(round(theta_contact_mean, 4))
cat("Posterior SDs  (contact):\n");  print(round(theta_contact_sd,   4))

# Implied contact-rate scales at posterior mean
scales_mean <- constrain_contact(theta_contact_mean)$C_contact_scale
cat("Implied C_contact_scale at posterior mean:\n")
print(round(scales_mean, 4))

saveRDS(
  list(chains = chains_s1, samples = samples_s1,
       theta_post_mean = theta_contact_mean,
       theta_post_sd   = theta_contact_sd,
       scales_mean     = scales_mean),
  file = file.path(OUT_DIR, "nc_stage1_results.rds")
)

# =============================================================================
# Derive Stage 2 baseline and starting values from Stage 1 posterior mean
#
# beta_nested: the unique beta that makes model age proportions = obs_prop
#              given the Stage 1 posterior-mean contact matrix.
# theta_b = 0 → keep this beta exactly; Stage 2 samples small residual
#              deviations on the log scale.
# =============================================================================
pm_s1_mean <- build_contact_params(theta_contact_mean, params)
adj_s1     <- nested_beta_adjust(pm_s1_mean, data, obs_prop)

if (!is.null(adj_s1)) {
  beta_nested      <- adj_s1$pm_adj$beta
  theta_beta_start <- rep(0.0, N_BETA)           # zero deviations from Stage-1 result
  names(theta_beta_start) <- BETA_PARAM_NAMES
  cat("\nStage-1 nested-adjusted beta (used as Stage-2 baseline):\n")
  print(round(beta_nested, 3))
  cat("Implied scale vs params$beta:\n")
  print(round(beta_nested / params$beta, 3))
} else {
  stop("Nested beta adjustment failed at Stage 1 posterior mean. Cannot start Stage 2.")
}

# Verify Stage 2 start is finite before launching
pm_contact_fixed <- build_contact_params(theta_contact_mean, params)
lp_s2_start <- log_posterior_stage2(theta_beta_start, pm_contact_fixed, beta_nested, data)
cat(sprintf("Stage 2 initial log-posterior: %.3f\n", lp_s2_start))
if (!is.finite(lp_s2_start)) {
  stop("Stage 2 initial log-posterior is non-finite. Check beta_nested and contact params.")
}

# =============================================================================
# STAGE 2: Beta refinement (contact rates fixed at Stage 1 posterior mean)
# =============================================================================
cat(sprintf("\n--- Stage 2: %d chains × %d iter (%d warmup) ---\n",
            N_CHAINS, N_ITER_S2, N_WARMUP_S2))

chains_s2 <- vector("list", N_CHAINS)
for (k in seq_len(N_CHAINS)) {
  theta_k <- theta_beta_start + rnorm(N_BETA, sd = 0.05)
  chains_s2[[k]] <- run_adaptive_mh(
    theta_init    = theta_k,
    log_post_fn   = function(th) log_posterior_stage2(th, pm_contact_fixed, beta_nested, data),
    n_iter        = N_ITER_S2,
    n_warmup      = N_WARMUP_S2,
    param_names   = BETA_PARAM_NAMES,
    chain_id      = k,
    verbose_every = 1000L
  )
}

samples_s2      <- do.call(rbind, lapply(chains_s2, `[[`, "samples"))
theta_beta_mean <- colMeans(samples_s2)
theta_beta_sd   <- apply(samples_s2, 2, sd)

# Stage 2 now uses nested adjustment internally, so beta_final is NOT simply
# beta_nested * exp(theta_beta_mean).  Run the nested adjustment at the
# posterior mean theta_b to recover the actual final beta.
pm_for_final       <- pm_contact_fixed
pm_for_final$beta  <- beta_nested * exp(theta_beta_mean)
adj_final          <- nested_beta_adjust(pm_for_final, data, obs_prop)
if (!is.null(adj_final)) {
  beta_final <- adj_final$pm_adj$beta
} else {
  warning("Final nested adjustment failed; using beta_nested * exp(theta_beta_mean).")
  beta_final <- beta_nested * exp(theta_beta_mean)
}

cat(sprintf("\nStage 2: %d post-warmup samples | acceptance rates: %s\n",
            nrow(samples_s2),
            paste(sprintf("%.3f", sapply(chains_s2, `[[`, "acc_rate")), collapse = ", ")))
cat("Posterior means (log_beta_delta):\n"); print(round(theta_beta_mean, 3))
cat("Posterior SDs  (log_beta_delta):\n");  print(round(theta_beta_sd,   3))
cat("Final beta (after nested adjustment at posterior mean):\n"); print(round(beta_final, 2))

saveRDS(
  list(chains = chains_s2, samples = samples_s2,
       theta_post_mean = theta_beta_mean,
       theta_post_sd   = theta_beta_sd,
       beta_nested     = beta_nested,
       beta_final      = beta_final),
  file = file.path(OUT_DIR, "nc_stage2_results.rds")
)

# =============================================================================
# COMBINED POSTERIOR PREDICTIVE CHECK
# =============================================================================
cat("\n--- Posterior predictive check at combined posterior means ---\n")

# pm_for_final already has the correct contact params; adj_final has the ODE result.
pm_final <- adj_final$pm_adj
q_final  <- adj_final$model_q
if (is.null(q_final)) {
  pm_final      <- build_contact_params(theta_contact_mean, params)
  pm_final$beta <- beta_final
  out_final     <- run_sim(pm_final, data)
  q_final       <- compute_age_quantities(as.numeric(out_final[nrow(out_final), -1L]))
}

ppc_df <- data.frame(
  age_group  = paste0("Age ", 1:9),
  obs_prop   = round(obs_prop,  4),
  model_prop = round(q_final$p_age, 4),
  prop_err   = round(abs(q_final$p_age - obs_prop), 4),
  obs_prev   = round(obs_prev,  4),
  model_prev = round(q_final$q_age, 4),
  prev_err   = round(abs(q_final$q_age - obs_prev), 4)
)
cat("Age proportions and HCV prevalences:\n")
print(ppc_df)

cat(sprintf(
  "\nMean absolute error — proportion: %.4f | prevalence: %.4f\n",
  mean(ppc_df$prop_err), mean(ppc_df$prev_err)
))

saveRDS(
  list(
    stage1_samples    = samples_s1,
    stage2_samples    = samples_s2,
    theta_contact     = theta_contact_mean,
    theta_beta        = theta_beta_mean,
    beta_nested       = beta_nested,
    beta_final        = beta_final,
    pm_final          = pm_final,
    q_final           = q_final,
    ppc_df            = ppc_df,
    scales_mean       = scales_mean
  ),
  file = file.path(OUT_DIR, "nc_combined_results.rds")
)

cat(sprintf("\nAll outputs saved to %s/\n", OUT_DIR))
cat("Run nested_calibration/nc_diagnostics.R for plots and convergence diagnostics.\n")
