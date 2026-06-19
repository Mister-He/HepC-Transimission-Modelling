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
#   Contact scales (manual): c(0.2, 0.06, 0.04, 0.04, 0.07, 0.2, 0.6, 0.1)
#                             [rows 1-8; row 9 is the reference = 1]
#
#   Beta scale factors (manual): c(1, 1, 0.8, 0.1, 0.1, 0.1, 1, 2, 40) * 0.7
# =============================================================================

# Contact starting values — non-centred parameterisation
contact_scales_manual <- c(0.2, 0.06, 0.04, 0.04, 0.07, 0.2, 0.6, 0.1)
log_s   <- log(contact_scales_manual)
mu_0    <- mean(log_s)
sigma_0 <- max(sd(log_s), 0.1)
eta_0   <- (log_s - mu_0) / sigma_0

THETA_CONTACT_INIT <- c(mu_0, log(sigma_0), eta_0)
names(THETA_CONTACT_INIT) <- CONTACT_PARAM_NAMES

# Beta starting values — log deviations from baseline
beta_scales_manual  <- c(1, 1, 0.8, 0.1, 0.1, 0.1, 1, 2, 40) * 0.7
THETA_BETA_INIT     <- log(beta_scales_manual)
names(THETA_BETA_INIT) <- BETA_PARAM_NAMES

# =============================================================================
# Sanity checks on initial log-posteriors
# =============================================================================
cat("=== Nested Calibration: HCV PWID Model ===\n")
cat(sprintf("Stage 1: %d contact-rate params | Stage 2: %d beta-deviation params\n",
            N_CONTACT, N_BETA))

lp0_s1 <- log_posterior_nested(THETA_CONTACT_INIT, params, data)
lp0_s2 <- log_posterior_stage2(
  THETA_BETA_INIT,
  build_contact_params(THETA_CONTACT_INIT, params),
  params, data
)
cat(sprintf("Initial log-posterior  Stage 1: %.3f\n", lp0_s1))
cat(sprintf("Initial log-posterior  Stage 2: %.3f\n", lp0_s2))

if (!is.finite(lp0_s1)) {
  stop("Stage 1 initial log-posterior is non-finite. Adjust THETA_CONTACT_INIT.")
}
if (!is.finite(lp0_s2)) {
  warning("Stage 2 initial log-posterior is non-finite — Stage 2 will fall back to zero-centred beta init.")
  THETA_BETA_INIT <- rep(0, N_BETA)
  names(THETA_BETA_INIT) <- BETA_PARAM_NAMES
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
# Derive Stage 2 starting values from Stage 1 posterior mean
# Run the nested back-calculation once at the Stage 1 posterior mean to obtain
# an adjusted beta; use its log-deviation from baseline as theta_b initialisation.
# =============================================================================
pm_s1_mean <- build_contact_params(theta_contact_mean, params)
adj_s1     <- nested_beta_adjust(pm_s1_mean, data, obs_prop)

if (!is.null(adj_s1)) {
  theta_beta_start <- log(adj_s1$pm_adj$beta / params$beta)
  cat("\nBack-calculated log_beta_delta at Stage 1 posterior mean:\n")
  print(round(theta_beta_start, 3))
  cat("Implied beta scale factors:\n")
  print(round(exp(theta_beta_start), 3))
} else {
  warning("Nested beta adjustment failed at Stage 1 mean; using manual beta init for Stage 2.")
  theta_beta_start <- THETA_BETA_INIT
}

# Verify Stage 2 start is finite before launching
pm_contact_fixed <- build_contact_params(theta_contact_mean, params)
lp_s2_start <- log_posterior_stage2(theta_beta_start, pm_contact_fixed, params, data)
if (!is.finite(lp_s2_start)) {
  warning("Stage 2 derived start is non-finite; falling back to THETA_BETA_INIT.")
  theta_beta_start <- THETA_BETA_INIT
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
    log_post_fn   = function(th) log_posterior_stage2(th, pm_contact_fixed, params, data),
    n_iter        = N_ITER_S2,
    n_warmup      = N_WARMUP_S2,
    param_names   = BETA_PARAM_NAMES,
    chain_id      = k,
    verbose_every = 1000L
  )
}

samples_s2       <- do.call(rbind, lapply(chains_s2, `[[`, "samples"))
theta_beta_mean  <- colMeans(samples_s2)
theta_beta_sd    <- apply(samples_s2, 2, sd)
beta_final       <- params$beta * exp(theta_beta_mean)

cat(sprintf("\nStage 2: %d post-warmup samples | acceptance rates: %s\n",
            nrow(samples_s2),
            paste(sprintf("%.3f", sapply(chains_s2, `[[`, "acc_rate")), collapse = ", ")))
cat("Posterior means (log_beta_delta):\n"); print(round(theta_beta_mean, 3))
cat("Posterior SDs  (log_beta_delta):\n");  print(round(theta_beta_sd,   3))
cat("Implied beta (exp scale × baseline):\n"); print(round(beta_final, 2))

saveRDS(
  list(chains = chains_s2, samples = samples_s2,
       theta_post_mean = theta_beta_mean,
       theta_post_sd   = theta_beta_sd,
       beta_final      = beta_final),
  file = file.path(OUT_DIR, "nc_stage2_results.rds")
)

# =============================================================================
# COMBINED POSTERIOR PREDICTIVE CHECK
# =============================================================================
cat("\n--- Posterior predictive check at combined posterior means ---\n")

pm_final       <- build_contact_params(theta_contact_mean, params)
pm_final$beta  <- beta_final
out_final      <- run_sim(pm_final, data)
q_final        <- compute_age_quantities(as.numeric(out_final[nrow(out_final), -1L]))

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
