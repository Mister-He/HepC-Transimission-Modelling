# =============================================================================
# DATA / SIMULATION SETTINGS
# =============================================================================
# Observation-model dispersion settings.
# These can be tuned or moved into `data` if you want to estimate them.
sigma_N <- 1.0
phi_overdisp <- 50.0

# =============================================================================
# MAIN: RUN HMC
# =============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)

# Observations
obs_pos <- c(11, 51, 99, 141, 209, 339, 437, 367, 351)
obs_tot <- c(223, 572, 790, 765, 747, 810, 770, 658, 803)
N_total_obs <- sum(obs_tot)

param_names_log <- c(
  "log_beta_scale", "mu_hier", "log_sigma_hier",
  paste0("eta_", 1:8)
)
param_names_orig <- c(
  "beta_scale", "mu_hier", "sigma_hier",
  paste0("C_contact_scale_", 1:8)
)

source("setup.R")  
source("HMC_core.R") 
source("HMC_conv.R")  

# ── Sampler settings ──────────────────────────────────────────────────────────
N_CHAINS    <- 4L      # parallel chains for R-hat / ESS
N_CORES     <- 4L      # cores for parallel gradient batches (set to 1 for debugging)
N_WARMUP    <- 200L    # adaptation (discarded)
N_ITER      <- 1000L   # total iterations per chain  (post-warmup = N_ITER - N_WARMUP)
N_PARAMS    <- 11L     # number of parameters being sampled (log scale)
EPS_INIT    <- 0.01    # initial step size (dual averaging will adapt)
L_STEPS     <- 10L     # leapfrog steps per proposal
ADAPT_DELTA <- 0.65    # target acceptance rate

# ── Initial points ────────────────────────────────────────────────────────────
set.seed(42)
inits <- lapply(seq_len(N_CHAINS), function(ch) {
  c(
    log(runif(1L, 0.01, 1.0)),    # log(beta_scale):  start in (0,1)
    rnorm(1L, 0.0, 0.5),          # mu_hier:          log-scale mean near 0
    rnorm(1L, log(0.5), 0.3),     # log(sigma_hier):  start near prior center
    rnorm(8L, 0.0, 0.5)           # eta[1:8]:         standardised deviates near 0
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
}, mc.cores = N_CHAINS)

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
orig_samples <- theta_to_orig(all_samples)
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
set.seed(114514)
cat("\n=== Generating Posterior Predictive Samples ===\n")
ppc_out <- generate_ppc_samples(all_samples, params, data, n_ppc = 600L)

cat("\nPosterior predictive p-values (near 0.5 = good fit):\n")
ppp_df <- data.frame(
  Age     = paste0("Age ", 1:9),
  ppp_pos = round(ppc_out$ppp_pos, 3),
  ppp_tot = round(ppc_out$ppp_tot, 3)
)
print(ppp_df, row.names = FALSE)

# ── Plots ─────────────────────────────────────────────────────────────────────
p_trace    <- plot_traces(hmc_chains, param_idx = 1:N_PARAMS)
p_density  <- plot_posterior_densities(post_warmup_list)
p_hist_pos <- plot_ppc_histograms(ppc_out, type = "pos")
p_hist_tot <- plot_ppc_histograms(ppc_out, type = "tot")
p_interval <- plot_ppc_intervals(ppc_out)
p_prev_int <- plot_ppc_prevalence_intervals(ppc_out)

print(p_trace)
print(p_density)
print(p_hist_pos)
print(p_hist_tot)
print(p_interval)
print(p_prev_int)


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
  file = "hmc_output_full.rds"
)
cat("\nAll results saved to hmc_output_full.rds\n")