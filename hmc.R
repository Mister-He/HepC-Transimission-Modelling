# =============================================================================
# MAIN: RUN HMC
# =============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)

# Observations
obs_prev <- c(0.1118421, 0.1470588, 0.1933842, 0.2507599, 0.2899083, 0.3596059, 0.5025295, 0.5061728, 0.4534314, 0.3544304)
obs_tot <- c(99, 552, 692, 763, 704, 847, 994, 847, 781, 409)
N_AGE <- length(obs_tot)
N_CONTACT <- N_AGE

param_names_log <- c(
  "log_beta_scale",
  paste0("log_C_contact_scale_", seq_len(N_CONTACT))
)
param_names_orig <- c(
  "beta_scale",
  paste0("C_contact_scale_", seq_len(N_CONTACT))
)

source("setup.R")  
source("HMC_core.r")
source("HMC_conv.R")  
data$count_log_sd <- 0.35
data$prev_logit_sd <- 0.25

# ── Sampler settings ──────────────────────────────────────────────────────────
N_CHAINS    <- 4L      # parallel chains for R-hat / ESS
N_CORES     <- 4L      # cores for parallel gradient batches (set to 1 for debugging)
N_WARMUP    <- 200L    # adaptation (discarded)
N_ITER      <- 1000L   # total iterations per chain  (post-warmup = N_ITER - N_WARMUP)
N_PARAMS    <- length(param_names_log) # number of parameters being sampled (log scale)
EPS_INIT    <- 0.01    # initial step size (dual averaging will adapt)
L_STEPS     <- 10L     # leapfrog steps per proposal
ADAPT_DELTA <- 0.65    # target acceptance rate

# ── Initial points ────────────────────────────────────────────────────────────
set.seed(114514)
inits <- lapply(seq_len(N_CHAINS), function(ch) {
  c(
    rnorm(1L, 0.0, 0.25),       # log_beta_scale:        shared inflow multiplier
    rnorm(N_CONTACT, 0.0, 0.25) # log_C_contact_scale_i: row-specific contact scales
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
  Age     = paste0("Age ", seq_len(N_AGE)),
  ppp_prev = round(ppc_out$ppp_prev, 3),
  ppp_tot = round(ppc_out$ppp_tot, 3)
)
print(ppp_df, row.names = FALSE)

# ── Plots ─────────────────────────────────────────────────────────────────────
p_trace    <- plot_traces(hmc_chains, param_idx = 1:N_PARAMS)
p_density  <- plot_posterior_densities(post_warmup_list)
p_hist_prev <- plot_ppc_histograms(ppc_out, type = "prev")
p_hist_tot <- plot_ppc_histograms(ppc_out, type = "tot")
p_interval <- plot_ppc_intervals(ppc_out)
p_prev_int <- plot_ppc_prevalence_intervals(ppc_out)

print(p_trace)
print(p_density)
print(p_hist_prev)
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
result = readRDS("hmc_output_full.rds")
hmc_chains = result$hmc_chains
post_warmup_list = result$post_warmup_list
all_samples = result$all_samples
diag_table = result$diag_table
post_summary = result$post_summary
ppc_out = result$ppc_out
