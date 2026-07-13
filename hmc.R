# =============================================================================
# MAIN: RUN HMC
#
# Run from this directory with:
#   Rscript hmc.R [output_dir] [nelder_mead_fit_rds]
# =============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)

# Observations
obs_prev <- c(0.1118421, 0.1470588, 0.1933842, 0.2507599, 0.2899083, 0.3596059, 0.5025295, 0.5061728, 0.4534314, 0.3544304)
obs_prev_se <- c(0.09940967, 0.12867003, 0.16065544, 0.18991135, 0.21719265, 0.23775531, 0.24736428, 0.24839586, 0.24995891, 0.24000000)
obs_tot <- c(99, 552, 692, 763, 704, 847, 994, 847, 781, 409)
obs_pos <- round(obs_prev * obs_tot)  # HCV positives per age group (binomial numerator)
N_AGE <- length(obs_tot)
N_CONTACT <- N_AGE

param_names_log <- c(
  paste0("log_C_contact_scale_", seq_len(N_CONTACT)),
  paste0("log_tot_in_scaling_fct_", seq_len(N_CONTACT))
)
param_names_orig <- c(
  paste0("C_contact_scale_", seq_len(N_CONTACT)),
  paste0("tot_in_scaling_fct_", seq_len(N_CONTACT))
)

source("setup.R")  
source("HMC_core.r")
source("HMC_conv.R")  

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1L && nzchar(args[[1L]])) args[[1L]] else "."
nm_fit_file <- if (length(args) >= 2L && nzchar(args[[2L]])) {
  args[[2L]]
} else {
  file.path(output_dir, "nelder_mead_fit.rds")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

nm_fit <- NULL
if (file.exists(nm_fit_file)) {
  nm_fit <- readRDS(nm_fit_file)
  nm_center <- if (!is.null(nm_fit$theta_hat_spline)) nm_fit$theta_hat_spline else nm_fit$theta_hat
  if (length(nm_center) != 2L * N_AGE || any(!is.finite(nm_center))) {
    stop("Nelder-Mead fit does not contain a valid finite theta estimate")
  }
  cat(sprintf("Using Nelder-Mead estimates from %s to construct HMC priors\n", nm_fit_file))
} else {
  nm_center <- c(CONTACT_PRIOR_MEANS, TOT_IN_PRIOR_MEANS)
  cat(sprintf("Nelder-Mead fit not found at %s; using baseline spline prior center\n", nm_fit_file))
}

data$theta_prior <- make_spline_theta_prior(
  center = nm_center,
  n_knots = 5L,
  coef_sd = 0.35,
  residual_sd = 0.08
)
data$prev_logit_sd <- 0.1
data$sigma_pop     <- rep(0.10, N_AGE)
data$sigma_shape   <- 0.18
data$sigma_prev_shape <- 0.14
data$nu_shape      <- 6L

# ── Sampler settings ──────────────────────────────────────────────────────────
N_CHAINS    <- 4L      # parallel chains for R-hat / ESS
N_CORES     <- 4L      # cores for parallel gradient batches (set to 1 for debugging)
N_WARMUP    <- 500L    # adaptation (discarded)
N_ITER      <- 600L   # total iterations per chain  (post-warmup = N_ITER - N_WARMUP)
N_PARAMS    <- length(param_names_log) # number of parameters being sampled (log scale)
EPS_INIT    <- 0.01    # initial step size (dual averaging will adapt)
L_STEPS     <- 10L     # leapfrog steps per proposal
ADAPT_DELTA <- 0.65    # target acceptance rate

# ── Initial points ────────────────────────────────────────────────────────────
set.seed(114514)
init_center <- data$theta_prior$center_smooth
inits <- lapply(seq_len(N_CHAINS), function(ch) {
  rnorm(2L * N_AGE, init_center, 0.05)
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
param_corr_log   <- cor(all_samples)

cat(sprintf("\nTotal post-warmup samples: %d (%d × %d chains)\n",
            nrow(all_samples), nrow(post_warmup_list[[1L]]), N_CHAINS))

cat("\n=== Parameter Correlation Matrix (log scale / fitted theta) ===\n")
print(round(param_corr_log, 3))

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

plot_dir <- file.path(output_dir, "hmc_plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(plot_dir, "trace_plot.png"), plot = p_trace)
ggsave(file.path(plot_dir, "density_plot.png"), plot = p_density)
ggsave(file.path(plot_dir, "hist_prev_plot.png"), plot = p_hist_prev)
ggsave(file.path(plot_dir, "hist_tot_plot.png"), plot = p_hist_tot)
ggsave(file.path(plot_dir, "interval_plot.png"), plot = p_interval)
ggsave(file.path(plot_dir, "prev_int_plot.png"), plot = p_prev_int)

# ── Save ─────────────────────────────────────────────────────────────────────
saveRDS(
  list(
    hmc_chains       = hmc_chains,
    post_warmup_list = post_warmup_list,
    all_samples      = all_samples,
    param_corr_log   = param_corr_log,
    diag_table       = diag_table,
    post_summary     = post_summary,
    ppc_out          = ppc_out,
    theta_prior      = data$theta_prior,
    nelder_mead_fit_file = nm_fit_file
  ),
  file = file.path(output_dir, "hmc_output_full.rds")
)
cat(sprintf("\nAll results saved to %s\n", file.path(output_dir, "hmc_output_full.rds")))
result = readRDS(file.path(output_dir, "hmc_output_full.rds"))
hmc_chains = result$hmc_chains
post_warmup_list = result$post_warmup_list
all_samples = result$all_samples
param_corr_log = result$param_corr_log
diag_table = result$diag_table
post_summary = result$post_summary
ppc_out = result$ppc_out
