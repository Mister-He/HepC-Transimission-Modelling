# =============================================================================
# MAIN: RUN HMC
#
# Default two-step use:
#   Rscript hmc.R --nm-fit=two-steps-calibration/nelder_mead_fit.rds \
#                 --out-dir=two-steps-calibration
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# Observations. Keep identical to nelder_mead.R: prevalence likelihood is
# obs_pos[a] ~ Binomial(obs_tot[a], q_age[a]).
obs_prev <- c(
  0.1118421, 0.1731044, 0.2684954, 0.4301165, 0.4821029, 0.3544304
)
obs_prev_se <- sqrt(c(
  0.09933344, 0.14313927, 0.19640562, 0.24511630, 0.24967969, 0.22880949
))
obs_tot <- c(99, 1244, 1467, 1841, 1628, 409)
obs_pos <- round(obs_prev * obs_tot)
N_AGE <- length(obs_tot)
N_CONTACT <- N_AGE

param_names_orig <- c(
  paste0("C_contact_scale_", seq_len(N_CONTACT)),
  paste0("inflow_scale_", seq_len(N_CONTACT))
)

source("setup.R")
source("HMC_core.r")
source("HMC_conv.R")

# Shared observation-model settings used by the NM stage. prev_logit_sd is only
# used to draw observed uncertainty bands; the fitted prevalence likelihood is
# binomial in HMC_core.r.
data$prev_logit_sd <- 0.10
data$sigma_pop <- c(0.10, rep(0.06, N_AGE - 1L))

parse_cli <- function(args) {
  get_arg <- function(name, default = NULL) {
    hit <- grep(paste0("^--", name, "="), args, value = TRUE)
    if (length(hit) == 0L) return(default)
    sub(paste0("^--", name, "="), "", hit[[1L]])
  }
  list(
    nm_fit = get_arg("nm-fit", NULL),
    out_dir = get_arg("out-dir", "two-steps-calibration"),
    chains = as.integer(get_arg("chains", "4")),
    cores = as.integer(get_arg("cores", "4")),
    warmup = as.integer(get_arg("warmup", "500")),
    iter = as.integer(get_arg("iter", "600")),
    ppc = as.integer(get_arg("ppc", "600"))
  )
}

configure_hmc_prior_from_nm <- function(nm_fit_path, prior_sd = 0.45) {
  if (is.null(nm_fit_path) || !nzchar(nm_fit_path)) return(NULL)
  if (!file.exists(nm_fit_path)) stop(sprintf("Nelder-Mead fit not found: %s", nm_fit_path))

  nm_fit <- readRDS(nm_fit_path)
  theta_raw <- nm_fit$theta_hat
  theta_hat <- as.numeric(theta_raw)
  names_match <- !is.null(names(theta_raw)) && identical(names(theta_raw), param_names_log)
  if (length(theta_hat) != 2L * N_AGE || any(!is.finite(theta_hat)) || !names_match) {
    stop(sprintf(
      paste0(
        "Nelder-Mead theta_hat is incompatible with the current %d-parameter ",
        "direct age-specific scaling-factor model; refit with nelder_mead.R"
      ),
      2L * N_AGE
    ))
  }

  configure_age_priors(
    contact_mean = theta_hat[seq_len(N_AGE)],
    inflow_scale_mean = theta_hat[N_AGE + seq_len(N_AGE)],
    contact_sd = prior_sd,
    inflow_scale_sd = prior_sd
  )
  nm_fit
}

predict_at_theta <- function(theta, base_params = params, sim_data = data) {
  pm <- build_params_from_theta(theta, base_params)
  out <- run_sim(pm, sim_data)
  if (!is.matrix(out) || nrow(out) < 1L || any(!is.finite(out))) {
    stop("Model simulation failed at the supplied theta")
  }
  ans <- compute_age_quantities(as.numeric(out[nrow(out), -1L]))
  if (is.null(ans)) stop("Simulation produced invalid age-specific quantities")
  ans
}

run_hmc_calibration <- function(nm_fit_path = NULL,
                                out_dir = "two-steps-calibration",
                                n_chains = 4L,
                                n_cores = 4L,
                                n_warmup = 500L,
                                n_iter = 600L,
                                n_ppc = 600L) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  plot_dir <- file.path(out_dir, "hmc_plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  nm_fit <- configure_hmc_prior_from_nm(nm_fit_path)
  N_PARAMS <- length(param_names_log)
  EPS_INIT <- 0.01
  L_STEPS <- 10L
  ADAPT_DELTA <- 0.65

  set.seed(114514)
  if (!is.null(nm_fit)) {
    theta_center <- as.numeric(nm_fit$theta_hat)
    inits <- lapply(seq_len(n_chains), function(ch) {
      theta_center + rnorm(length(theta_center), 0, 0.10)
    })
  } else {
    inits <- lapply(seq_len(n_chains), function(ch) {
      c(
        rnorm(N_AGE, CONTACT_PRIOR_MEANS, 0.25),
        rnorm(N_AGE, INFLOW_SCALE_PRIOR_MEANS, 0.25)
      )
    })
  }

  cat("=== HMC Calibration: HCV PWID Model ===\n")
  cat(sprintf("Chains: %d | Warmup: %d | Sampling: %d | L: %d\n",
              n_chains, n_warmup, n_iter - n_warmup, L_STEPS))
  if (!is.null(nm_fit_path)) {
    cat(sprintf("Using Nelder-Mead priors from: %s\n", nm_fit_path))
  }
  cat(sprintf("Approx ODE runs per HMC step: %d\n", (L_STEPS + 1L) * 2L * N_PARAMS))

  hmc_chains <- parallel::mclapply(seq_len(n_chains), function(ch) {
    run_hmc_chain(
      theta_init = inits[[ch]],
      n_iter = n_iter,
      n_warmup = n_warmup,
      eps_init = EPS_INIT,
      L = L_STEPS,
      adapt_delta = ADAPT_DELTA,
      base_params = params,
      data = data,
      seed = ch * 314L,
      chain_id = ch
    )
  }, mc.cores = min(n_cores, n_chains))

  post_warmup_list <- lapply(hmc_chains, function(ch) ch$samples)
  all_samples <- do.call(rbind, post_warmup_list)
  param_corr_log <- cor(all_samples)

  cat(sprintf("\nTotal post-warmup samples: %d (%d x %d chains)\n",
              nrow(all_samples), nrow(post_warmup_list[[1L]]), n_chains))
  cat("\n=== Parameter Correlation Matrix (log scale / fitted theta) ===\n")
  print(round(param_corr_log, 3))

  diag_table <- print_diagnostics(post_warmup_list, chains_raw = hmc_chains)

  cat("\nAcceptance rates per chain (post-warmup):\n")
  acc_df <- data.frame(
    Chain = seq_len(n_chains),
    Acc_rate = round(vapply(hmc_chains, function(ch) ch$acc_rate, numeric(1L)), 3),
    Eps_final = round(vapply(hmc_chains, function(ch) ch$eps_final, numeric(1L)), 5)
  )
  print(acc_df, row.names = FALSE)

  orig_samples <- theta_to_orig(all_samples)
  post_summary <- data.frame(
    Parameter = param_names_orig,
    Mean = round(colMeans(orig_samples), 4),
    SD = round(apply(orig_samples, 2, sd), 4),
    Q2.5 = round(apply(orig_samples, 2, quantile, 0.025), 4),
    Q50 = round(apply(orig_samples, 2, quantile, 0.500), 4),
    Q97.5 = round(apply(orig_samples, 2, quantile, 0.975), 4)
  )
  cat("\n=== Posterior Summary (original scale) ===\n")
  print(post_summary, row.names = FALSE)

  set.seed(114514)
  cat("\n=== Generating Posterior Predictive Samples ===\n")
  ppc_out <- generate_ppc_samples(all_samples, params, data, n_ppc = n_ppc)

  ppp_df <- data.frame(
    Age = paste0("Age ", seq_len(N_AGE)),
    ppp_prev = round(ppc_out$ppp_prev, 3),
    ppp_tot = round(ppc_out$ppp_tot, 3)
  )
  cat("\nPosterior predictive p-values (near 0.5 = good fit):\n")
  print(ppp_df, row.names = FALSE)

  p_trace <- plot_traces(hmc_chains, param_idx = 1:N_PARAMS)
  p_density <- plot_posterior_densities(post_warmup_list)
  p_hist_prev <- plot_ppc_histograms(ppc_out, type = "prev")
  p_hist_tot <- plot_ppc_histograms(ppc_out, type = "tot")
  p_interval <- plot_ppc_intervals(ppc_out)
  p_prev_int <- plot_ppc_prevalence_intervals(ppc_out)
  p_smooth <- plot_posterior_fitted_curves(ppc_out)

  ggsave(file.path(plot_dir, "trace_plot.png"), plot = p_trace, width = 8, height = 8, dpi = 300)
  ggsave(file.path(plot_dir, "density_plot.png"), plot = p_density, width = 8, height = 10, dpi = 300)
  ggsave(file.path(plot_dir, "hist_prev_plot.png"), plot = p_hist_prev, width = 8, height = 8, dpi = 300)
  ggsave(file.path(plot_dir, "hist_tot_plot.png"), plot = p_hist_tot, width = 8, height = 8, dpi = 300)
  ggsave(file.path(plot_dir, "interval_plot.png"), plot = p_interval, width = 8, height = 7, dpi = 300)
  ggsave(file.path(plot_dir, "prev_int_plot.png"), plot = p_prev_int, width = 8, height = 5, dpi = 300)
  ggsave(file.path(plot_dir, "posterior_fitted_curves.png"), plot = p_smooth, width = 8, height = 7, dpi = 300)

  output <- list(
    nm_fit_path = nm_fit_path,
    nm_fit = nm_fit,
    hmc_chains = hmc_chains,
    post_warmup_list = post_warmup_list,
    all_samples = all_samples,
    param_corr_log = param_corr_log,
    diag_table = diag_table,
    post_summary = post_summary,
    ppc_out = ppc_out,
    ppp_df = ppp_df,
    age_specific_prior = list(
      contact_mean = CONTACT_PRIOR_MEANS,
      inflow_scale_mean = INFLOW_SCALE_PRIOR_MEANS,
      contact_sd = CONTACT_PRIOR_SDS,
      inflow_scale_sd = INFLOW_SCALE_PRIOR_SDS
    )
  )
  saveRDS(output, file = file.path(out_dir, "hmc_output_full.rds"))
  cat(sprintf("\nAll HMC results saved to %s\n", file.path(out_dir, "hmc_output_full.rds")))
  invisible(output)
}

if (sys.nframe() == 0L) {
  cli <- parse_cli(commandArgs(trailingOnly = TRUE))
  run_hmc_calibration(
    nm_fit_path = cli$nm_fit,
    out_dir = cli$out_dir,
    n_chains = cli$chains,
    n_cores = cli$cores,
    n_warmup = cli$warmup,
    n_iter = cli$iter,
    n_ppc = cli$ppc
  )
}
