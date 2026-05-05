# =============================================================================
# setup_stan.R  —  HCV PWID model calibration via cmdstanr
# =============================================================================
# Uses cmdstanr instead of rstan to avoid arm64/llvm/sink() issues.
#
# One-time setup (run once, then never again):
#   install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/",
#                                           getOption("repos")))
#   cmdstanr::install_cmdstan()   # installs standalone Stan binary (~200 MB)
# =============================================================================

library(Rcpp)
library(cmdstanr) # replaces rstan
library(posterior) # tidy posterior extraction (installed with cmdstanr)
library(bayesplot) # mcmc_trace, ppc plots
library(dplyr)
library(ggplot2)

sourceCpp("sim.cpp")

# ── Check cmdstan is installed and print version ──────────────────────────────
cmdstan_version() # should print e.g. "2.35.0"

# =============================================================================
# HELPER
# =============================================================================
idx <- function(s, k, h, i) s * 4L * 4L * 9L + (k - 1L) * 4L * 9L + h * 9L + i + 1L

# =============================================================================
# FIXED PARAMETERS (unchanged from setup.R)
# =============================================================================
params <- list(
  q = 0.0057,
  kappa = 0.26,
  iota1 = 180 / 365,
  iota2 = 84 / 365,
  rho = 0.78,
  alpha_NC = 0.9787,
  alpha_DC_pos = 0.8818,
  alpha_DC_neg = 0.5490,
  alpha_HCC = 0.8818,
  tau = c(0.0, 0.0, 0.0, 0.0),
  p_NC_CC = 0.027,
  p_CC_DC = 0.039,
  p_CC_HCC = 0.014,
  p_DC_HCC = 0.014,
  r3_NC_CC = 1.36,
  r3_CC_DC = 1.36,
  r3_CC_HCC = 1.93,
  r3_DC_HCC = 1.93,
  phi_CC_DC = 0.07,
  phi_CC_HCC = 0.23,
  phi_DC_HCC = 1.00,
  mu = c(
    0.001267, 0.000300, 0.000300, 0.000400, 0.000500,
    0.000700, 0.001400, 0.002300, 0.016100
  ),
  omega = 14.68,
  mu_DC = 0.130,
  mu_HCC = 0.430,
  psi_DC = 0.45,
  psi_HCC = 0.37,
  lambda1 = c(
    0.3985248, 0.5686156, 0.5485466,
    0.6830963, 0.6971583, 1.1823825,
    1.6449108, 1.2289346, 0.6400324
  ),
  lambda2 = c(
    0.489, 0.620, 0.663, 0.628, 0.533,
    0.475, 0.472, 0.441, 0.451
  ),
  lambda3 = c(
    0.3985248, 0.5686156, 0.5485466,
    0.6830963, 0.6971583, 1.1823825,
    1.6449108, 1.2289346, 0.6400324
  ),
  pi_recid = 0.65,
  C_contact = rbind(
    c(7, 4, 1, 1, 0, 1, 1, 0, 0),
    c(11, 34, 21, 11, 6, 2, 1, 0.5, 0.5),
    c(7, 30, 80, 62, 30, 10, 2, 1.5, 1.5),
    c(2, 10, 60, 121, 65, 38, 15, 2, 2),
    c(1, 11, 22, 67, 107, 41, 18, 2.5, 2.5),
    c(0, 4, 6, 22, 32, 31, 10, 2, 2),
    c(0, 1, 1, 8, 10, 15, 11, 0.5, 0.5),
    c(0, 5, 5, 5.5, 6.5, 6.5, 3.5, 1.5, 1.5),
    c(0, 5, 5, 5.5, 6.5, 6.5, 3.5, 1.5, 1.5)
  ) * 4,
  beta = c(235, 565 / 2, 565 / 2, 301 / 2, 301 / 2, 111 / 2, 111 / 2, 33 / 2, 33 / 2 + 4)
)

# =============================================================================
# OBSERVED CALIBRATION TARGETS
# =============================================================================
obs_J_total <- c(307L, 797L, 829L, 633L, 598L, 642L, 481L, 439L, 366L)
obs_J_susc <- c(55L, 145L, 183L, 164L, 212L, 299L, 222L, 190L, 133L)

# =============================================================================
# BASE INITIAL CONDITIONS
# =============================================================================
y0_base <- rep(0.0, 576)
age_wt <- params$beta / sum(params$beta) * 1000.0
for (i in 0:8) {
  y0_base[idx(0, 1, 0, i)] <- age_wt[i + 1]
}

# =============================================================================
# STAN DATA LIST
# =============================================================================
stan_data <- list(
  obs_J_total    = obs_J_total,
  obs_J_susc     = obs_J_susc,
  n_years        = 30L,
  steps_per_year = 12L,
  y0_base        = y0_base,
  q              = params$q,
  kappa          = params$kappa,
  iota1          = params$iota1,
  iota2          = params$iota2,
  rho            = params$rho, # needed by transformed data

  alpha_NC       = params$alpha_NC,
  alpha_DC_pos   = params$alpha_DC_pos,
  alpha_DC_neg   = params$alpha_DC_neg,
  alpha_HCC      = params$alpha_HCC,
  tau            = params$tau,
  p_NC_CC        = params$p_NC_CC,
  p_CC_DC        = params$p_CC_DC,
  p_CC_HCC       = params$p_CC_HCC,
  p_DC_HCC       = params$p_DC_HCC,
  r3_NC_CC       = params$r3_NC_CC,
  r3_CC_DC       = params$r3_CC_DC,
  r3_CC_HCC      = params$r3_CC_HCC,
  r3_DC_HCC      = params$r3_DC_HCC,
  phi_CC_DC      = params$phi_CC_DC,
  phi_CC_HCC     = params$phi_CC_HCC,
  phi_DC_HCC     = params$phi_DC_HCC,
  mu             = params$mu,
  omega          = params$omega,
  mu_DC          = params$mu_DC,
  mu_HCC         = params$mu_HCC,
  psi_DC         = params$psi_DC,
  psi_HCC        = params$psi_HCC,
  lambda1        = params$lambda1,
  lambda2        = params$lambda2,
  lambda3        = params$lambda3,
  pi_recid       = params$pi_recid,
  C_contact      = params$C_contact,
  beta_base      = params$beta
)

# =============================================================================
# COMPILE STAN MODEL
# cmdstanr:
#   - prints the FULL error with line numbers if compilation fails
#   - uses a standalone stanc binary (no rstan/Rcpp/llvm conflicts)
#   - compiled model is cached to <stan_file>.exe automatically
# =============================================================================
message("Compiling Stan model (first run only, then cached)...")
mod <- cmdstan_model(
  stan_file = "sim.stan",
  compile   = TRUE
)

# ── Optional: syntax check before sampling (prints exact line of error) ───────
# mod$check_syntax(pedantic = TRUE)

# =============================================================================
# INITIAL VALUES
# One list per chain, both parameters start near 1 with small jitter.
# =============================================================================
init_list <- lapply(seq_len(4), function(chain_id) {
  list(
    beta_scale = runif(1, 0, 1),
    y0_scale   = runif(1, 0, 10)
  )
})

# =============================================================================
# RUN HMC SAMPLER
# =============================================================================
message("Running HMC sampler...")
t0 <- proc.time()

fit <- mod$sample(
  data = stan_data,
  init = init_list,
  seed = 114514,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.95,
  max_treedepth = 12,
  refresh = 100,

  # ── persistence: keep CSV files alongside the saved fit object ─────────────
  output_dir = "stan_output", # folder to save CSVs (created if missing)
  output_basename = "hcv_fit", # prefix for CSV filenames
  save_cmdstan_config = TRUE # saves the config JSON needed for reload
)

# Save the R object (now points to files in stan_output/)
fit$save_object("hcv_cmdstan_fit.rds")

elapsed <- (proc.time() - t0)["elapsed"]
message(sprintf("Sampling completed in %.1f minutes.", elapsed / 60))

# =============================================================================
# DIAGNOSTICS
# =============================================================================

cat("\n=== POSTERIOR SUMMARY ===\n")
print(fit$summary(variables = c("beta_scale", "y0_scale")))

diag <- fit$diagnostic_summary()
cat(sprintf(
  "\nDivergent transitions: %d  |  Max treedepth exceeded: %d\n",
  sum(diag$num_divergent), sum(diag$num_max_treedepth)
))
if (sum(diag$num_divergent) > 0) {
  warning("Divergences detected — increase adapt_delta (e.g. 0.99).")
}

# =============================================================================
# EXTRACT POSTERIOR DRAWS
# =============================================================================
post_beta <- as_draws_matrix(fit$draws("beta_scale"))[, 1]
post_y0 <- as_draws_matrix(fit$draws("y0_scale"))[, 1]

cat(sprintf(
  "\nbeta_scale: median = %.4f  95%% CI [%.4f, %.4f]\n",
  median(post_beta), quantile(post_beta, 0.025), quantile(post_beta, 0.975)
))
cat(sprintf(
  "y0_scale:   median = %.4f  95%% CI [%.4f, %.4f]\n",
  median(post_y0), quantile(post_y0, 0.025), quantile(post_y0, 0.975)
))

# =============================================================================
# PLOTS
# =============================================================================
age_labels <- paste0("age", 1:9)

# ── Trace plots ──────────────────────────────────────────────────────────────
p_trace <- mcmc_trace(fit$draws(c("beta_scale", "y0_scale"))) +
  theme_minimal(base_size = 13)
print(p_trace)

# ── Joint posterior ──────────────────────────────────────────────────────────
p_joint <- ggplot(
  data.frame(beta_scale = post_beta, y0_scale = post_y0),
  aes(x = beta_scale, y = y0_scale)
) +
  geom_point(alpha = 0.12, colour = "steelblue", size = 0.7) +
  geom_density_2d(colour = "navy", linewidth = 0.5) +
  labs(
    title = "Joint posterior", x = expression(beta[scale]),
    y = expression(y[0][",scale"])
  ) +
  theme_minimal(base_size = 13)
print(p_joint)

# ── PPC helper ───────────────────────────────────────────────────────────────
ppc_plot <- function(var_prefix, obs_vec, title_str) {
  ppc_mat <- as_draws_matrix(
    fit$draws(paste0(var_prefix, "[", 1:9, "]"))
  )
  df <- data.frame(
    age      = rep(age_labels, each = nrow(ppc_mat)),
    ppc      = as.vector(ppc_mat),
    observed = rep(obs_vec, each = nrow(ppc_mat))
  )
  ggplot(df, aes(x = age)) +
    stat_summary(aes(y = ppc),
      fun = median,
      geom = "point", colour = "steelblue", size = 3
    ) +
    stat_summary(aes(y = ppc),
      fun.min = \(x) quantile(x, 0.025),
      fun.max = \(x) quantile(x, 0.975),
      geom = "errorbar", colour = "steelblue", width = 0.3
    ) +
    geom_point(aes(y = observed), colour = "firebrick", size = 3.5, shape = 17) +
    labs(
      title = title_str,
      subtitle = "Blue: median + 95% PPI  |  Red \u25b2: observed",
      x = "Age group", y = "Count"
    ) +
    theme_minimal(base_size = 13)
}

print(ppc_plot("ppc_J_total", obs_J_total, "PPC \u2014 Total J stratum by age"))
print(ppc_plot("ppc_J_susc", obs_J_susc, "PPC \u2014 Susceptible J stratum by age"))

# ── Residuals ────────────────────────────────────────────────────────────────
resid_tot_mat <- as_draws_matrix(fit$draws(paste0("resid_J_total[", 1:9, "]")))
resid_susc_mat <- as_draws_matrix(fit$draws(paste0("resid_J_susc[", 1:9, "]")))

resid_df <- data.frame(
  age = rep(age_labels, 2),
  target = rep(c("J total", "J susceptible"), each = 9),
  mean_res = c(colMeans(resid_tot_mat), colMeans(resid_susc_mat)),
  lo95 = c(
    apply(resid_tot_mat, 2, quantile, 0.025),
    apply(resid_susc_mat, 2, quantile, 0.025)
  ),
  hi95 = c(
    apply(resid_tot_mat, 2, quantile, 0.975),
    apply(resid_susc_mat, 2, quantile, 0.975)
  )
)

p_resid <- ggplot(resid_df, aes(x = age, y = mean_res, colour = target)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3, position = position_dodge(0.4)) +
  geom_errorbar(aes(ymin = lo95, ymax = hi95),
    width = 0.3, position = position_dodge(0.4)
  ) +
  scale_colour_manual(values = c(
    "J total" = "steelblue",
    "J susceptible" = "coral"
  )) +
  labs(
    title = "Posterior residuals (observed \u2212 predicted)",
    x = "Age group", y = "Residual", colour = NULL
  ) +
  theme_minimal(base_size = 13)
print(p_resid)

# =============================================================================
# SAVE
# =============================================================================
fit$save_object("hcv_cmdstan_fit.rds")
message("Fit saved to hcv_cmdstan_fit.rds")
# Reload: fit <- readRDS("hcv_cmdstan_fit.rds")

# LOO
log_lik_all <- cbind(
  as_draws_matrix(fit$draws(paste0("log_lik_total[", 1:9, "]"))),
  as_draws_matrix(fit$draws(paste0("log_lik_susc[", 1:9, "]")))
)
if (requireNamespace("loo", quietly = TRUE)) {
  print(loo::loo(log_lik_all))
}

# =============================================================================
# VALIDATE AGAINST C++ SIMULATOR
# =============================================================================
beta_scale_hat <- median(post_beta)
y0_scale_hat <- median(post_y0)

y0_cal <- rep(0.0, 576)
for (i in 0:8) {
  y0_cal[idx(0, 1, 0, i)] <- y0_scale_hat * age_wt[i + 1]
}

out_cal <- run_sim(
  modifyList(params, list(beta = params$beta * beta_scale_hat)),
  list(t_start = 0.0, t_end = 30.0, dt = 1 / 12, y0 = y0_cal)
)

final <- nrow(out_cal)
J_total_cpp <- vapply(0:8, function(i) {
  sum(vapply(1:4, function(k) {
    sum(vapply(0:3, function(h) {
      out_cal[final, idx(1, k, h, i) + 1L]
    }, numeric(1)))
  }, numeric(1)))
}, numeric(1))
J_susc_cpp <- vapply(0:8, function(i) {
  sum(vapply(1:4, function(k) {
    out_cal[final, idx(1, k, 0, i) + 1L]
  }, numeric(1)))
}, numeric(1))

cat("\n=== C++ VALIDATION ===\n")
print(
  data.frame(
    age = age_labels,
    obs_tot = obs_J_total, pred_tot = round(J_total_cpp, 1),
    obs_susc = obs_J_susc, pred_susc = round(J_susc_cpp, 1)
  ),
  row.names = FALSE
)

# =============================================================================
# IF COMPILATION STILL FAILS — full diagnostics
# =============================================================================
# 1. Get the exact Stan error with line number:
#      mod$check_syntax(pedantic = TRUE)
#
# 2. Run stanc directly from shell:
#      system(paste0(cmdstanr::cmdstan_path(), "/bin/stanc hcv_model.stan"))
#
# 3. If sticking with rstan, capture the full error:
#      tryCatch(
#        rstan::stan_model("hcv_model.stan", verbose = TRUE),
#        error = function(e) writeLines(conditionMessage(e), "err.txt")
#      )