# =============================================================================
# nc_diagnostics.R  —  Convergence diagnostics and plots for nested calibration
#
# Run AFTER nc_run.R has completed.  Source from the repo root:
#   source("nested_calibration/nc_diagnostics.R")
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(coda)   # gelman.diag, effectiveSize

source("nested_calibration/nc_core.R")  # constants (param names, obs_*)

OUT_DIR <- "nested_calibration/output"
results_s1 <- readRDS(file.path(OUT_DIR, "nc_stage1_results.rds"))
results_s2 <- readRDS(file.path(OUT_DIR, "nc_stage2_results.rds"))
combined   <- readRDS(file.path(OUT_DIR, "nc_combined_results.rds"))

# ── 95% Wilson CIs for observed data points ───────────────────────────────────
wilson_ci <- function(x, n, alpha = 0.05) {
  z     <- qnorm(1 - alpha / 2)
  p     <- x / n
  denom <- 1 + z^2 / n
  ctr   <- (p + z^2 / (2 * n)) / denom
  half  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  list(lo = pmax(0, ctr - half), hi = pmin(1, ctr + half))
}

N_total      <- sum(obs_tot)
prev_ci95    <- wilson_ci(obs_pos, obs_tot)      # per-age prevalence CIs
prop_ci95    <- wilson_ci(obs_tot, N_total)      # per-age proportion CIs

age_labels <- paste0("Age ", 1:9)

# =============================================================================
# 1.  CONVERGENCE DIAGNOSTICS
# =============================================================================

summarise_convergence <- function(chains, stage_label) {
  coda_chains <- mcmc.list(lapply(chains, function(ch) mcmc(ch$samples)))
  rhat <- gelman.diag(coda_chains, multivariate = FALSE)
  ess  <- effectiveSize(coda_chains)

  cat(sprintf("\n=== %s Convergence ===\n", stage_label))
  cat("Gelman-Rubin R-hat (target < 1.05 for all params):\n")
  print(round(rhat$psrf[, 1, drop = FALSE], 3))
  cat("Effective sample size:\n")
  print(round(ess, 0))

  invisible(list(rhat = rhat, ess = ess))
}

diag_s1 <- summarise_convergence(results_s1$chains, "Stage 1 — contact rates")
diag_s2 <- summarise_convergence(results_s2$chains, "Stage 2 — beta deviations")

# =============================================================================
# 2.  TRACE PLOTS
# =============================================================================

make_trace_df <- function(chains) {
  bind_rows(lapply(seq_along(chains), function(k) {
    ch <- chains[[k]]
    as.data.frame(ch$samples_all) %>%
      mutate(
        iter   = seq_len(nrow(.)),
        chain  = factor(k),
        warmup = iter <= ch$n_warmup
      )
  }))
}

plot_traces <- function(chains, title) {
  df   <- make_trace_df(chains)
  long <- pivot_longer(df, -c(iter, chain, warmup), names_to = "param", values_to = "value")

  ggplot(long, aes(x = iter, y = value, colour = chain)) +
    geom_line(data = filter(long, warmup),  linewidth = 0.25, alpha = 0.25) +
    geom_line(data = filter(long, !warmup), linewidth = 0.25, alpha = 0.80) +
    geom_vline(
      data = distinct(long, param) %>% mutate(warmup_end = chains[[1]]$n_warmup),
      aes(xintercept = warmup_end), linetype = "dashed", colour = "grey50", linewidth = 0.4
    ) +
    facet_wrap(~param, scales = "free_y", ncol = 2) +
    labs(title = title, x = "Iteration", y = "Value", colour = "Chain") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(size = 8))
}

p_trace_s1 <- plot_traces(results_s1$chains, "Stage 1 trace plots — contact-rate parameters")
p_trace_s2 <- plot_traces(results_s2$chains, "Stage 2 trace plots — log_beta_delta parameters")

# =============================================================================
# 3.  POSTERIOR DENSITY PLOTS
# =============================================================================

plot_posterior_densities <- function(samples, title) {
  pivot_longer(as.data.frame(samples), everything(),
               names_to = "param", values_to = "value") %>%
    ggplot(aes(x = value)) +
    geom_density(fill = "#2166AC", alpha = 0.35, colour = "#2166AC", linewidth = 0.6) +
    facet_wrap(~param, scales = "free", ncol = 2) +
    labs(title = title, x = NULL, y = "Density") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(size = 8))
}

p_dens_s1 <- plot_posterior_densities(
  results_s1$samples,
  "Stage 1 posterior densities — contact-rate parameters"
)

# For Stage 2, plot on exp scale (beta scale factors) for interpretability
beta_scale_samples <- exp(results_s2$samples)
colnames(beta_scale_samples) <- paste0("beta_scale_", 1:9)
p_dens_s2 <- plot_posterior_densities(
  beta_scale_samples,
  "Stage 2 posterior densities — beta scale factors (exp of log_beta_delta)"
)

# Implied contact scales (original scale) — exclude reference row (constant = 1.0)
contact_orig     <- contact_to_orig(results_s1$samples)
keep_cols        <- grepl("C_contact_scale", colnames(contact_orig)) &
                    colnames(contact_orig) != "C_contact_scale_1"  # row 1 is fixed reference
contact_scale_df <- as.data.frame(contact_orig[, keep_cols, drop = FALSE])
p_dens_contact_scales <- plot_posterior_densities(
  contact_scale_df,
  "Contact-rate scale posteriors (original scale; row 1 = 1.0 reference, excluded)"
)

# =============================================================================
# 4.  POSTERIOR PREDICTIVE CHECK (PPC)
# =============================================================================

ppc_df <- combined$ppc_df

obs_ci_prop_df <- data.frame(
  age_group = paste0("Age ", 1:9),
  obs_prop  = obs_prop,
  lo95      = prop_ci95$lo,
  hi95      = prop_ci95$hi
)
obs_ci_prev_df <- data.frame(
  age_group = paste0("Age ", 1:9),
  obs_prev  = obs_prev,
  lo95      = prev_ci95$lo,
  hi95      = prev_ci95$hi
)

p_ppc_prop <- ggplot(ppc_df, aes(x = age_group)) +
  geom_col(aes(y = model_prop), fill = "#4292C6", alpha = 0.75, width = 0.55) +
  geom_errorbar(data = obs_ci_prop_df, aes(ymin = lo95, ymax = hi95),
                width = 0.22, colour = "#D7191C", linewidth = 0.9) +
  geom_point(aes(y = obs_prop),
             shape = 21, size = 3.5, fill = "#D7191C", colour = "black", stroke = 1.2) +
  scale_y_continuous(labels = function(x) sprintf("%.1f%%", 100 * x)) +
  labs(
    title    = "PPC: age proportions — bar = model posterior mean, dot = observed",
    subtitle = sprintf("Mean absolute error = %.4f | Error bars: 95%% CI on observed", mean(ppc_df$prop_err)),
    x = "Age group", y = "Proportion of PWID"
  ) +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())

p_ppc_prev <- ggplot(ppc_df, aes(x = age_group)) +
  geom_col(aes(y = model_prev), fill = "#4292C6", alpha = 0.75, width = 0.55) +
  geom_errorbar(data = obs_ci_prev_df, aes(ymin = lo95, ymax = hi95),
                width = 0.22, colour = "#D7191C", linewidth = 0.9) +
  geom_point(aes(y = obs_prev),
             shape = 21, size = 3.5, fill = "#D7191C", colour = "black", stroke = 1.2) +
  scale_y_continuous(labels = function(x) sprintf("%.1f%%", 100 * x)) +
  labs(
    title    = "PPC: HCV prevalence — bar = model posterior mean, dot = observed",
    subtitle = sprintf("Mean absolute error = %.4f | Error bars: 95%% CI on observed", mean(ppc_df$prev_err)),
    x = "Age group", y = "HCV prevalence"
  ) +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())

# =============================================================================
# 5.  CREDIBLE INTERVAL BANDS FROM POSTERIOR SAMPLES
# =============================================================================
# Draw a random subset of posterior samples, run the ODE for each, and
# compute model-implied proportions and prevalences to build Bayesian CIs.

source("setup.R")       # params, data, run_sim
source("HMC_core.r")    # compute_age_quantities

N_PPC_DRAWS <- min(200L, nrow(results_s1$samples))
set.seed(42L)
ppc_idx <- sample(nrow(results_s1$samples), N_PPC_DRAWS)

ppc_prop_mat <- matrix(NA_real_, N_PPC_DRAWS, 9)
ppc_prev_mat <- matrix(NA_real_, N_PPC_DRAWS, 9)

cat(sprintf("Running %d PPC simulations from posterior samples...\n", N_PPC_DRAWS))
beta_nested_ppc <- combined$beta_nested   # Stage-1 nested beta; Stage-2 theta_b is relative to this

for (ii in seq_len(N_PPC_DRAWS)) {
  # Stage 1 contact params
  theta_c <- results_s1$samples[ppc_idx[ii], ]
  # Stage 2 beta perturbation (cycling if sizes differ)
  theta_b <- results_s2$samples[((ppc_idx[ii] - 1L) %% nrow(results_s2$samples)) + 1L, ]

  pm      <- build_contact_params(theta_c, params)
  pm$beta <- beta_nested_ppc * exp(theta_b)

  # Stage 2 now uses nested adjustment internally — apply it here too so the
  # PPC samples are drawn from the same model as the Stage 2 posterior.
  result <- tryCatch(nested_beta_adjust(pm, data, obs_prop), error = function(e) NULL)
  if (is.null(result)) next

  ppc_prop_mat[ii, ] <- result$model_q$p_age   # enforced ≈ obs_prop by construction
  ppc_prev_mat[ii, ] <- result$model_q$q_age
}

# Summarise into credible intervals (includes 95% Wilson CIs on observations)
ci_summary <- function(mat, obs_vals, obs_lo, obs_hi, label) {
  data.frame(
    Age      = factor(age_labels, levels = age_labels),
    obs      = obs_vals,
    obs_lo95 = obs_lo,
    obs_hi95 = obs_hi,
    med      = apply(mat, 2, median,   na.rm = TRUE),
    lo50     = apply(mat, 2, quantile, 0.25,  na.rm = TRUE),
    hi50     = apply(mat, 2, quantile, 0.75,  na.rm = TRUE),
    lo95     = apply(mat, 2, quantile, 0.025, na.rm = TRUE),
    hi95     = apply(mat, 2, quantile, 0.975, na.rm = TRUE),
    label    = label,
    stringsAsFactors = FALSE
  )
}

ci_prop <- ci_summary(ppc_prop_mat, obs_prop, prop_ci95$lo, prop_ci95$hi, "Age proportion")
ci_prev <- ci_summary(ppc_prev_mat, obs_prev, prev_ci95$lo, prev_ci95$hi, "HCV prevalence")

plot_ci_bands <- function(df, y_lab) {
  ggplot(df, aes(x = Age)) +
    geom_linerange(aes(ymin = lo95, ymax = hi95),
                   linewidth = 1.0, colour = "#4292C6", alpha = 0.50) +
    geom_linerange(aes(ymin = lo50, ymax = hi50),
                   linewidth = 3.0, colour = "#08519C", alpha = 0.80) +
    geom_point(aes(y = med),
               shape = 18, size = 3.5, colour = "#08306B") +
    geom_errorbar(aes(ymin = obs_lo95, ymax = obs_hi95),
                  width = 0.22, colour = "#D7191C", linewidth = 0.9) +
    geom_point(aes(y = obs),
               shape = 21, size = 3.5, fill = "#D7191C", colour = "black", stroke = 1.2) +
    scale_y_continuous(labels = function(x) sprintf("%.1f%%", 100 * x)) +
    labs(
      subtitle = "Diamond: posterior median | Thick: 50% CI | Thin: 95% CI | Red dot+bar: observed ± 95% CI",
      x = "Age group", y = y_lab
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.subtitle = element_text(size = 8.5, colour = "grey40"))
}

p_ci_prop <- plot_ci_bands(ci_prop, "Proportion of PWID") +
  labs(title = "Posterior predictive intervals: age proportions")
p_ci_prev <- plot_ci_bands(ci_prev, "HCV prevalence") +
  labs(title = "Posterior predictive intervals: HCV prevalence")

# =============================================================================
# 6.  ABSOLUTE POPULATION COUNT COMPARISON
# =============================================================================
# The proportion plots show relative age distribution, but absolute counts
# reveal whether the model's total population scale is also calibrated.
# Model counts are from the J-stratum at the posterior-mean parameters.
# They are scaled to sum(obs_tot) so the two series are directly comparable.
# A poor absolute fit (while proportions are correct) indicates the total
# J-stratum population size needs further calibration.

model_total_by_age <- combined$q_final$total_by_age
obs_survey_total   <- sum(obs_tot)
model_scaled       <- model_total_by_age / sum(model_total_by_age) * obs_survey_total

count_df <- data.frame(
  Age   = factor(rep(age_labels, 2), levels = age_labels),
  Count = c(model_scaled, obs_tot),
  Type  = rep(c("Model (scaled to survey N)", "Observed survey count"), each = 9)
)

p_abs_counts <- ggplot(count_df, aes(x = Age, y = Count, fill = Type)) +
  geom_col(position = position_dodge(width = 0.65), width = 0.6, alpha = 0.85) +
  scale_fill_manual(values = c("Model (scaled to survey N)" = "#4292C6",
                               "Observed survey count"      = "#D7191C")) +
  labs(
    title    = "Absolute PWID counts by age group",
    subtitle = sprintf(
      "Model J-stratum total = %d | Survey total = %d | Model scaled to survey N for comparison",
      round(sum(model_total_by_age)), obs_survey_total
    ),
    x = "Age group", y = "PWID count", fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position  = "top",
        plot.subtitle    = element_text(size = 8.5, colour = "grey40"))

# =============================================================================
# 7.  PRINT AND SAVE
# =============================================================================
print(p_trace_s1)
print(p_trace_s2)
print(p_dens_s1)
print(p_dens_s2)
print(p_dens_contact_scales)
print(p_ppc_prop)
print(p_ppc_prev)
print(p_ci_prop)
print(p_ci_prev)
print(p_abs_counts)

ggsave(file.path(OUT_DIR, "nc_trace_stage1.png"),         p_trace_s1,           width = 10, height = 8,   dpi = 200)
ggsave(file.path(OUT_DIR, "nc_trace_stage2.png"),         p_trace_s2,           width = 10, height = 8,   dpi = 200)
ggsave(file.path(OUT_DIR, "nc_density_contact_log.png"),  p_dens_s1,            width = 10, height = 8,   dpi = 200)
ggsave(file.path(OUT_DIR, "nc_density_contact_orig.png"), p_dens_contact_scales, width = 10, height = 8,   dpi = 200)
ggsave(file.path(OUT_DIR, "nc_density_beta_scale.png"),   p_dens_s2,            width = 10, height = 8,   dpi = 200)
ggsave(file.path(OUT_DIR, "nc_ppc_proportion_mean.png"),  p_ppc_prop,           width = 8,  height = 5,   dpi = 200)
ggsave(file.path(OUT_DIR, "nc_ppc_prevalence_mean.png"),  p_ppc_prev,           width = 8,  height = 5,   dpi = 200)
ggsave(file.path(OUT_DIR, "nc_ci_proportion.png"),        p_ci_prop,            width = 8,  height = 5,   dpi = 200)
ggsave(file.path(OUT_DIR, "nc_ci_prevalence.png"),        p_ci_prev,            width = 8,  height = 5,   dpi = 200)
ggsave(file.path(OUT_DIR, "nc_abs_counts.png"),           p_abs_counts,         width = 9,  height = 5,   dpi = 200)

saveRDS(list(ci_prop = ci_prop, ci_prev = ci_prev,
             ppc_prop_mat = ppc_prop_mat, ppc_prev_mat = ppc_prev_mat),
        file.path(OUT_DIR, "nc_ppc_intervals.rds"))

cat(sprintf("\nAll diagnostic plots saved to %s/\n", OUT_DIR))
