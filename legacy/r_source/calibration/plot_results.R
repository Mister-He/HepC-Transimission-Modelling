# =============================================================================
# plot_results.R — ALL figures via ggplot2 (publication quality)
#
# Functions only; call plot_all(fit, out_dir) to render every figure.
# CLI usage:
#   Rscript scripts/plot_results.R --fit=<fit.rds> --out-dir=<dir>
# =============================================================================

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

okabe_ito <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
               "#E69F00", "#56B4E9", "#F0E442", "#000000")

theme_pub <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, colour = "grey60", linewidth = 0.4),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      axis.ticks = element_line(colour = "grey60"),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(colour = "grey30"),
      legend.position = "top",
      legend.title = element_blank(),
      strip.background = element_rect(fill = "white", colour = NA),
      strip.text = element_text(face = "bold")
    )
}

fig_prevalence_fit <- function(fit) {
  pred <- fit$predictions
  int  <- fit$laplace$intervals
  age_levels <- fit$targets$age_groups
  d <- pred %>%
    mutate(age = factor(age_group, levels = age_levels))
  int$age <- factor(int$age_group, levels = age_levels)

  ggplot(d, aes(x = age)) +
    geom_errorbar(data = int, aes(ymin = p_obs_lo, ymax = p_obs_hi),
                  width = 0.15, linewidth = 0.8, colour = okabe_ito[2]) +
    geom_errorbar(data = int, aes(ymin = p_lo, ymax = p_hi),
                  width = 0.08, linewidth = 1.0, colour = okabe_ito[1]) +
    geom_point(aes(y = p_obs_binom, colour = "Observed (Binomial 95% CI)"),
               shape = 17, size = 2.6) +
    geom_point(aes(y = p_hat_sero, colour = "Model seroprevalence (Laplace 95% CI)"),
               shape = 19, size = 2.8) +
    scale_colour_manual(values = c("Observed (Binomial 95% CI)" = okabe_ito[2],
                                   "Model seroprevalence (Laplace 95% CI)" = okabe_ito[1])) +
    scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
    labs(x = "Prison age group", y = "HCV seroprevalence",
         title = "Prison HCV seroprevalence: model vs target (2015)",
         subtitle = paste0("Run: ", fit$run_id)) +
    theme_pub()
}

fig_population_fit <- function(fit) {
  pred <- fit$predictions
  int  <- fit$laplace$intervals
  age_levels <- fit$targets$age_groups
  d <- pred %>% mutate(age = factor(age_group, levels = age_levels))
  int$age <- factor(int$age_group, levels = age_levels)

  ggplot(d, aes(x = age)) +
    geom_errorbar(data = int, aes(ymin = N_obs_lo, ymax = N_obs_hi),
                  width = 0.15, linewidth = 0.8, colour = okabe_ito[2]) +
    geom_errorbar(data = int, aes(ymin = N_lo, ymax = N_hi),
                  width = 0.08, linewidth = 1.0, colour = okabe_ito[1]) +
    geom_point(aes(y = N_obs, colour = "Observed (log-Normal 95% CI)"),
               shape = 17, size = 2.6) +
    geom_point(aes(y = N_hat, colour = "Model (Laplace 95% CI)"),
               shape = 19, size = 2.8) +
    scale_colour_manual(values = c("Observed (log-Normal 95% CI)" = okabe_ito[2],
                                   "Model (Laplace 95% CI)" = okabe_ito[1])) +
    labs(x = "Prison age group", y = "Prison population",
         title = "Prison population: model vs target (2015)",
         subtitle = paste0("Run: ", fit$run_id)) +
    theme_pub()
}

fig_residuals <- function(fit) {
  res <- fit$residuals
  age_levels <- fit$targets$age_groups
  d <- res %>%
    mutate(age = factor(age_group, levels = age_levels)) %>%
    tidyr::pivot_longer(c(prev_residual, log_pop_residual),
                        names_to = "type", values_to = "residual") %>%
    mutate(type = recode(type,
                         prev_residual = "Prevalence residual (p_hat - p_obs)",
                         log_pop_residual = "log-population residual"))

  ggplot(d, aes(x = age, y = residual, colour = type)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(size = 2.6) +
    geom_line(aes(group = type), linewidth = 0.8) +
    scale_colour_manual(values = c("Prevalence residual (p_hat - p_obs)" = okabe_ito[1],
                                   "log-population residual" = okabe_ito[2])) +
    labs(x = "Prison age group", y = "Residual",
         title = "Fit residuals by age group",
         subtitle = paste0("Run: ", fit$run_id)) +
    theme_pub()
}

fig_parameter_scales <- function(fit) {
  sol <- fit$solutions[which.min(fit$solutions$objective), ]
  age_levels <- fit$targets$age_groups

  d_contact <- data.frame(
    age = factor(age_levels, levels = age_levels),
    scale = as.numeric(sol[paste0("contact_scale", 1:6)]),
    type = "Contact row scale"
  )
  d_beta <- data.frame(
    age = factor(age_levels, levels = age_levels),
    scale = as.numeric(sol[paste0("beta_scale", 1:6)]),
    type = "beta inflow scale"
  )
  d <- rbind(d_contact, d_beta)

  ggplot(d, aes(x = age, y = scale, fill = type)) +
    geom_col(position = "dodge", width = 0.7, alpha = 0.9) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    scale_fill_manual(values = c("Contact row scale" = okabe_ito[1],
                                 "beta inflow scale" = okabe_ito[2])) +
    scale_y_log10() +
    labs(x = "Age group", y = "Fitted scale (log axis)",
         title = "Fitted contact and beta scaling factors",
         subtitle = paste0("Run: ", fit$run_id,
                           "; dashed line = 1 (beta_scale > 1 is the soft target)")) +
    theme_pub()
}

fig_objective_starts <- function(fit) {
  sol <- fit$solutions
  d <- sol %>%
    mutate(start_id = forcats::fct_reorder(start_id, objective))
  ggplot(d, aes(x = start_id, y = objective, fill = objective == min(objective))) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c("FALSE" = "grey75", "TRUE" = okabe_ito[1]),
                      guide = "none") +
    coord_flip() +
    labs(x = NULL, y = "Final objective (NLL)",
         title = "Multi-start Nelder-Mead: final objective per start",
         subtitle = paste0("Run: ", fit$run_id,
                           "; highlighted = best start")) +
    theme_pub() + theme(panel.grid.major.y = element_blank())
}

fig_equilibrium_stability <- function(fit) {
  eq <- fit$equilibrium$by_age
  age_levels <- fit$targets$age_groups
  d <- eq %>%
    mutate(age = factor(age_group, levels = age_levels))

  p1 <- ggplot(d, aes(x = age)) +
    geom_point(aes(y = N_hat_T, colour = "T"), size = 2.6) +
    geom_point(aes(y = N_hat_T5, colour = "T-5"), shape = 17, size = 2.6) +
    scale_colour_manual(values = c("T" = okabe_ito[1], "T-5" = okabe_ito[2])) +
    labs(x = "Age group", y = "J population", title = "J population: T vs T-5") +
    theme_pub()

  p2 <- ggplot(d, aes(x = age)) +
    geom_point(aes(y = p_hat_T, colour = "T"), size = 2.6) +
    geom_point(aes(y = p_hat_T5, colour = "T-5"), shape = 17, size = 2.6) +
    scale_colour_manual(values = c("T" = okabe_ito[1], "T-5" = okabe_ito[2])) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Age group", y = "J seroprevalence", title = "J seroprevalence: T vs T-5") +
    theme_pub()

  p1 + p2 + plot_annotation(
    title = "Equilibrium stability (T vs T-5)",
    subtitle = paste0("Run: ", fit$run_id, "; max |prev change| = ",
                      sprintf("%.4f", fit$equilibrium$max_prev_change),
                      " (criterion 0.005)"))
}

fig_dashboard <- function(fit) {
  (fig_prevalence_fit(fit) + fig_population_fit(fit)) /
    (fig_residuals(fit) + fig_parameter_scales(fit)) +
    plot_annotation(
      title = paste0("Calibration dashboard — ", fit$run_id),
      subtitle = sprintf("NLL %.2f | prev RMSE %.4f | pop MAPE %.3f",
                         fit$metrics$nll_prev + fit$metrics$nll_pop,
                         fit$metrics$prevalence_rmse,
                         fit$metrics$population_mape),
      theme = theme(plot.title = element_text(face = "bold", size = 16))
    )
}

plot_all <- function(fit, out_dir, width = 9, height = 6, dpi = 300) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(out_dir, "fig01_prevalence_fit.png"),
         fig_prevalence_fit(fit), width = width, height = height, dpi = dpi)
  ggsave(file.path(out_dir, "fig02_population_fit.png"),
         fig_population_fit(fit), width = width, height = height, dpi = dpi)
  ggsave(file.path(out_dir, "fig03_residuals.png"),
         fig_residuals(fit), width = width, height = height, dpi = dpi)
  ggsave(file.path(out_dir, "fig04_parameter_scales.png"),
         fig_parameter_scales(fit), width = width, height = height, dpi = dpi)
  ggsave(file.path(out_dir, "fig05_objective_starts.png"),
         fig_objective_starts(fit), width = width, height = height, dpi = dpi)
  ggsave(file.path(out_dir, "fig06_equilibrium_stability.png"),
         fig_equilibrium_stability(fit), width = width + 2, height = height, dpi = dpi)
  ggsave(file.path(out_dir, "fig07_fit_dashboard.png"),
         fig_dashboard(fit), width = 12, height = 10, dpi = dpi)
  invisible(file.path(out_dir, paste0("fig0", 1:7, ".png")))
}

# CLI guard: only runs when invoked directly with --fit.
if (length(commandArgs(trailingOnly = TRUE)) > 0 &&
    any(grepl("^--fit=", commandArgs(trailingOnly = TRUE)))) {
  cli_args <- commandArgs(trailingOnly = TRUE)
  arg_val <- function(name, default = NULL) {
    hit <- cli_args[grepl(paste0("^", name, "="), cli_args)]
    if (length(hit) == 0) return(default)
    sub(paste0("^", name, "="), "", hit[1])
  }
  fit <- readRDS(arg_val("--fit"))
  out_dir <- arg_val("--out-dir", "figures")
  plot_all(fit, out_dir)
  cat("Figures written to", out_dir, "\n")
}
