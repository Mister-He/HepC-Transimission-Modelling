# Compare observed age-specific HCV prevalence and prison population with
# posterior-predictive MCMC estimates.
#
# Run from the project root with:
#   Rscript src/fig_est_hcv_tot.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

age_groups <- c("<20", "20-29", "30-39", "40-49", "50-59", "60+")

theme_estimation <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, colour = "grey60", linewidth = 0.4),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      axis.ticks = element_line(colour = "grey60"),
      plot.title = element_text(face = "bold", hjust = 0),
      legend.position = "right",
      legend.direction = "vertical",
      legend.title.position = "top",
      legend.box = "vertical"
    )
}

root <- getwd()
fit_path <- file.path(root, "output", "calibration", "run1_4strata", "fit.rds")
ci_path <- file.path(root, "output", "calibration", "npe_bayes", "credible_intervals.csv")
mcmc_path <- file.path(root, "output", "calibration", "npe_bayes",
                       "posterior_samples_mcmc.csv")
out_path <- file.path(root, "figures", "estimation", "fig_est_hcv_tot.png")

required_paths <- c(fit_path, ci_path, mcmc_path)
if (!all(file.exists(required_paths))) {
  stop("Required calibration output is missing: ",
       paste(required_paths[!file.exists(required_paths)], collapse = ", "))
}

fit <- readRDS(fit_path)
ci <- read.csv(ci_path, check.names = FALSE)
mcmc_draws <- read.csv(mcmc_path, check.names = FALSE)

if (!identical(ci$age_group, age_groups)) {
  ci <- ci[match(age_groups, ci$age_group), , drop = FALSE]
}
if (anyNA(ci$age_group)) stop("credible_intervals.csv does not contain all age groups")

required_ci <- c(
  "p_mcmc_median", "p_mcmc_lo", "p_mcmc_hi",
  "N_mcmc_median", "N_mcmc_lo", "N_mcmc_hi"
)
if (!all(required_ci %in% names(ci))) {
  stop("credible_intervals.csv is missing required MCMC columns")
}

# Observed prevalence uses the reconstructed Binomial proportion. Its 95% CI
# is the Jeffreys interval used throughout the calibration reports. Population
# uncertainty uses the reported log-Normal observation model (sigma = 0.10).
observed_prevalence <- fit$targets$prev_binom
observed_population <- fit$targets$prison_total
observed_prevalence_lo <- qbeta(
  0.025,
  fit$targets$x_prev + 0.5,
  fit$targets$n_prev - fit$targets$x_prev + 0.5
)
observed_prevalence_hi <- qbeta(
  0.975,
  fit$targets$x_prev + 0.5,
  fit$targets$n_prev - fit$targets$x_prev + 0.5
)
observed_population_lo <- exp(
  log(observed_population) - qnorm(0.975) * 0.10
)
observed_population_hi <- exp(
  log(observed_population) + qnorm(0.975) * 0.10
)

comparison_data <- function(observed, observed_lo, observed_hi,
                            model, model_lo, model_hi) {
  bind_rows(
    tibble(
      age_group = age_groups,
      source = "Observed",
      estimate = observed,
      lower = observed_lo,
      upper = observed_hi
    ),
    tibble(
      age_group = age_groups,
      source = "MCMC posterior",
      estimate = model,
      lower = model_lo,
      upper = model_hi
    )
  ) %>%
    mutate(
      age_group = factor(age_group, levels = age_groups),
      source = factor(source, levels = c("Observed", "MCMC posterior"))
    )
}

plot_comparison <- function(data, y_label, title, y_labels = waiver()) {
  dodge <- position_dodge(width = 0.55)

  ggplot(data, aes(x = age_group, colour = age_group, shape = source)) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper),
      position = dodge,
      width = 0.13,
      linewidth = 0.8
    ) +
    geom_point(
      aes(y = estimate),
      position = dodge,
      size = 2.8
    ) +
    scale_colour_viridis_d(option = "C", direction = 1, end = 0.9,
                           name = "Age Group") +
    scale_shape_manual(
      values = c("Observed" = 17, "MCMC posterior" = 16),
      name = "Source"
    ) +
    scale_y_continuous(labels = y_labels,
                       limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(
      x = "Age group",
      y = y_label,
      title = title
    ) +
    theme_estimation()
}

prevalence_data <- comparison_data(
  observed = observed_prevalence,
  observed_lo = observed_prevalence_lo,
  observed_hi = observed_prevalence_hi,
  model = ci$p_mcmc_median,
  model_lo = ci$p_mcmc_lo,
  model_hi = ci$p_mcmc_hi
)

population_data <- comparison_data(
  observed = observed_population,
  observed_lo = observed_population_lo,
  observed_hi = observed_population_hi,
  model = ci$N_mcmc_median,
  model_lo = ci$N_mcmc_lo,
  model_hi = ci$N_mcmc_hi
)

prevalence_plot <- plot_comparison(
  prevalence_data,
  y_label = "HCV seroprevalence",
  title = "Age-specific HCV seroprevalence",
  y_labels = scales::label_percent(accuracy = 1)
)

population_plot <- plot_comparison(
  population_data,
  y_label = "Prison population",
  title = "Age-specific prison population",
  y_labels = scales::label_comma(accuracy = 1)
)

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
ggsave(
  out_path,
  (prevalence_plot / population_plot) +
    plot_layout(guides = "collect") &
    theme(legend.position = "right", legend.direction = "vertical"),
  width = 10,
  height = 11,
  dpi = 300
)

message("Used ", nrow(mcmc_draws), " retained MCMC posterior draws.")
message("Wrote ", out_path)
