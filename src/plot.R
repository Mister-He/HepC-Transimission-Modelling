# Age-specific estimation plots
#
# Draws the age-specific one-year recidivism probability implied by the
# rearrest model and the age-specific release rate from the release model.
# Run from the project root with:
#   Rscript src/plot.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(posterior)
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
      legend.title = element_text(face = "bold")
    )
}

extract_indexed_draws <- function(draws, parameter, n_age = length(age_groups)) {
  parameter_names <- paste0(parameter, "[", seq_len(n_age), "]")
  draws_df <- as.data.frame(posterior::as_draws_df(draws))
  draws_df %>% select(all_of(parameter_names))
}

summarise_curve <- function(values, age_group, time_grid) {
  tibble(
    time = time_grid,
    age_group = factor(age_group, levels = age_groups),
    estimate = apply(values, 2, median),
    lower = apply(values, 2, quantile, probs = 0.025),
    upper = apply(values, 2, quantile, probs = 0.975)
  )
}

plot_survival <- function(data, title) {
  ggplot(data, aes(x = time, y = estimate, colour = age_group, fill = age_group)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA,
                show.legend = FALSE) +
    geom_line(linewidth = 0.9) +
    scale_colour_viridis_d(option = "C", direction = 1, end = 0.9) +
    scale_fill_viridis_d(option = "C", direction = 1, end = 0.9, guide = "none") +
    labs(
      x = ifelse(title == "Recidivism event", "Time since release (years)", "Time since rearrest (years)"),
      y = "Survival probability",
      title = title,
      colour = "Age Group"
    ) +
    theme_estimation()
}

dir.create("figures/estimation", recursive = TRUE, showWarnings = FALSE)

release_fit <- readRDS("output/release_rate_fit.rds")
rearrest_fit <- readRDS("output/rearrest_rate_fit.rds")

# The rearrest model defines the probability of a rearrest over an exposure
# interval with hazard h as pi * (1 - exp(-h)). Thus, the survival probability
# through time t in age group i is 1 - pi + pi * exp(-lambda[i] * t).
rearrest_draws <- rearrest_fit[[2]]
rearrest_parameters <- as.data.frame(posterior::as_draws_df(rearrest_draws))
pi_draws <- plogis(rearrest_parameters$pi_raw)
rearrest_lambda <- extract_indexed_draws(rearrest_draws, "log_lambda") %>%
  mutate(across(everything(), exp))

release_lambda <- extract_indexed_draws(release_fit$draws(variables = "lambda"), "lambda")

time_grid <- seq(0, 10, length.out = 201)

recidivism_survival <- bind_rows(lapply(seq_along(age_groups), function(i) {
  lambda <- rearrest_lambda[[i]]
  survival <- outer(lambda, time_grid, function(rate, time) {
    1 - pi_draws + pi_draws * exp(-rate * time)
  })
  summarise_curve(survival, age_groups[i], time_grid)
}))

release_survival <- bind_rows(lapply(seq_along(age_groups), function(i) {
  lambda <- release_lambda[[i]]
  survival <- outer(lambda, time_grid, function(rate, time) exp(-rate * time))
  summarise_curve(survival, age_groups[i], time_grid)
}))

recidivism_plot <- plot_survival(
  recidivism_survival,
  title = "Recidivism event"
)

release_plot <- plot_survival(
  release_survival,
  title = "Release event"
)

ggsave(
  "figures/estimation/estimation_survival_curves.png",
  (recidivism_plot + release_plot) +
    plot_layout(nrow = 1, guides = "collect") &
    theme(legend.position = "right", legend.direction = "vertical"),
  width = 12,
  height = 5.5,
  dpi = 300
)

message("Wrote figures/estimation/estimation_survival_curves.png")
