# =============================================================================
# plot_mcmc.R — ggplot2 figures for the Bayesian (MCMC) phase
# =============================================================================

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

okabe_ito_mcmc <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
                    "#E69F00", "#56B4E9", "#F0E442", "#000000")

theme_mcmc <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, colour = "grey60", linewidth = 0.4),
      axis.ticks = element_line(colour = "grey60"),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(colour = "grey30"),
      legend.position = "top",
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text = element_text(face = "bold", size = 9)
    )
}

param_labels <- function() {
  c(paste0("contact[", 1:6, "]"), paste0("beta[", 1:6, "]"))
}

plot_trace <- function(chains, out_file = NULL, thin = 2, width = 12, height = 10) {
  # chains: named list of matrices (n x 12), post-burnin
  d <- do.call(rbind, lapply(names(chains), function(nm) {
    m <- chains[[nm]]
    colnames(m) <- paste0("V", 1:12)
    idx <- seq(1, nrow(m), by = thin)
    data.frame(chain = nm, iter = idx,
               m[idx, , drop = FALSE])
  }))
  long <- d %>%
    tidyr::pivot_longer(c(-chain, -iter), names_to = "param", values_to = "value") %>%
    mutate(param = factor(param, levels = paste0("V", 1:12)))
  p <- ggplot(long, aes(x = iter, y = value, colour = chain)) +
    geom_line(linewidth = 0.25, alpha = 0.8) +
    facet_wrap(~param, scales = "free_y", ncol = 3,
               labeller = labeller(param = setNames(param_labels(), paste0("V", 1:12)))) +
    scale_colour_manual(values = okabe_ito_mcmc) +
    labs(x = "Iteration (post-burn-in)", y = "theta",
         title = "MCMC trace plots (all chains)",
         subtitle = "12 fitted log-scale parameters") +
    theme_mcmc()
  if (!is.null(out_file)) {
    ggsave(out_file, p, width = width, height = height, dpi = 300)
  }
  p
}

plot_density <- function(chains, priors, out_file = NULL,
                         width = 12, height = 10) {
  pooled <- do.call(rbind, chains)
  d <- as.data.frame(pooled)
  names(d) <- paste0("V", 1:12)
  long <- d %>%
    tidyr::pivot_longer(everything(), names_to = "param", values_to = "value")

  grid_pts <- seq(-6, 6, length.out = 400)
  prior_df <- do.call(rbind, lapply(1:6, function(i) {
    data.frame(param = paste0("V", i),
               x = grid_pts,
               y = dnorm(grid_pts, priors$contact_mean, priors$contact_sd))
  }))
  prior_df <- rbind(prior_df, do.call(rbind, lapply(7:12, function(i) {
    z <- (grid_pts - priors$beta_mean) / priors$beta_sd
    data.frame(param = paste0("V", i), x = grid_pts,
               y = dt(z, df = priors$beta_df) / priors$beta_sd)
  })))

  p <- ggplot(long, aes(x = value)) +
    geom_density(fill = okabe_ito_mcmc[1], alpha = 0.45, colour = okabe_ito_mcmc[1]) +
    geom_line(data = prior_df, aes(x = x, y = y, colour = "prior"),
              linetype = "dashed", linewidth = 0.8) +
    facet_wrap(~param, scales = "free", ncol = 3,
               labeller = labeller(param = setNames(param_labels(), paste0("V", 1:12)))) +
    scale_colour_manual(values = c(prior = okabe_ito_mcmc[2])) +
    labs(x = "theta", y = "density",
         title = "Marginal posterior densities vs priors",
         subtitle = "blue = posterior (3 chains pooled); red dashed = prior") +
    theme_mcmc()
  if (!is.null(out_file)) {
    ggsave(out_file, p, width = width, height = height, dpi = 300)
  }
  p
}

plot_predictive_density <- function(target_posterior, fit, out_file = NULL,
                                    width = 11, height = 8) {
  # target_posterior: data.frame with p1..p6, N1..N6 columns (equilibrium-feasible draws)
  age_levels <- fit$targets$age_groups
  obs_p <- fit$targets$prev_binom
  obs_p_lo <- qbeta(0.025, fit$targets$x_prev + 0.5,
                    fit$targets$n_prev - fit$targets$x_prev + 0.5)
  obs_p_hi <- qbeta(0.975, fit$targets$x_prev + 0.5,
                    fit$targets$n_prev - fit$targets$x_prev + 0.5)
  obs_N_lo <- exp(log(fit$targets$prison_total) - qnorm(0.975) * 0.10)
  obs_N_hi <- exp(log(fit$targets$prison_total) + qnorm(0.975) * 0.10)

  p_long <- target_posterior %>%
    select(p1:p6) %>%
    tidyr::pivot_longer(everything(), names_to = "age", values_to = "p") %>%
    mutate(age = factor(sub("p", "", age), levels = as.character(1:6)))
  obs_p_df <- data.frame(
    age = factor(1:6, levels = 1:6),
    lo = obs_p_lo, hi = obs_p_hi, mid = obs_p
  )
  p1 <- ggplot(p_long, aes(x = p)) +
    geom_histogram(bins = 50, fill = okabe_ito_mcmc[1], alpha = 0.55) +
    geom_errorbarh(data = obs_p_df, aes(xmin = lo, xmax = hi, y = Inf),
                   inherit.aes = FALSE, colour = okabe_ito_mcmc[2], linewidth = 1.0, height = 0) +
    facet_wrap(~age, scales = "free", ncol = 3,
               labeller = labeller(age = setNames(age_levels, 1:6))) +
    labs(x = "Posterior predictive seroprevalence", y = "count",
         title = "Posterior predictive prevalence by age group",
         subtitle = "red bar = observed 95% interval") +
    theme_mcmc()

  N_long <- target_posterior %>%
    select(N1:N6) %>%
    tidyr::pivot_longer(everything(), names_to = "age", values_to = "N") %>%
    mutate(age = factor(sub("N", "", age), levels = as.character(1:6)))
  obs_N_df <- data.frame(
    age = factor(1:6, levels = 1:6),
    lo = obs_N_lo, hi = obs_N_hi, mid = fit$targets$prison_total
  )
  p2 <- ggplot(N_long, aes(x = N)) +
    geom_histogram(bins = 50, fill = okabe_ito_mcmc[3], alpha = 0.55) +
    geom_errorbarh(data = obs_N_df, aes(xmin = lo, xmax = hi, y = Inf),
                   inherit.aes = FALSE, colour = okabe_ito_mcmc[2], linewidth = 1.0, height = 0) +
    facet_wrap(~age, scales = "free", ncol = 3,
               labeller = labeller(age = setNames(age_levels, 1:6))) +
    labs(x = "Posterior predictive prison population", y = "count",
         title = "Posterior predictive population by age group",
         subtitle = "red bar = observed 95% interval") +
    theme_mcmc()

  p <- p1 / p2
  if (!is.null(out_file)) {
    ggsave(out_file, p, width = width, height = height, dpi = 300)
  }
  p
}

plot_ci_compare <- function(ci, fit, out_file = NULL, width = 11, height = 8) {
  # ci: credible_intervals.csv-style data.frame with columns
  # age_group, p_median, p_lo, p_hi, p_la_lo, p_la_hi, p_obs_lo, p_obs_hi,
  # N_median, N_lo, N_hi, N_la_lo, N_la_hi, N_obs_lo, N_obs_hi
  age_levels <- fit$targets$age_groups
  d_p <- data.frame(
    age = factor(rep(ci$age_group, 3), levels = age_levels),
    source = factor(rep(c("MCMC", "Laplace", "Observed"), each = nrow(ci)),
                    levels = c("MCMC", "Laplace", "Observed")),
    mid = c(ci$p_median, fit$laplace$intervals$p_hat,
            fit$targets$prev_binom),
    lo = c(ci$p_lo, fit$laplace$intervals$p_lo, ci$p_obs_lo),
    hi = c(ci$p_hi, fit$laplace$intervals$p_hi, ci$p_obs_hi)
  )
  p1 <- ggplot(d_p, aes(x = age, colour = source)) +
    geom_point(aes(y = mid), position = position_dodge(width = 0.55), size = 2.2) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12,
                  position = position_dodge(width = 0.55), linewidth = 0.9) +
    scale_colour_manual(values = c(MCMC = okabe_ito_mcmc[1],
                                   Laplace = okabe_ito_mcmc[4],
                                   Observed = okabe_ito_mcmc[2])) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Age group", y = "Seroprevalence",
         title = "Seroprevalence: MCMC CrI vs Laplace CI vs observed") +
    theme_mcmc()

  d_N <- data.frame(
    age = factor(rep(ci$age_group, 3), levels = age_levels),
    source = factor(rep(c("MCMC", "Laplace", "Observed"), each = nrow(ci)),
                    levels = c("MCMC", "Laplace", "Observed")),
    mid = c(ci$N_median, fit$laplace$intervals$N_hat,
            fit$targets$prison_total),
    lo = c(ci$N_lo, fit$laplace$intervals$N_lo, ci$N_obs_lo),
    hi = c(ci$N_hi, fit$laplace$intervals$N_hi, ci$N_obs_hi)
  )
  p2 <- ggplot(d_N, aes(x = age, colour = source)) +
    geom_point(aes(y = mid), position = position_dodge(width = 0.55), size = 2.2) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12,
                  position = position_dodge(width = 0.55), linewidth = 0.9) +
    scale_colour_manual(values = c(MCMC = okabe_ito_mcmc[1],
                                   Laplace = okabe_ito_mcmc[4],
                                   Observed = okabe_ito_mcmc[2])) +
    labs(x = "Age group", y = "Prison population",
         title = "Population: MCMC CrI vs Laplace CI vs observed") +
    theme_mcmc()

  p <- p1 / p2
  if (!is.null(out_file)) {
    ggsave(out_file, p, width = width, height = height, dpi = 300)
  }
  p
}

# NPE version: posterior densities of NPE seed1, NPE seed2, MCMC, prior
plot_density_npe <- function(npe, npe2, mcmc, priors, out_file = NULL,
                             width = 13, height = 11) {
  to_long <- function(df, src) {
    m <- as.matrix(df[, paste0("theta", 1:12)])
    colnames(m) <- paste0("V", 1:12)
    as.data.frame(m) %>%
      mutate(source = src) %>%
      tidyr::pivot_longer(c(-source), names_to = "param", values_to = "value")
  }
  d <- rbind(to_long(npe, "NPE seed 1"),
             to_long(npe2, "NPE seed 2"),
             to_long(mcmc, "MCMC validation"))
  d$param <- factor(d$param, levels = paste0("V", 1:12))

  grid_pts <- seq(-6, 6, length.out = 400)
  prior_df <- do.call(rbind, lapply(1:12, function(i) {
    if (i <= 6) {
      y <- dnorm(grid_pts, priors$contact_mean, priors$contact_sd)
    } else {
      z <- (grid_pts - priors$beta_mean) / priors$beta_sd
      y <- dt(z, df = priors$beta_df) / priors$beta_sd
    }
    data.frame(param = paste0("V", i), x = grid_pts, y = y)
  }))

  p <- ggplot(d, aes(x = value, colour = source)) +
    geom_density(linewidth = 0.7) +
    geom_line(data = prior_df, aes(x = x, y = y, colour = "prior"),
              linetype = "dashed", linewidth = 0.7) +
    facet_wrap(~param, scales = "free", ncol = 3,
               labeller = labeller(param = setNames(param_labels(), paste0("V", 1:12)))) +
    scale_colour_manual(values = c("NPE seed 1" = okabe_ito_mcmc[1],
                                   "NPE seed 2" = okabe_ito_mcmc[5],
                                   "MCMC validation" = okabe_ito_mcmc[2],
                                   "prior" = "grey45")) +
    labs(x = "theta", y = "density",
         title = "Posterior densities: NPE (2 seeds) vs MCMC validation vs prior") +
    theme_mcmc()
  if (!is.null(out_file)) ggsave(out_file, p, width = width, height = height, dpi = 300)
  p
}

# MCMC-only marginal densities: 2 rows x 6 columns.
# Row 1 = contact scaling factors (theta1..theta6 -> exp), row 2 = beta
# scaling factors (theta7..theta12 -> exp), x on log10 scale, dashed line
# at scale = 1. Three chains overlaid.
plot_density_mcmc_chains <- function(mcmc, out_file = NULL,
                                     width = 13, height = 7) {
  sc <- exp(as.matrix(mcmc[, paste0("theta", 1:12)]))
  colnames(sc) <- paste0("V", 1:12)
  d <- as.data.frame(sc) %>%
    mutate(chain = mcmc$chain) %>%
    tidyr::pivot_longer(c(-chain), names_to = "param", values_to = "scale")
  d$param <- factor(d$param, levels = paste0("V", 1:12))

  lab <- c(paste0("contact scale ", 1:6), paste0("beta scale ", 1:6))
  p <- ggplot(d, aes(x = scale, colour = chain)) +
    geom_density(linewidth = 0.7) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    scale_x_log10() +
    facet_wrap(~param, ncol = 6, scales = "free",
               labeller = labeller(param = setNames(lab, paste0("V", 1:12)))) +
    scale_colour_manual(values = okabe_ito_mcmc) +
    labs(x = "fitted scale factor (log10 axis)", y = "density",
         title = "MCMC marginal posterior densities (3 chains)",
         subtitle = paste0("row 1: contact scaling factors | ",
                           "row 2: beta scaling factors | dashed line = 1")) +
    theme_mcmc()
  if (!is.null(out_file)) {
    ggsave(out_file, p, width = width, height = height, dpi = 300)
  }
  p
}

# NPE version: four interval sources (NPE, MCMC, Laplace, observed)
plot_ci_compare_npe <- function(ci, fit, out_file = NULL, width = 11, height = 8) {
  age_levels <- fit$targets$age_groups
  d_p <- data.frame(
    age = factor(rep(ci$age_group, 4), levels = age_levels),
    source = factor(rep(c("NPE", "MCMC", "Laplace", "Observed"),
                        each = nrow(ci)),
                    levels = c("NPE", "MCMC", "Laplace", "Observed")),
    mid = c(ci$p_npe_median, ci$p_mcmc_median, ci$p_la_hat,
            fit$targets$prev_binom),
    lo = c(ci$p_npe_lo, ci$p_mcmc_lo, ci$p_la_lo, ci$p_obs_lo),
    hi = c(ci$p_npe_hi, ci$p_mcmc_hi, ci$p_la_hi, ci$p_obs_hi)
  )
  p1 <- ggplot(d_p, aes(x = age, colour = source)) +
    geom_point(aes(y = mid), position = position_dodge(width = 0.6), size = 2.0) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12,
                  position = position_dodge(width = 0.6), linewidth = 0.85) +
    scale_colour_manual(values = c(NPE = okabe_ito_mcmc[1],
                                   MCMC = okabe_ito_mcmc[3],
                                   Laplace = okabe_ito_mcmc[4],
                                   Observed = okabe_ito_mcmc[2])) +
    scale_y_continuous(labels = scales::percent) +
    labs(x = "Age group", y = "Seroprevalence",
         title = "Seroprevalence: NPE vs MCMC vs Laplace vs observed") +
    theme_mcmc()

  d_N <- data.frame(
    age = factor(rep(ci$age_group, 4), levels = age_levels),
    source = factor(rep(c("NPE", "MCMC", "Laplace", "Observed"),
                        each = nrow(ci)),
                    levels = c("NPE", "MCMC", "Laplace", "Observed")),
    mid = c(ci$N_npe_median, ci$N_mcmc_median, ci$N_la_hat,
            fit$targets$prison_total),
    lo = c(ci$N_npe_lo, ci$N_mcmc_lo, ci$N_la_lo, ci$N_obs_lo),
    hi = c(ci$N_npe_hi, ci$N_mcmc_hi, ci$N_la_hi, ci$N_obs_hi)
  )
  p2 <- ggplot(d_N, aes(x = age, colour = source)) +
    geom_point(aes(y = mid), position = position_dodge(width = 0.6), size = 2.0) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12,
                  position = position_dodge(width = 0.6), linewidth = 0.85) +
    scale_colour_manual(values = c(NPE = okabe_ito_mcmc[1],
                                   MCMC = okabe_ito_mcmc[3],
                                   Laplace = okabe_ito_mcmc[4],
                                   Observed = okabe_ito_mcmc[2])) +
    labs(x = "Age group", y = "Prison population",
         title = "Population: NPE vs MCMC vs Laplace vs observed") +
    theme_mcmc()

  p <- p1 / p2
  if (!is.null(out_file)) ggsave(out_file, p, width = width, height = height, dpi = 300)
  p
}

plot_sbc <- function(sbc, out_file = NULL, width = 10, height = 5) {
  d <- sbc %>%
    mutate(param = factor(parameter, levels = paste0("theta", 1:12))) %>%
    tidyr::pivot_longer(c(coverage_95, coverage_90),
                        names_to = "level", values_to = "coverage") %>%
    mutate(target = ifelse(level == "coverage_95", 0.95, 0.90),
           level = recode(level, coverage_95 = "95% band", coverage_90 = "90% band"))

  p <- ggplot(d, aes(x = param, y = coverage, fill = level)) +
    geom_col(position = "dodge", width = 0.7) +
    geom_hline(aes(yintercept = target, colour = level), linetype = "dashed") +
    scale_fill_manual(values = c("95% band" = okabe_ito_mcmc[1],
                                 "90% band" = okabe_ito_mcmc[5])) +
    scale_colour_manual(values = c("95% band" = okabe_ito_mcmc[1],
                                   "90% band" = okabe_ito_mcmc[5])) +
    coord_flip() +
    labs(x = NULL, y = "Empirical SBC coverage",
         title = "Simulation-based calibration: posterior coverage per parameter",
         subtitle = "dashed lines = nominal 95% / 90% coverage") +
    theme_mcmc() + theme(panel.grid.major.y = element_blank())
  if (!is.null(out_file)) ggsave(out_file, p, width = width, height = height, dpi = 300)
  p
}
