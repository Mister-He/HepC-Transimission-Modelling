# =============================================================================
# Nelder-Mead MAP calibration for the 20-parameter HCV PWID model
#
# Run from this directory with:
#   Rscript nelder_mead.R
#
# Or source the file and call fit_nelder_mead() with custom settings.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
})

# Observations and parameter metadata required by HMC_core.r. Keep these values
# identical to hmc.R so both stages use the same binomial prevalence likelihood.
obs_prev <- c(
  0.1118421, 0.1470588, 0.1933842, 0.2507599, 0.2899083,
  0.3596059, 0.5025295, 0.5061728, 0.4534314, 0.3544304
)
obs_prev_se <- sqrt(c(
  0.09940967, 0.12867003, 0.16065544, 0.18991135, 0.21719265,
  0.23775531, 0.24736428, 0.24839586, 0.24995891, 0.24000000
))
obs_tot <- c(99, 552, 692, 763, 704, 847, 994, 847, 781, 409)
obs_pos <- round(obs_prev * obs_tot)  # HCV positives per age group (binomial numerator)
N_AGE <- length(obs_tot)

# Interpretable (age-level) parameter names, used for theta_to_orig() output.
param_names_orig <- c(
  paste0("C_contact_scale_", seq_len(N_AGE)),
  paste0("tot_in_scaling_fct_", seq_len(N_AGE))
)

source("setup.R")
source("HMC_core.r")

# Fitted parameters are the spline coefficients (2*SPLINE_K), not age values.
param_names_log <- c(
  paste0("alpha_contact_", seq_len(SPLINE_K)),
  paste0("gamma_tot_in_", seq_len(SPLINE_K))
)

# Shared observation-model settings used in hmc.R. Prevalence is binomial;
# prev_logit_sd only affects plot intervals, not the likelihood.
data$prev_logit_sd <- 0.10
data$sigma_pop <- c(0.20, rep(0.12, N_AGE - 1L))
data$sigma_shape <- 0.20
data$nu_shape <- 7L

# Return all model quantities used in fitting at one theta value.
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

# Fit one or more Nelder-Mead starts and retain every objective evaluation.
fit_nelder_mead <- function(
    theta_init = c(CONTACT_COEF_PRIOR_MEANS, TOT_IN_COEF_PRIOR_MEANS),
    base_params = params,
    sim_data = data,
    n_starts = 6L,
    start_sd = 0.25,
    maxit = 1e+05,
    reltol = 1e-6,
    seed = 42L,
    report_every = 500L) {

  stopifnot(length(theta_init) == 2L * SPLINE_K, n_starts >= 1L)
  set.seed(seed)
  starts <- vector("list", n_starts)
  starts[[1L]] <- theta_init
  if (n_starts > 1L) {
    for (i in 2:n_starts) starts[[i]] <- theta_init + rnorm(length(theta_init), 0, start_sd)
  }
  if (n_starts >= 3L) {
    # Deterministic early-age prevalence starts. These keep the 5-knot spline
    # constraints but help avoid the local optimum with a near-zero age-1 curve.
    starts[[2L]][1L] <- starts[[2L]][1L] + 1.0
    starts[[3L]][1L] <- starts[[3L]][1L] + 1.8
  }

  fits <- vector("list", n_starts)
  for (s in seq_len(n_starts)) {
    history <- list()
    n_eval <- 0L
    objective <- function(theta) {
      n_eval <<- n_eval + 1L
      lp <- log_posterior(theta, base_params, sim_data)
      value <- if (is.finite(lp)) -lp else .Machine$double.xmax / 100
      history[[n_eval]] <<- c(eval = n_eval, neg_log_posterior = value, theta)
      if (report_every > 0L && n_eval %% report_every == 0L) {
        cat(sprintf("start %d: evaluation %d, best negative log-posterior = %.6f\n",
                    s, n_eval,
                    min(vapply(history, function(z) z[[2L]], numeric(1L)))))
      }
      value
    }

    cat(sprintf("Starting Nelder-Mead run %d/%d\n", s, n_starts))
    opt <- optim(
      par = starts[[s]], fn = objective, method = "Nelder-Mead",
      control = list(maxit = maxit, reltol = reltol)
    )
    hist_mat <- do.call(rbind, history)
    colnames(hist_mat) <- c("evaluation", "neg_log_posterior", param_names_log)
    fits[[s]] <- list(optim = opt, history = as.data.frame(hist_mat), start = starts[[s]])
  }

  values <- vapply(fits, function(x) x$optim$value, numeric(1L))
  best_id <- which.min(values)
  best <- fits[[best_id]]
  theta_hat <- setNames(best$optim$par, param_names_log)
  pred <- predict_at_theta(theta_hat, base_params, sim_data)

  result <- list(
    theta_hat = theta_hat,
    # Reconstruct interpretable age-level scales from the spline coefficients.
    par_hat = setNames(as.numeric(theta_to_orig(theta_hat)), param_names_orig),
    log_posterior = -best$optim$value,
    log_likelihood = log_likelihood(theta_hat, base_params, sim_data),
    log_prior = log_prior(theta_hat),
    convergence = best$optim$convergence,
    message = best$optim$message,
    counts = best$optim$counts,
    best_start = best_id,
    prediction = pred,
    spline_prior = list(
      contact_mean = CONTACT_COEF_PRIOR_MEANS,
      tot_in_mean = TOT_IN_COEF_PRIOR_MEANS,
      contact_sd = CONTACT_COEF_PRIOR_SDS,
      tot_in_sd = TOT_IN_COEF_PRIOR_SDS,
      rw_sd = SPLINE_RW_SD
    ),
    fits = fits,
    base_params = base_params,
    sim_data = sim_data
  )
  class(result) <- "nm_hcv_fit"
  result
}

print.nm_hcv_fit <- function(x, ...) {
  cat("Nelder-Mead MAP fit\n")
  cat(sprintf("  best start: %d\n", x$best_start))
  cat(sprintf("  convergence code: %d (0 means successful)\n", x$convergence))
  cat(sprintf("  log posterior: %.6f\n", x$log_posterior))
  cat(sprintf("  log likelihood: %.6f\n", x$log_likelihood))
  cat(sprintf("  log prior: %.6f\n", x$log_prior))
  print(round(x$par_hat, 5))
  invisible(x)
}

plot_nm_convergence <- function(fit) {
  dat <- do.call(rbind, lapply(seq_along(fit$fits), function(i) {
    x <- fit$fits[[i]]$history
    x$best_so_far <- cummin(x$neg_log_posterior)
    x$start <- factor(i)
    x
  }))
  ggplot(dat, aes(evaluation, best_so_far, colour = start)) +
    geom_line(linewidth = 0.7) +
    labs(x = "Objective evaluations", y = "Best negative log-posterior",
         colour = "Start", title = "Nelder-Mead convergence") +
    theme_bw()
}

plot_nm_parameters <- function(fit) {
  prior_mean <- exp(c(CONTACT_PRIOR_MEANS, TOT_IN_PRIOR_MEANS))
  dat <- data.frame(
    parameter = factor(param_names_orig, levels = param_names_orig),
    prior_mean = prior_mean,
    estimate = as.numeric(fit$par_hat),
    family = rep(c("Contact scale", "Total-in scale"), each = N_AGE)
  )
  ggplot(dat, aes(prior_mean, estimate, colour = family, label = seq_len(nrow(dat)))) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey45") +
    geom_point(size = 2.5) +
    geom_text(nudge_y = 0.03, size = 3, show.legend = FALSE) +
    scale_x_log10() + scale_y_log10() +
    labs(x = "Prior median", y = "MAP estimate", colour = NULL,
         title = "Parameter estimates versus prior medians") +
    theme_bw()
}

plot_nm_fit <- function(fit) {
  pred <- fit$prediction
  n_obs <- sum(obs_tot)
  pop_model <- pred$p_age * n_obs
  dat <- rbind(
    data.frame(age = seq_len(N_AGE), outcome = "Population count",
               series = "Observed", value = obs_tot),
    data.frame(age = seq_len(N_AGE), outcome = "Population count",
               series = "Fitted", value = pop_model),
    data.frame(age = seq_len(N_AGE), outcome = "HCV prevalence",
               series = "Observed", value = obs_prev),
    data.frame(age = seq_len(N_AGE), outcome = "HCV prevalence",
               series = "Fitted", value = pred$q_age)
  )
  ggplot(dat, aes(age, value, colour = series, group = series)) +
    geom_line(linewidth = 0.8) + geom_point(size = 2) +
    facet_wrap(~outcome, scales = "free_y", ncol = 1) +
    scale_x_continuous(breaks = seq_len(N_AGE)) +
    labs(x = "Age group", y = NULL, colour = NULL, title = "Observed and fitted outcomes") +
    theme_bw()
}

plot_nm_residuals <- function(fit) {
  pred <- fit$prediction
  # Binomial Pearson residual: (obs_pos - N * q) / sqrt(N * q * (1 - q))
  q_prev <- pmin(pmax(pred$q_age, 1e-12), 1 - 1e-12)
  prev_resid <- (obs_pos - obs_tot * q_prev) / sqrt(obs_tot * q_prev * (1 - q_prev))
  pop_sd <- if (length(fit$sim_data$sigma_pop) == 1L) {
    rep(fit$sim_data$sigma_pop, N_AGE)
  } else fit$sim_data$sigma_pop
  q_obs <- obs_tot / sum(obs_tot)
  ref <- N_AGE
  pop_resid <- rep(NA_real_, N_AGE)
  pop_resid[-ref] <- (
    log(q_obs[-ref] / q_obs[ref]) - log(pred$p_age[-ref] / pred$p_age[ref])
  ) / pop_sd[-ref]
  dat <- rbind(
    data.frame(age = seq_len(N_AGE), component = "Prevalence",
               residual = prev_resid),
    data.frame(age = seq_len(N_AGE), component = "Population ALR",
               residual = pop_resid)
  )
  ggplot(dat, aes(age, residual, colour = component)) +
    geom_hline(yintercept = 0, colour = "grey45") +
    geom_hline(yintercept = c(-2, 2), linetype = 2, colour = "grey65") +
    geom_segment(aes(xend = age, y = 0, yend = residual), linewidth = 0.7) +
    geom_point(size = 2.5) +
    facet_wrap(~component, ncol = 1) +
    scale_x_continuous(breaks = seq_len(N_AGE)) +
    labs(x = "Age group", y = "Standardized residual", colour = NULL,
         title = "Fitting residuals") +
    theme_bw() + theme(legend.position = "none")
}

save_nm_plots <- function(fit, directory = file.path("two-steps-calibration", "nelder_mead_plots"), width = 8, height = 6) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  plots <- list(
    convergence = plot_nm_convergence(fit),
    parameters = plot_nm_parameters(fit),
    fitted = plot_nm_fit(fit),
    residuals = plot_nm_residuals(fit)
  )
  for (nm in names(plots)) {
    ggsave(file.path(directory, paste0(nm, ".png")), plots[[nm]],
           width = width, height = height, dpi = 300)
  }
  invisible(plots)
}

if (sys.nframe() == 0L) {
  out_dir <- "two-steps-calibration"
  args <- commandArgs(trailingOnly = TRUE)
  out_arg <- grep("^--out-dir=", args, value = TRUE)
  if (length(out_arg) > 0L) out_dir <- sub("^--out-dir=", "", out_arg[[1L]])
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  fit <- fit_nelder_mead()
  print(fit)
  saveRDS(fit, file.path(out_dir, "nelder_mead_fit.rds"))
  save_nm_plots(fit, file.path(out_dir, "nelder_mead_plots"))
  cat(sprintf("Saved Nelder-Mead outputs in %s\n", out_dir))
}
