library(ggplot2)
# Trial
fit = readRDS("two-steps-calibration/nelder_mead_fit.rds")
fit1 = fit
contact_scale <- c(7.68983, 0.06629, 0.08773, 0.61861, 1.74367, -3.24803)
inflow_scale <- c(0.63544, 6.33889, 10.05326, 0.33554, 1.16247, 1.55939)

source('setup.R')
obs_prev <- c(
    0.1118421, 0.1731044, 0.2684954, 0.4301165, 0.4821029, 0.3544304
)
obs_prev_se <- sqrt(c(
    0.09933344, 0.14313927, 0.19640562, 0.24511630, 0.24967969, 0.22880949
))
obs_tot <- c(99, 1244, 1467, 1841, 1628, 409)
obs_pos <- round(obs_prev * obs_tot) # HCV positives per age group (binomial numerator)
N_AGE <- length(obs_tot)
sigma_pop <- c(0.10, rep(0.06, N_AGE - 1L))
pm = params
for (j in seq_along(contact_scale)) {
    pm$C_contact[j, ] <- params$C_contact[j, ] * contact_scale[j]
}

# Scaling total inflow beta
pm$beta <- params$beta * inflow_scale

c_true <- params$c_composite / inflow_scale
pm$lambda1 <- params$lambda3 * c_true

out <- run_sim(pm, data)
fit1$prediction <- compute_age_quantities(as.numeric(out[nrow(out), -1L]))
plot_nm_fit(
    fit = fit1,
    obs_tot = obs_tot,
    obs_prev = obs_prev,
    sigma_pop = sigma_pop,
    obs_pos = obs_pos,
    level = 0.95
)


# # =============================================================================
# # Laplace approximation around the Nelder-Mead MAP estimate
# # =============================================================================

# fit_laplace_approximation <- function(
#   fit,
#   n_draws = 5000L,
#   level = 0.95,
#   hessian_step = 1e-3,
#   eigen_floor_rel = 1e-8,
#   seed = 1234L
# ) {
#     if (!inherits(fit, "nm_hcv_fit")) {
#         stop("fit must be an nm_hcv_fit object")
#     }

#     if (fit$convergence != 0L) {
#         warning(
#             "Nelder-Mead did not report successful convergence. ",
#             "The Laplace approximation may not be reliable."
#         )
#     }

#     theta_hat <- as.numeric(fit$theta_hat)
#     names(theta_hat) <- names(fit$theta_hat)

#     n_par <- length(theta_hat)

#     if (n_par != 2L * SPLINE_K) {
#         stop(sprintf(
#             "Expected %d parameters but found %d",
#             2L * SPLINE_K,
#             n_par
#         ))
#     }

#     # Negative log-posterior.
#     #
#     # The Hessian must be calculated for -log posterior, not log posterior,
#     # so that a well-defined local MAP gives a positive-definite Hessian.
#     neg_log_posterior <- function(theta) {
#         lp <- tryCatch(
#             log_posterior(
#                 theta = theta,
#                 base_params = fit$base_params,
#                 data = fit$sim_data
#             ),
#             error = function(e) -Inf
#         )

#         if (!is.finite(lp)) {
#             return(.Machine$double.xmax / 100)
#         }

#         -lp
#     }

#     # Numerical differentiation step for each parameter.
#     #
#     # Second derivatives are usually less numerically stable than first
#     # derivatives, so 1e-3 is often safer than 1e-4 for an ODE model.
#     ndeps <- hessian_step * pmax(abs(theta_hat), 1)

#     cat("Calculating numerical Hessian at the MAP estimate...\n")

#     hessian_raw <- optimHess(
#         par = theta_hat,
#         fn = neg_log_posterior,
#         control = list(ndeps = ndeps)
#     )

#     if (any(!is.finite(hessian_raw))) {
#         stop("The numerical Hessian contains non-finite values")
#     }

#     # Numerical differentiation may produce slight asymmetry.
#     hessian_raw <- 0.5 * (hessian_raw + t(hessian_raw))

#     rownames(hessian_raw) <- names(theta_hat)
#     colnames(hessian_raw) <- names(theta_hat)

#     eig <- eigen(hessian_raw, symmetric = TRUE)

#     raw_eigenvalues <- eig$values
#     max_abs_eigenvalue <- max(abs(raw_eigenvalues))

#     eigen_floor <- max(
#         eigen_floor_rel * max_abs_eigenvalue,
#         sqrt(.Machine$double.eps)
#     )

#     adjusted_eigenvalues <- pmax(
#         raw_eigenvalues,
#         eigen_floor
#     )

#     n_adjusted <- sum(raw_eigenvalues < eigen_floor)

#     if (n_adjusted > 0L) {
#         warning(sprintf(
#             paste0(
#                 "%d Hessian eigenvalues were below %.3e and were regularized. ",
#                 "Inspect the Hessian diagnostics carefully."
#             ),
#             n_adjusted,
#             eigen_floor
#         ))
#     }

#     # Positive-definite approximation to the Hessian.
#     hessian_pd <-
#         eig$vectors %*%
#         diag(adjusted_eigenvalues, nrow = n_par) %*%
#         t(eig$vectors)

#     # Laplace covariance: inverse Hessian.
#     covariance <-
#         eig$vectors %*%
#         diag(1 / adjusted_eigenvalues, nrow = n_par) %*%
#         t(eig$vectors)

#     covariance <- 0.5 * (covariance + t(covariance))

#     rownames(covariance) <- names(theta_hat)
#     colnames(covariance) <- names(theta_hat)

#     standard_errors <- sqrt(diag(covariance))

#     alpha <- 1 - level
#     z_critical <- qnorm(1 - alpha / 2)

#     coefficient_summary <- data.frame(
#         parameter = names(theta_hat),
#         estimate = theta_hat,
#         standard_error = standard_errors,
#         lower = theta_hat - z_critical * standard_errors,
#         upper = theta_hat + z_critical * standard_errors,
#         row.names = NULL
#     )

#     # Draw from N(theta_hat, covariance).
#     #
#     # If R = chol(covariance), then Z %*% R has covariance R'R = covariance.
#     set.seed(seed)

#     chol_covariance <- tryCatch(
#         chol(covariance),
#         error = function(e) {
#             stop(
#                 "Adjusted covariance matrix could not be Cholesky decomposed: ",
#                 conditionMessage(e)
#             )
#         }
#     )

#     z_draws <- matrix(
#         rnorm(n_draws * n_par),
#         nrow = n_draws,
#         ncol = n_par
#     )

#     theta_draws <- sweep(
#         z_draws %*% chol_covariance,
#         MARGIN = 2L,
#         STATS = theta_hat,
#         FUN = "+"
#     )

#     colnames(theta_draws) <- names(theta_hat)

#     diagnostics <- list(
#         raw_eigenvalues = raw_eigenvalues,
#         adjusted_eigenvalues = adjusted_eigenvalues,
#         minimum_raw_eigenvalue = min(raw_eigenvalues),
#         maximum_raw_eigenvalue = max(raw_eigenvalues),
#         n_adjusted_eigenvalues = n_adjusted,
#         eigenvalue_floor = eigen_floor,
#         adjusted_condition_number =
#             max(adjusted_eigenvalues) / min(adjusted_eigenvalues),
#         negative_log_posterior_at_map =
#             neg_log_posterior(theta_hat),
#         hessian_step = hessian_step
#     )

#     result <- list(
#         theta_hat = theta_hat,
#         hessian_raw = hessian_raw,
#         hessian_pd = hessian_pd,
#         covariance = covariance,
#         standard_errors = standard_errors,
#         coefficient_summary = coefficient_summary,
#         theta_draws = theta_draws,
#         diagnostics = diagnostics,
#         level = level
#     )

#     class(result) <- "laplace_hcv_fit"
#     result
# }

# laplace_fit <- fit_laplace_approximation(
#     fit = fit,
#     n_draws = 5000L,
#     level = 0.95,
#     hessian_step = 1e-3,
#     seed = 1234L
# )

# summarise_laplace_age_parameters <- function(
#   laplace_fit,
#   level = 0.95
# ) {
#     theta_draws <- laplace_fit$theta_draws

#     age_parameter_draws <- theta_to_orig(theta_draws)

#     alpha <- 1 - level

#     probs <- c(
#         alpha / 2,
#         0.50,
#         1 - alpha / 2
#     )

#     q <- apply(
#         age_parameter_draws,
#         MARGIN = 2L,
#         FUN = quantile,
#         probs = probs,
#         na.rm = TRUE
#     )

#     data.frame(
#         parameter = colnames(age_parameter_draws),
#         median = q[2L, ],
#         lower = q[1L, ],
#         upper = q[3L, ],
#         row.names = NULL
#     )
# }

# generate_laplace_fitted_draws <- function(
#   fit,
#   laplace_fit,
#   n_prediction_draws = 1000L,
#   seed = 4321L,
#   report_every = 100L
# ) {
#     theta_draws <- as.matrix(laplace_fit$theta_draws)

#     if (n_prediction_draws < 1L) {
#         stop("n_prediction_draws must be at least 1")
#     }

#     set.seed(seed)

#     if (nrow(theta_draws) > n_prediction_draws) {
#         selected_rows <- sort(sample(
#             seq_len(nrow(theta_draws)),
#             size = n_prediction_draws,
#             replace = FALSE
#         ))

#         theta_used <- theta_draws[
#             selected_rows, ,
#             drop = FALSE
#         ]
#     } else {
#         selected_rows <- seq_len(nrow(theta_draws))
#         theta_used <- theta_draws
#     }

#     n_used <- nrow(theta_used)

#     p_age_draws <- matrix(
#         NA_real_,
#         nrow = n_used,
#         ncol = N_AGE
#     )

#     q_age_draws <- matrix(
#         NA_real_,
#         nrow = n_used,
#         ncol = N_AGE
#     )

#     valid <- logical(n_used)

#     for (s in seq_len(n_used)) {
#         pred_s <- tryCatch(
#             predict_at_theta(
#                 theta = as.numeric(theta_used[s, ]),
#                 base_params = fit$base_params,
#                 sim_data = fit$sim_data
#             ),
#             error = function(e) NULL
#         )

#         if (!is.null(pred_s) &&
#             all(is.finite(pred_s$p_age)) &&
#             all(is.finite(pred_s$q_age))) {
#             p_age_draws[s, ] <- pred_s$p_age
#             q_age_draws[s, ] <- pred_s$q_age
#             valid[s] <- TRUE
#         }

#         if (report_every > 0L &&
#             s %% report_every == 0L) {
#             cat(sprintf(
#                 "Laplace ODE predictions: %d/%d completed\n",
#                 s,
#                 n_used
#             ))
#         }
#     }

#     if (!any(valid)) {
#         stop("All Laplace ODE prediction runs failed")
#     }

#     if (any(!valid)) {
#         warning(sprintf(
#             "%d of %d Laplace draws produced invalid ODE predictions",
#             sum(!valid),
#             n_used
#         ))
#     }

#     list(
#         theta = theta_used[valid, , drop = FALSE],
#         p_age = p_age_draws[valid, , drop = FALSE],
#         q_age = q_age_draws[valid, , drop = FALSE],
#         selected_rows = selected_rows[valid],
#         n_attempted = n_used,
#         n_valid = sum(valid)
#     )
# }

# laplace_fitted_draws <-
#     generate_laplace_fitted_draws(
#         fit = fit,
#         laplace_fit = laplace_fit,
#         n_prediction_draws = 1000L,
#         seed = 4321L
#     )


# summarise_laplace_fitted_values <- function(
#   fit,
#   laplace_fitted_draws,
#   obs_tot,
#   level = 0.95
# ) {
#     alpha <- 1 - level

#     probs <- c(
#         alpha / 2,
#         0.50,
#         1 - alpha / 2
#     )

#     population_draws <-
#         laplace_fitted_draws$p_age * sum(obs_tot)

#     prevalence_draws <-
#         laplace_fitted_draws$q_age

#     population_quantiles <- apply(
#         population_draws,
#         MARGIN = 2L,
#         FUN = quantile,
#         probs = probs,
#         na.rm = TRUE
#     )

#     prevalence_quantiles <- apply(
#         prevalence_draws,
#         MARGIN = 2L,
#         FUN = quantile,
#         probs = probs,
#         na.rm = TRUE
#     )

#     map_population <-
#         fit$prediction$p_age * sum(obs_tot)

#     map_prevalence <-
#         fit$prediction$q_age

#     rbind(
#         data.frame(
#             age = seq_len(N_AGE),
#             outcome = "Population count",
#             map_fitted = map_population,
#             median = population_quantiles[2L, ],
#             lower = population_quantiles[1L, ],
#             upper = population_quantiles[3L, ]
#         ),
#         data.frame(
#             age = seq_len(N_AGE),
#             outcome = "HCV prevalence",
#             map_fitted = map_prevalence,
#             median = prevalence_quantiles[2L, ],
#             lower = prevalence_quantiles[1L, ],
#             upper = prevalence_quantiles[3L, ]
#         )
#     )
# }

# laplace_fitted_summary <-
#     summarise_laplace_fitted_values(
#         fit = fit,
#         laplace_fitted_draws = laplace_fitted_draws,
#         obs_tot = obs_tot,
#         level = 0.95
#     )

# print(laplace_fitted_summary)