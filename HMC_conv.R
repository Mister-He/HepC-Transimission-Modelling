# =============================================================================
# 6.  CONVERGENCE DIAGNOSTICS: R-HAT AND ESS
# =============================================================================

# ── 6a. Gelman-Rubin R-hat ───────────────────────────────────────────────────
#
# R-hat measures mixing across chains.  R-hat < 1.01 indicates convergence
# (Vehtari et al. 2021); values 1.01–1.05 are borderline.
#
# Formula (split-chain pooled variance):
#   B = (N/(M-1)) Σ_m (ψ̄_m - ψ̄)²    between-chain variance
#   W = (1/M)     Σ_m s²_m             within-chain variance
#   var̂(ψ) = (N-1)/N · W + B/N
#   R-hat   = sqrt(var̂(ψ) / W)
#
compute_rhat <- function(chains) {
    M <- length(chains)
    N <- nrow(chains[[1L]])
    J <- ncol(chains[[1L]])

    rhat <- setNames(numeric(J), param_names_log)

    for (j in seq_len(J)) {
        chain_j <- lapply(chains, function(ch) ch[, j])
        chain_means <- vapply(chain_j, mean, numeric(1L))
        chain_vars <- vapply(chain_j, var, numeric(1L))

        grand_mean <- mean(chain_means)
        B <- N / (M - 1) * sum((chain_means - grand_mean)^2)
        W <- mean(chain_vars)
        var_hat <- (N - 1) / N * W + B / N

        rhat[j] <- sqrt(var_hat / W)
    }
    rhat
}

# ── 6b. Effective Sample Size (ESS) ──────────────────────────────────────────
#
# ESS based on Geyer's initial monotone sequence estimator applied to the
# pooled autocorrelation function across all chains.
#
# ESS = (M × N) / τ,   where τ = 1 + 2 Σ_{t=1}^{T*} ρ_t
# T* = last lag where ρ_{2t} + ρ_{2t+1} ≥ 0 (monotone pairs truncation).
#
# Target: ESS > 400, or ESS/N_total > 0.1 per parameter.
#
compute_ess <- function(chains) {
    M <- length(chains)
    N <- nrow(chains[[1L]])
    J <- ncol(chains[[1L]])

    ess <- setNames(numeric(J), param_names_log)

    for (j in seq_len(J)) {
        all_samps <- unlist(lapply(chains, function(ch) ch[, j]))
        n_total <- length(all_samps)
        max_lag <- min(n_total - 1L, 500L)

        acf_vals <- acf(all_samps,
            lag.max = max_lag,
            plot = FALSE, type = "correlation"
        )$acf[-1L]

        # Geyer monotone pairs:  Γ_t = ρ_{2t} + ρ_{2t+1}
        n_pairs <- floor(length(acf_vals) / 2L)
        if (n_pairs < 1L) {
            ess[j] <- n_total
            next
        }

        pairs <- acf_vals[2L * seq_len(n_pairs) - 1L] +
            acf_vals[2L * seq_len(n_pairs)]
        cutoff <- which(pairs <= 0)[1L]
        if (is.na(cutoff)) cutoff <- n_pairs

        # Sum autocorrelations up to cut-off pair
        rho_sum <- 1 + 2 * sum(acf_vals[seq_len(2L * (cutoff - 1L) + 1L)])
        ess[j] <- n_total / max(rho_sum, 1)
    }
    ess
}

# ── 6c. Divergence check ─────────────────────────────────────────────────────
#
# A "divergence" occurs when the leapfrog trajectory explodes (H changes by
# > threshold).  Here we flag iterations where |ΔH| > 1000 as divergent.
# In practice, this is visible as sudden jumps in the lp trace.
#
count_divergences <- function(chains_raw, threshold = 1000) {
    lapply(chains_raw, function(ch) {
        lp <- ch$lp_trace[ch$n_warmup + seq_len(nrow(ch$samples))]
        sum(abs(diff(lp)) > threshold, na.rm = TRUE)
    })
}

# ── 6d. Print formatted diagnostics table ───────────────────────────────────
print_diagnostics <- function(post_warmup_list, chains_raw = NULL) {
    rhat_vals <- compute_rhat(post_warmup_list)
    ess_vals <- compute_ess(post_warmup_list)
    all_samps <- do.call(rbind, post_warmup_list)
    orig_samps <- exp(all_samps)

    diag_df <- data.frame(
        Parameter = param_names_log,
        Mean_log = round(colMeans(all_samps), 4),
        SD_log = round(apply(all_samps, 2, sd), 4),
        Mean_orig = round(colMeans(orig_samps), 4),
        Q2.5_orig = round(apply(orig_samps, 2, quantile, 0.025), 4),
        Q97.5_orig = round(apply(orig_samps, 2, quantile, 0.975), 4),
        Rhat = round(rhat_vals, 4),
        ESS = round(ess_vals, 1),
        Conv_OK = ifelse(rhat_vals < 1.05 & ess_vals > 100, "YES", "NO"),
        stringsAsFactors = FALSE
    )

    cat("\n======================================================\n")
    cat("          Convergence Diagnostics Summary\n")
    cat("  Target: R-hat < 1.05   |   ESS > 100 per parameter\n")
    cat("======================================================\n")
    print(diag_df, row.names = FALSE)

    if (!is.null(chains_raw)) {
        div_counts <- count_divergences(chains_raw)
        cat("\nDivergences per chain (post-warmup):\n")
        print(setNames(unlist(div_counts), paste0("Chain ", seq_along(div_counts))))
    }

    invisible(diag_df)
}


# =============================================================================
# 7.  POSTERIOR PREDICTIVE CHECKS (PPC)
# =============================================================================

# ── 7a. Generate PPC replicates ──────────────────────────────────────────────
#
# For each sampled theta_s, run the ODE model and draw:
#   y_rep_pos[i] ~ Poisson(lambda_pos[i] | theta_s)
#   y_rep_tot[i] ~ Poisson(lambda_tot[i] | theta_s)
#
# Returns replicated counts (integer) and the Poisson means (real),
# plus posterior predictive p-values (ppp) per age group.
#
generate_ppc_samples <- function(post_samples, base_params, data,
                                 n_ppc = 600L) {
    n_avail <- nrow(post_samples)
    draw_idx <- sort(sample(n_avail, min(n_ppc, n_avail)))
    n_draws <- length(draw_idx)

    ppc_pos <- matrix(NA_integer_, n_draws, 9L)
    ppc_tot <- matrix(NA_integer_, n_draws, 9L)
    lam_pos <- matrix(NA_real_, n_draws, 9L)
    lam_tot <- matrix(NA_real_, n_draws, 9L)

    cat(sprintf(
        "Generating %d PPC draws from %d posterior samples...\n",
        n_draws, n_avail
    ))

    for (ii in seq_len(n_draws)) {
        theta <- post_samples[draw_idx[ii], ]

        result <- tryCatch(
            {
                pm <- build_params_from_theta(theta, base_params)
                out <- run_sim(pm, data)
                y_fin <- as.numeric(out[nrow(out), -1L])
                compute_age_quantities(y_fin)
            },
            error = function(e) NULL
        )

        if (!is.null(result)) {
            sigma_N <- if (!is.null(data$sigma_N)) data$sigma_N else 0.10
            phi_overdisp <- if (!is.null(data$phi_overdisp)) data$phi_overdisp else 50.0

            # Draw a replicated total count from the LogNormal observation model
            N_rep <- as.integer(round(stats::rlnorm(1L, meanlog = log(result$n_model_total), sdlog = sigma_N)))
            if (!is.finite(N_rep) || N_rep < 1L) N_rep <- 1L

            # Poisson means (expected counts) under this replicate total
            lam_tot[ii, ] <- N_rep * result$p_age
            lam_pos[ii, ] <- lam_tot[ii, ] * result$q_age

            # Draw multinomial totals conditional on the replicated total
            ppc_tot[ii, ] <- as.integer(rmultinom(1L, size = N_rep, prob = result$p_age))

            # For positives, incorporate Beta-Binomial overdispersion via a beta draw
            p_draw <- rbeta(9L, shape1 = result$q_age * phi_overdisp, shape2 = (1 - result$q_age) * phi_overdisp)
            ppc_pos[ii, ] <- rbinom(9L, size = ppc_tot[ii, ], prob = p_draw)
        }

        if (ii %% 50L == 0L) {
            cat(sprintf("  PPC draw %d / %d\n", ii, n_draws))
        }
    }

    # Posterior predictive p-values: Pr(y_rep >= y_obs | y_obs)
    # Calibrated model: ppp near 0.5; poor fit: near 0 or 1.
    ppp_pos <- vapply(1:9, function(i) {
        mean(ppc_pos[, i] >= obs_pos[i], na.rm = TRUE)
    }, numeric(1L))
    ppp_tot <- vapply(1:9, function(i) {
        mean(ppc_tot[, i] >= obs_tot[i], na.rm = TRUE)
    }, numeric(1L))

    list(
        ppc_pos = ppc_pos, ppc_tot = ppc_tot,
        lam_pos = lam_pos, lam_tot = lam_tot,
        ppp_pos = ppp_pos, ppp_tot = ppp_tot
    )
}

# ── 7b. PPC histogram plot ───────────────────────────────────────────────────
#
# One panel per age group: histogram of y_rep (9 × n_ppc), vertical line at
# observed count, and posterior predictive p-value annotated per panel.
#
plot_ppc_histograms <- function(ppc, type = c("pos", "tot")) {
    type <- match.arg(type)

    ppc_mat <- if (type == "pos") ppc$ppc_pos else ppc$ppc_tot
    obs_vals <- if (type == "pos") obs_pos else obs_tot
    ppp_vals <- if (type == "pos") ppc$ppp_pos else ppc$ppp_tot

    fill_col <- if (type == "pos") "#2166AC" else "#D6604D"
    main_ttl <- if (type == "pos") {
        "PPC: HCV-Positive Counts (J-stratum)"
    } else {
        "PPC: Total PWID Counts (J-stratum)"
    }

    age_lbl <- paste0("Age group ", 1:9)

    # Drop failed draws
    ppc_mat <- ppc_mat[complete.cases(ppc_mat), , drop = FALSE]

    ppc_long <- as.data.frame(ppc_mat) %>%
        setNames(age_lbl) %>%
        pivot_longer(everything(), names_to = "Age", values_to = "y_rep") %>%
        mutate(Age = factor(Age, levels = age_lbl))

    obs_df <- data.frame(
        Age = factor(age_lbl, levels = age_lbl),
        obs = obs_vals
    )
    ppp_df <- data.frame(
        Age = factor(age_lbl, levels = age_lbl),
        ppp = round(ppp_vals, 2)
    )

    ggplot(ppc_long, aes(x = y_rep)) +
        geom_histogram(aes(y = after_stat(density)),
            bins = 30,
            fill = fill_col,
            alpha = 0.65,
            colour = "white",
            linewidth = 0.2
        ) +
        geom_vline(
            data = obs_df, aes(xintercept = obs),
            colour = "black", linewidth = 1.1
        ) +
        geom_text(
            data = ppp_df,
            aes(label = paste0("ppp=", ppp)),
            x = Inf, y = Inf, hjust = 1.15, vjust = 1.6,
            size = 3.0, colour = "grey25"
        ) +
        facet_wrap(~Age, scales = "free", ncol = 3) +
        labs(
            title    = main_ttl,
            subtitle = "Histogram: posterior predictive  |  Vertical line: observed  |  ppp near 0.5 = good fit",
            x        = "Count",
            y        = "Density"
        ) +
        theme_bw(base_size = 11) +
        theme(
            strip.background  = element_rect(fill = "grey92", colour = "grey70"),
            panel.grid.minor  = element_blank(),
            plot.subtitle     = element_text(size = 9, colour = "grey40")
        )
}

# ── 7c. PPC interval plot ─────────────────────────────────────────────────────
#
# Dot-and-whisker: posterior predictive median ± 50% and 90% intervals per age
# group, with observed count overlaid as a red point.
#
plot_ppc_intervals <- function(ppc) {
    age_lbl <- paste0("Age ", 1:9)
    z_crit <- stats::qnorm(0.975)

    summarise_ppc_mat <- function(ppc_mat, obs, type_label) {
        mat <- ppc_mat[complete.cases(ppc_mat), , drop = FALSE]
        data.frame(
            Age = factor(age_lbl, levels = age_lbl),
            obs = obs,
            med = apply(mat, 2, median),
            lo50 = apply(mat, 2, quantile, 0.25),
            hi50 = apply(mat, 2, quantile, 0.75),
            lo95 = apply(mat, 2, quantile, 0.025),
            hi95 = apply(mat, 2, quantile, 0.975),
            type = type_label,
            stringsAsFactors = FALSE
        )
    }

    summarise_obs_ci <- function(ppc_mat, obs, type_label) {
        keep <- complete.cases(ppc_mat)

        if (!any(keep)) {
            return(data.frame(
                Age = factor(age_lbl, levels = age_lbl),
                obs_lo = rep(NA_real_, length(obs)),
                obs_hi = rep(NA_real_, length(obs)),
                type = type_label,
                stringsAsFactors = FALSE
            ))
        }

        if (type_label == "HCV-positive") {
            q_mat <- ppc$lam_pos[keep, , drop = FALSE] / ppc$lam_tot[keep, , drop = FALSE]
            sd_vals <- vapply(seq_len(ncol(q_mat)), function(j) {
                q_draws <- q_mat[, j]
                q_draws <- q_draws[is.finite(q_draws)]
                if (!length(q_draws)) {
                    return(NA_real_)
                }
                stats::median(
                    sqrt(obs[j] * q_draws * (1 - q_draws) * (obs[j] + phi_overdisp) / (1 + phi_overdisp)),
                    na.rm = TRUE
                )
            }, numeric(1L))
        } else {
            p_mat <- ppc$lam_tot[keep, , drop = FALSE] / N_total_obs
            sd_vals <- vapply(seq_len(ncol(p_mat)), function(j) {
                p_draws <- p_mat[, j]
                p_draws <- p_draws[is.finite(p_draws)]
                if (!length(p_draws)) {
                    return(NA_real_)
                }
                stats::median(sqrt(N_total_obs * p_draws * (1 - p_draws)), na.rm = TRUE)
            }, numeric(1L))
        }

        data.frame(
            Age = factor(age_lbl, levels = age_lbl),
            obs_lo = pmax(0, obs - z_crit * sd_vals),
            obs_hi = obs + z_crit * sd_vals,
            type = type_label,
            stringsAsFactors = FALSE
        )
    }

    plot_df <- bind_rows(
        summarise_ppc_mat(ppc$ppc_pos, obs_pos, "HCV-positive"),
        summarise_ppc_mat(ppc$ppc_tot, obs_tot, "Total PWID")
    ) %>%
        left_join(
            bind_rows(
                summarise_obs_ci(ppc$ppc_pos, obs_pos, "HCV-positive"),
                summarise_obs_ci(ppc$ppc_tot, obs_tot, "Total PWID")
            ),
            by = c("Age", "type")
        )

    ggplot(plot_df, aes(x = Age)) +
        geom_linerange(aes(ymin = lo95, ymax = hi95),
            linewidth = 1.0, colour = "#4292C6", alpha = 0.45
        ) +
        geom_linerange(aes(ymin = lo50, ymax = hi50),
            linewidth = 3.0, colour = "#08519C", alpha = 0.75
        ) +
        geom_linerange(aes(ymin = obs_lo, ymax = obs_hi),
            linewidth = 0.9, colour = "#D7191C", alpha = 0.75
        ) +
        geom_point(aes(y = med),
            shape = 18, size = 3.5, colour = "#08306B"
        ) +
        geom_point(aes(y = obs),
            shape = 21, size = 3.5,
            fill = "#D7191C", colour = "black", stroke = 1.2
        ) +
        facet_wrap(~type, scales = "free_y", ncol = 1) +
        labs(
            title    = "PPC: Posterior Predictive Intervals vs Observed",
            subtitle = "Diamond: predictive median  |  Thick bar: 50% PI  |  Thin bar: 90% PI  |  Red whisker: likelihood-based 95% CI  |  Red dot: observed",
            x        = "Age group",
            y        = "Count"
        ) +
        theme_bw(base_size = 11) +
        theme(
            axis.text.x       = element_text(angle = 35, hjust = 1),
            strip.background  = element_rect(fill = "grey92", colour = "grey70"),
            panel.grid.minor  = element_blank(),
            plot.subtitle     = element_text(size = 9, colour = "grey40")
        )
}

# ── 7d. HCV prevalence interval plot ────────────────────────────────────────
#
# Age-group prevalence = HCV-positive / total PWID in the J-stratum.
# Blue interval: posterior predictive prevalence draws (BCI).
# Red interval: observed binomial confidence interval.
#
plot_ppc_prevalence_intervals <- function(ppc) {
    age_lbl <- paste0("Age ", 1:9)

    keep <- complete.cases(ppc$ppc_pos, ppc$ppc_tot)
    if (!any(keep)) {
        stop("No complete PPC draws available for prevalence plotting")
    }

    prev_mat <- ppc$ppc_pos[keep, , drop = FALSE] / ppc$ppc_tot[keep, , drop = FALSE]

    sample_df <- data.frame(
        Age = factor(age_lbl, levels = age_lbl),
        obs = obs_pos / obs_tot,
        med = apply(prev_mat, 2, median, na.rm = TRUE),
        lo95 = apply(prev_mat, 2, quantile, 0.025, na.rm = TRUE),
        hi95 = apply(prev_mat, 2, quantile, 0.975, na.rm = TRUE),
        type = "HCV prevalence",
        stringsAsFactors = FALSE
    )

    obs_ci <- vapply(seq_along(obs_pos), function(i) {
        stats::binom.test(obs_pos[i], obs_tot[i])$conf.int
    }, numeric(2L))

    obs_df <- data.frame(
        Age = factor(age_lbl, levels = age_lbl),
        obs_lo = obs_ci[1L, ],
        obs_hi = obs_ci[2L, ],
        type = "HCV prevalence",
        stringsAsFactors = FALSE
    )

    plot_df <- left_join(sample_df, obs_df, by = c("Age", "type"))

    ggplot(plot_df, aes(x = Age)) +
        geom_linerange(aes(ymin = lo95, ymax = hi95),
            linewidth = 1.1, colour = "#08519C", alpha = 0.55
        ) +
        geom_point(aes(y = med),
            shape = 18, size = 3.5, colour = "#08306B"
        ) +
        geom_linerange(aes(ymin = obs_lo, ymax = obs_hi),
            linewidth = 0.9, colour = "#D7191C", alpha = 0.75
        ) +
        geom_point(aes(y = obs),
            shape = 21, size = 3.5,
            fill = "#D7191C", colour = "black", stroke = 1.2
        ) +
        labs(
            title    = "PPC: Age-group HCV prevalence",
            subtitle = "Blue interval: posterior predictive 95% BCI  |  Red whisker: observed exact 95% CI  |  Diamond: sample median  |  Red dot: observed prevalence",
            x        = "Age group",
            y        = "Prevalence"
        ) +
        scale_y_continuous(labels = function(x) sprintf("%.1f%%", 100 * x)) +
        theme_bw(base_size = 11) +
        theme(
            strip.background  = element_rect(fill = "grey92", colour = "grey70"),
            panel.grid.minor  = element_blank(),
            plot.subtitle     = element_text(size = 9, colour = "grey40")
        )
}

# ── 7d. Trace plots ───────────────────────────────────────────────────────────
#
# Full-chain trace for selected parameters (all iterations, all chains),
# with a dashed vertical line marking the end of warmup.
#
plot_traces <- function(chains_raw, param_idx = 1:min(5L, N_PARAMS)) {
    trace_list <- lapply(seq_along(chains_raw), function(ch) {
        df <- as.data.frame(chains_raw[[ch]]$samples_all)
        colnames(df) <- param_names_log
        df$iter <- seq_len(nrow(df))
        df$chain <- factor(ch)
        df
    })
    trace_df <- bind_rows(trace_list)
    n_warmup <- chains_raw[[1L]]$n_warmup
    params_to <- param_names_log[param_idx]

    trace_df %>%
        select(iter, chain, all_of(params_to)) %>%
        pivot_longer(all_of(params_to),
            names_to  = "param",
            values_to = "value"
        ) %>%
        mutate(param = factor(param, levels = params_to)) %>%
        ggplot(aes(x = iter, y = value, colour = chain)) +
        geom_line(alpha = 0.70, linewidth = 0.35) +
        geom_vline(
            xintercept = n_warmup,
            linetype = "dashed", colour = "grey35", linewidth = 0.8
        ) +
        facet_wrap(~param, scales = "free_y", ncol = 1) +
        scale_colour_brewer(palette = "Set1") +
        labs(
            title    = "Trace plots — all chains (log scale)",
            subtitle = "Dashed line = end of warmup",
            x        = "Iteration",
            y        = "Value (log-transformed parameter)",
            colour   = "Chain"
        ) +
        theme_bw(base_size = 11) +
        theme(
            strip.background  = element_rect(fill = "grey92", colour = "grey70"),
            panel.grid.minor  = element_blank(),
            legend.position   = "right",
            plot.subtitle     = element_text(size = 9, colour = "grey40")
        )
}

# ── 7e. Posterior density plot (original scale) ────────────────────────────
plot_posterior_densities <- function(post_warmup_list) {
    all_samps <- do.call(rbind, lapply(seq_along(post_warmup_list), function(ch) {
        df <- as.data.frame(exp(post_warmup_list[[ch]]))
        colnames(df) <- param_names_orig
        df$chain <- factor(ch)
        df
    }))

    all_samps %>%
        pivot_longer(-chain, names_to = "param", values_to = "value") %>%
        mutate(param = factor(param, levels = param_names_orig)) %>%
        ggplot(aes(x = value, fill = chain, colour = chain)) +
        geom_density(alpha = 0.30, linewidth = 0.6) +
        facet_wrap(~param, scales = "free", ncol = 2) +
        scale_fill_brewer(palette = "Set1") +
        scale_colour_brewer(palette = "Set1") +
        labs(
            title    = "Posterior densities (original scale, all chains)",
            x        = "Parameter value",
            y        = "Density",
            fill     = "Chain",
            colour   = "Chain"
        ) +
        theme_bw(base_size = 11) +
        theme(
            strip.background = element_rect(fill = "grey92", colour = "grey70"),
            panel.grid.minor = element_blank(),
            legend.position  = "right"
        )
}
