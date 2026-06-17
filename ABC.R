# =============================================================================
# MAIN: ABC CALIBRATION (rewritten)
#
# Design rationale (see accompanying notes):
#   - run_sim() is a DETERMINISTIC ODE solver: same theta -> same trajectory.
#   - The binomial / beta-binomial observation LINK used in HMC_core.R is
#     exactly the thing we do NOT want to assume here. So this script never
#     calls dmultinom(), log_beta_binomial(), qbinom(), or qbeta() on the
#     path that decides acceptance.
#   - Because run_sim() is deterministic, there is nothing to "simulate
#     stochastically" — one call per theta gives the exact model summary
#     statistics. The only legitimate likelihood-free comparison is:
#         distance( summary_stats(model | theta), summary_stats(obs) ) <= eps
#   - Summary statistics are standardized (robust scale from prior predictive
#     draws) before computing the distance, so age-totals (counts ~200-800)
#     don't dominate prevalences (proportions in [0,1]) for free.
#   - eps is chosen via an explicit, inspectable quantile of the prior
#     predictive distance distribution (a standard ABC tuning device), not
#     hard-coded, and the full distance distribution is saved so you can
#     re-tune without rerunning simulations.
#   - A hard stop (not just a warning) fires if too few samples are accepted,
#     so downstream PPC/plotting code never silently runs on a degenerate or
#     empty sample.
# =============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)

# ── Observations ─────────────────────────────────────────────────────────────
obs_pos <- c(11, 51, 99, 141, 209, 339, 437, 367, 351)
obs_tot <- c(223, 572, 790, 765, 747, 810, 770, 658, 803)
obs_prev <- obs_pos / obs_tot

param_names_log <- c(
    "mu_hier", "log_sigma_hier",
    paste0("eta_", 1:8)
)
param_names_orig <- c(
    "mu_hier", "sigma_hier",
    paste0("C_contact_scale_", 1:8)
)
N_PARAMS <- length(param_names_log)

source("setup.R")
source("HMC_core.R") # reused as-is: build_params_from_theta, theta_to_orig,
# compute_age_quantities (all link-agnostic; do NOT
# reuse log_likelihood/log_posterior, which assume the
# observation link we're avoiding here)

# ── ABC settings ──────────────────────────────────────────────────────────────
N_PRIOR_PRED <- 3000L # prior predictive draws used ONLY to set the scale
# of each summary statistic and to calibrate eps
N_ACCEPT_TARGET <- 4000L
N_PROPOSAL_MAX <- 500000L
N_PPC <- 600L
ABC_DIR <- "abc"
EPS_QUANTILE <- 0.01 # accept the closest 1% of prior predictive distances
# by default; override EPS_FIXED below to set a
# specific distance threshold instead
EPS_FIXED <- NA_real_ # set to a numeric value to bypass quantile-based eps

set.seed(114514)
dir.create(ABC_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Prior proposal ────────────────────────────────────────────────────────────
# Same prior as HMC_core.R, used here as the ABC proposal distribution.
sample_theta_prior <- function() {
    c(
        rnorm(1L, 0.0, 1.0), # mu_hier
        rnorm(1L, log(0.5), 0.5), # log(sigma_hier)
        rnorm(8L, 0.0, 1.0) # eta[1:8]
    )
}

# ── Deterministic summary statistics from one theta ──────────────────────────
# Returns a named numeric vector of length 18 (9 age-proportions p_age +
# 9 age-prevalences q_age), or NULL if the simulation failed / produced
# non-finite or otherwise invalid output.
simulate_summary_stats <- function(theta, base_params, data) {
    pm <- build_params_from_theta(theta, base_params)
    out <- tryCatch(run_sim(pm, data), error = function(e) NULL)

    if (is.null(out) || !is.matrix(out) || nrow(out) == 0L || any(!is.finite(out))) {
        return(NULL)
    }

    y_final <- as.numeric(out[nrow(out), -1L])
    model_q <- compute_age_quantities(y_final)
    if (is.null(model_q)) {
        return(NULL)
    }

    c(model_q$p_age, model_q$q_age)
}

obs_stats <- c(obs_tot / sum(obs_tot), obs_prev)
stat_names <- c(paste0("p_age_", 1:9), paste0("q_age_", 1:9))
names(obs_stats) <- stat_names

# ── Step 1: prior predictive draws to set robust scales + calibrate eps ───────
cat("=== ABC Calibration: HCV PWID Model ===\n")
cat(sprintf("Step 1/3: %d prior predictive draws (scale + eps calibration)\n", N_PRIOR_PRED))

prior_pred_stats <- matrix(NA_real_, N_PRIOR_PRED, length(stat_names),
    dimnames = list(NULL, stat_names)
)
prior_pred_theta <- matrix(NA_real_, N_PRIOR_PRED, N_PARAMS,
    dimnames = list(NULL, param_names_log)
)
n_valid_pp <- 0L

for (ii in seq_len(N_PRIOR_PRED)) {
    th <- sample_theta_prior()
    s <- simulate_summary_stats(th, params, data)
    if (!is.null(s)) {
        n_valid_pp <- n_valid_pp + 1L
        prior_pred_stats[n_valid_pp, ] <- s
        prior_pred_theta[n_valid_pp, ] <- th
    }
    if (ii %% 500L == 0L) {
        cat(sprintf("  Prior predictive: %d / %d (valid: %d)\n", ii, N_PRIOR_PRED, n_valid_pp))
    }
}

if (n_valid_pp < 50L) {
    stop(sprintf(
        "Only %d / %d prior predictive simulations were valid. The prior is too
     diffuse relative to what run_sim()/compute_age_quantities() can handle,
     or there is a bug upstream. Fix this before calibrating eps.",
        n_valid_pp, N_PRIOR_PRED
    ))
}

prior_pred_stats <- prior_pred_stats[seq_len(n_valid_pp), , drop = FALSE]
prior_pred_theta <- prior_pred_theta[seq_len(n_valid_pp), , drop = FALSE]

# Robust per-statistic scale (MAD; falls back to sd if MAD is 0, e.g. a stat
# that is nearly constant across prior draws).
stat_scale <- apply(prior_pred_stats, 2, function(col) {
    m <- mad(col, na.rm = TRUE)
    if (!is.finite(m) || m <= 0) {
        s <- sd(col, na.rm = TRUE)
        if (!is.finite(s) || s <= 0) {
            return(1.0)
        }
        return(s)
    }
    m
})

cat("Per-statistic robust scale (MAD, fallback SD):\n")
print(round(stat_scale, 5))

# ── Weighted joint distance ───────────────────────────────────────────────────
# Standardized Euclidean distance: each of the 18 summary stats is divided by
# its prior-predictive robust scale before differencing, so no single
# statistic's raw units dominate. (Equivalent to a diagonal Mahalanobis
# distance using MAD^2 in place of the covariance diagonal; ignoring
# off-diagonal covariance between p_age stats, which are not independent
# since they sum to 1, is a deliberate simplification — see README note below
# the script if you want a full covariance version.)
abc_distance <- function(sim_stats, obs_stats, scale) {
    sqrt(sum(((sim_stats - obs_stats) / scale)^2))
}

prior_pred_dist <- apply(prior_pred_stats, 1, abc_distance, obs_stats = obs_stats, scale = stat_scale)

eps <- if (is.finite(EPS_FIXED)) {
    EPS_FIXED
} else {
    as.numeric(quantile(prior_pred_dist, EPS_QUANTILE, na.rm = TRUE))
}

cat(sprintf(
    "Distance threshold eps = %.5f (%s)\n",
    eps,
    if (is.finite(EPS_FIXED)) "fixed" else sprintf("%.1f%% quantile of prior predictive distances", 100 * EPS_QUANTILE)
))
cat(sprintf(
    "Prior predictive distance summary: min=%.4f median=%.4f max=%.4f\n",
    min(prior_pred_dist), median(prior_pred_dist), max(prior_pred_dist)
))

# ── ABC acceptance rule ───────────────────────────────────────────────────────
# Accept iff the single joint standardized distance is <= eps. No assumed
# sampling distribution anywhere in this function.
abc_screen_theta <- function(theta, base_params, data, obs_stats, scale, eps) {
    s <- simulate_summary_stats(theta, base_params, data)
    if (is.null(s)) {
        return(list(accepted = FALSE, reason = "simulation_failed", dist = NA_real_))
    }
    d <- abc_distance(s, obs_stats, scale)
    list(
        accepted = is.finite(d) && d <= eps,
        stats    = s,
        dist     = d,
        reason   = if (is.finite(d) && d <= eps) "accepted" else "distance_exceeds_eps"
    )
}

# ── Storage ───────────────────────────────────────────────────────────────────
accepted_theta <- matrix(
    NA_real_,
    nrow = N_ACCEPT_TARGET, ncol = N_PARAMS,
    dimnames = list(NULL, param_names_log)
)
accepted_ids <- integer(N_ACCEPT_TARGET)
accepted_dist <- numeric(N_ACCEPT_TARGET)
accepted_stats <- matrix(NA_real_,
    nrow = N_ACCEPT_TARGET, ncol = length(stat_names),
    dimnames = list(NULL, stat_names)
)

# Also keep ALL proposal distances (not just accepted) for diagnostics: lets
# you re-tune eps after the fact without rerunning simulations.
all_dist <- numeric(N_PROPOSAL_MAX)
all_failed <- logical(N_PROPOSAL_MAX)

# ── Step 2: rejection ABC loop ────────────────────────────────────────────────
cat(sprintf(
    "\nStep 2/3: rejection ABC | Accepted target: %d | Max proposals: %d\n",
    N_ACCEPT_TARGET, N_PROPOSAL_MAX
))

n_accepted <- 0L
n_proposed <- 0L
n_failed <- 0L

while (n_accepted < N_ACCEPT_TARGET && n_proposed < N_PROPOSAL_MAX) {
    n_proposed <- n_proposed + 1L
    theta <- sample_theta_prior()

    res <- abc_screen_theta(theta, params, data, obs_stats, stat_scale, eps)

    all_dist[n_proposed] <- res$dist
    all_failed[n_proposed] <- identical(res$reason, "simulation_failed")
    if (all_failed[n_proposed]) n_failed <- n_failed + 1L

    if (isTRUE(res$accepted)) {
        n_accepted <- n_accepted + 1L
        accepted_theta[n_accepted, ] <- theta
        accepted_ids[n_accepted] <- n_proposed
        accepted_dist[n_accepted] <- res$dist
        accepted_stats[n_accepted, ] <- res$stats
    }

    if (n_proposed %% 500L == 0L || (n_accepted > 0L && n_accepted %% 250L == 0L)) {
        cat(sprintf(
            "Proposals=%d | Accepted=%d | Failed=%d (sim errors) | AcceptRate=%.4f\n",
            n_proposed, n_accepted, n_failed,
            if (n_proposed > 0L) n_accepted / n_proposed else 0.0
        ))
    }
}

all_dist <- all_dist[seq_len(n_proposed)]
all_failed <- all_failed[seq_len(n_proposed)]

# ── Hard stop on inadequate sample (not just a warning) ──────────────────────
MIN_USABLE_ACCEPT <- 200L
if (n_accepted < MIN_USABLE_ACCEPT) {
    stop(sprintf(
        "ABC produced only %d accepted samples (target was %d) after %d proposals.
     This is below the minimum usable threshold (%d) and the script is
     stopping rather than silently proceeding to PPC/plots on an inadequate
     sample. Likely causes: eps is too tight (current eps=%.5f; prior
     predictive median distance was %.4f), the prior is a poor proposal
     distribution for this posterior, or there is a structural mismatch
     between the model and the data. Inspect 'abc/abc_diagnostics.rds' (saved
     below) and consider relaxing EPS_QUANTILE or increasing N_PROPOSAL_MAX.",
        n_accepted, N_ACCEPT_TARGET, n_proposed, MIN_USABLE_ACCEPT,
        eps, median(prior_pred_dist)
    ))
}
if (n_accepted < N_ACCEPT_TARGET) {
    warning(sprintf(
        "ABC stopped at the proposal budget with %d / %d target samples accepted
     (%d proposals used). Proceeding with the partial sample since it clears
     the minimum usable threshold (%d), but treat downstream posterior
     summaries as having a smaller effective sample size than requested.",
        n_accepted, N_ACCEPT_TARGET, n_proposed, MIN_USABLE_ACCEPT
    ))
} else {
    cat(sprintf("Reached accepted-sample target with %d proposals.\n", n_proposed))
}

accepted_theta <- accepted_theta[seq_len(n_accepted), , drop = FALSE]
accepted_ids <- accepted_ids[seq_len(n_accepted)]
accepted_dist <- accepted_dist[seq_len(n_accepted)]
accepted_stats <- accepted_stats[seq_len(n_accepted), , drop = FALSE]

cat(sprintf(
    "ABC acceptance rate: %d / %d = %.5f | sim failure rate: %d / %d = %.5f\n",
    n_accepted, n_proposed, n_accepted / n_proposed,
    n_failed, n_proposed, n_failed / n_proposed
))

# ── Build accepted sample table ───────────────────────────────────────────────
accepted_orig <- theta_to_orig(accepted_theta)
accepted_log_df <- as.data.frame(accepted_theta)
colnames(accepted_log_df) <- param_names_log
accepted_orig_df <- as.data.frame(accepted_orig)

accepted_df <- bind_cols(
    data.frame(
        proposal_id = accepted_ids,
        distance = accepted_dist,
        stringsAsFactors = FALSE
    ),
    accepted_log_df,
    accepted_orig_df,
    as.data.frame(accepted_stats, stringsAsFactors = FALSE)
)
accepted_df$eps <- eps
accepted_df$accepted_rate <- n_accepted / n_proposed
accepted_df$proposal_total <- n_proposed

# ── Save accepted samples + full diagnostics ──────────────────────────────────
saveRDS(
    list(
        accepted_theta   = accepted_theta,
        accepted_ids     = accepted_ids,
        accepted_dist    = accepted_dist,
        accepted_df      = accepted_df,
        n_accepted       = n_accepted,
        n_proposed       = n_proposed,
        n_failed         = n_failed,
        eps              = eps,
        eps_quantile     = EPS_QUANTILE,
        stat_scale       = stat_scale,
        prior_pred_dist  = prior_pred_dist,
        prior_pred_stats = prior_pred_stats,
        prior_pred_theta = prior_pred_theta,
        all_dist         = all_dist,
        all_failed       = all_failed,
        obs_pos          = obs_pos,
        obs_tot          = obs_tot,
        obs_stats        = obs_stats
    ),
    file = file.path(ABC_DIR, "abc_diagnostics.rds")
)

saveRDS(accepted_theta, file = file.path(ABC_DIR, "abc_accepted_theta.rds"))
write.csv(accepted_df, file = file.path(ABC_DIR, "abc_accepted_samples.csv"), row.names = FALSE)

# ── Step 3: posterior-predictive comparison (link-agnostic) ──────────────────
# IMPORTANT: this is NOT the same PPC as HMC_conv.R's generate_ppc_samples(),
# which draws synthetic counts from the multinomial/beta-binomial link you do
# not trust. Here, "predictive uncertainty" comes entirely from the spread of
# deterministic model outputs (p_age, q_age) across accepted theta draws —
# i.e., posterior uncertainty in the simulator's output, not from an assumed
# sampling distribution on top of it. This is the correct PPC notion for a
# deterministic simulator under likelihood-free inference.
cat(sprintf("\nStep 3/3: building posterior summaries from %d accepted draws\n", n_accepted))

ppc_n <- min(N_PPC, n_accepted)
ppc_idx <- sort(sample(n_accepted, ppc_n))
ppc_p_age <- accepted_stats[ppc_idx, 1:9, drop = FALSE]
ppc_q_age <- accepted_stats[ppc_idx, 10:18, drop = FALSE]

age_lbl <- paste0("Age ", 1:9)

summarise_stat <- function(mat, obs_vals, label) {
    data.frame(
        Age = factor(age_lbl, levels = age_lbl),
        obs = obs_vals,
        med = apply(mat, 2, median, na.rm = TRUE),
        lo50 = apply(mat, 2, quantile, 0.25, na.rm = TRUE),
        hi50 = apply(mat, 2, quantile, 0.75, na.rm = TRUE),
        lo95 = apply(mat, 2, quantile, 0.025, na.rm = TRUE),
        hi95 = apply(mat, 2, quantile, 0.975, na.rm = TRUE),
        type = label,
        stringsAsFactors = FALSE
    )
}

prop_df <- summarise_stat(ppc_p_age, obs_tot / sum(obs_tot), "Age proportion (model)")
prev_df <- summarise_stat(ppc_q_age, obs_prev, "HCV prevalence (model)")

plot_model_intervals <- function(df, y_lab, pct = FALSE) {
    p <- ggplot(df, aes(x = Age)) +
        geom_linerange(aes(ymin = lo95, ymax = hi95), linewidth = 1.0, colour = "#4292C6", alpha = 0.55) +
        geom_linerange(aes(ymin = lo50, ymax = hi50), linewidth = 3.0, colour = "#08519C", alpha = 0.80) +
        geom_point(aes(y = med), shape = 18, size = 3.5, colour = "#08306B") +
        geom_point(aes(y = obs), shape = 21, size = 3.5, fill = "#D7191C", colour = "black", stroke = 1.2) +
        labs(
            title = "ABC posterior: deterministic model output vs observed",
            subtitle = "Diamond: posterior median of simulator output | Thick bar: 50% interval | Thin bar: 95% interval | Red dot: observed",
            x = "Age group", y = y_lab
        ) +
        theme_bw(base_size = 11) +
        theme(panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 9, colour = "grey40"))
    if (pct) p <- p + scale_y_continuous(labels = function(x) sprintf("%.1f%%", 100 * x))
    p
}

p_prop_interval <- plot_model_intervals(prop_df, "Proportion of PWID population", pct = TRUE)
p_prev_interval <- plot_model_intervals(prev_df, "HCV prevalence", pct = TRUE)

p_density <- as.data.frame(theta_to_orig(accepted_theta)) %>%
    setNames(param_names_orig) %>%
    pivot_longer(everything(), names_to = "param", values_to = "value") %>%
    mutate(param = factor(param, levels = param_names_orig)) %>%
    ggplot(aes(x = value)) +
    geom_density(fill = "#2166AC", alpha = 0.35, linewidth = 0.6, colour = "#2166AC") +
    facet_wrap(~param, scales = "free", ncol = 2) +
    labs(title = "ABC posterior densities (original scale)", x = "Parameter value", y = "Density") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())

p_dist_hist <- data.frame(distance = prior_pred_dist) %>%
    ggplot(aes(x = distance)) +
    geom_histogram(bins = 50, fill = "grey60", colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = eps, colour = "#D7191C", linewidth = 1.0, linetype = "dashed") +
    labs(
        title = "Prior predictive distance distribution",
        subtitle = sprintf("Dashed line: eps = %.4f (accepts closest %.1f%% of prior draws)", eps, 100 * EPS_QUANTILE),
        x = "Standardized joint distance to observed data", y = "Count"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 9, colour = "grey40"))

print(p_density)
print(p_prop_interval)
print(p_prev_interval)
print(p_dist_hist)

ggsave(file.path(ABC_DIR, "abc_posterior_densities.png"), p_density, width = 11, height = 8, dpi = 300)
ggsave(file.path(ABC_DIR, "abc_proportion_intervals.png"), p_prop_interval, width = 8, height = 5, dpi = 300)
ggsave(file.path(ABC_DIR, "abc_prevalence_intervals.png"), p_prev_interval, width = 8, height = 5, dpi = 300)
ggsave(file.path(ABC_DIR, "abc_distance_distribution.png"), p_dist_hist, width = 8, height = 5, dpi = 300)

saveRDS(
    list(
        p_density        = p_density,
        p_prop_interval  = p_prop_interval,
        p_prev_interval  = p_prev_interval,
        p_dist_hist      = p_dist_hist,
        prop_df          = prop_df,
        prev_df          = prev_df
    ),
    file = file.path(ABC_DIR, "abc_ppc_results.rds")
)

cat("\nAll ABC results saved to abc/\n")

# =============================================================================
# NOTES ON WHAT CHANGED AND WHY (kept here rather than in a separate README so
# they travel with the script)
# =============================================================================
# 1. Single joint distance + tunable eps replaces independent per-statistic
#    interval screening. eps is calibrated from the prior predictive distance
#    distribution (Step 1) rather than hard-coded, and is fully inspectable /
#    re-tunable from the saved 'prior_pred_dist' vector without rerunning any
#    simulations.
# 2. No qbinom/qbeta/dmultinom/beta-binomial anywhere in the acceptance path.
#    Comparison is purely between deterministic simulator summary stats and
#    observed summary stats (age proportions, age prevalences) — appropriate
#    for a deterministic ODE simulator with an untrusted observation link.
# 3. Each of the 18 summary statistics is standardized by its robust prior
#    predictive scale (MAD) before differencing, so age-totals and
#    prevalences contribute comparably to the distance instead of one
#    dominating by raw units.
# 4. phi_overdisp / data$phi_overdisp is no longer referenced anywhere in this
#    script: that parameter belonged exclusively to the beta-binomial link
#    being avoided here. If you still want it for a *different*, trusted-link
#    analysis (e.g. HMC_core.R), keep setting it there, not here.
# 5. All proposal distances (accepted and rejected) are stored in
#    'all_dist'/'prior_pred_dist', not just a binary accept/reject flag, so
#    you can do post-hoc regression adjustment (Beaumont et al. 2002) or
#    re-derive a different eps later.
# 6. A hard stop (via stop()) fires if accepted samples fall below a minimum
#    usable threshold (200), instead of a warning() that lets the rest of the
#    script run on a degenerate or empty sample.
# 7. Acceptance-rate accounting separates simulation failures (n_failed) from
#    the eps-based reject decision, so a low acceptance rate from "ODE solver
#    instability for bad theta draws" can be told apart from "the prior is a
#    poor match to the data."
# 8. Simplification kept deliberately: the distance treats the 18 statistics
#    as having independent scales (diagonal weighting), not a full 18x18
#    covariance matrix. p_age sums to 1 across age groups so the statistics
#    are not actually independent; a full Mahalanobis distance using the
#    empirical covariance of prior_pred_stats would be a defensible next
#    step if diagonal weighting turns out to give poor coverage, but it adds
#    estimation noise from inverting an 18x18 covariance matrix with only
#    N_PRIOR_PRED ~ 3000 draws, so MAD-diagonal weighting was chosen as the
#    safer default here.