#!/usr/bin/env Rscript
# =============================================================================
# compare_npe_mcmc.R — independent cross-validation: compares the single-round NPE
# posterior (50,000 LHS prior simulations) with the reference MCMC posterior for
# parameters, joint structure and policy outcomes; writes output/calibration/npe_vs_mcmc/.
#
# Tolerances (fixed before comparison): |standardised mean difference| <= 0.1,
# |standardised 95% endpoint difference| <= 0.2, 95% width ratio in [0.9, 1.1].
# NPE draws outside the MCMC support (bounds, equilibrium) are dropped first. The NPE policy
# run (1,000 of the remaining states) is generated if policy/npe_lhs/ is missing or
# older than the NPE posterior or the production scenario draws, and is compared with the
# production scenario draws in output/scenarios_2026/. Nothing here feeds the production path.
# Usage: Rscript scripts/compare_npe_mcmc.R [--root .] [--policy-cores 14]
# =============================================================================
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) { i <- match(name, args); if (is.na(i)) default else args[i + 1] }
ROOT <- normalizePath(arg("--root", getwd()), mustWork = TRUE); setwd(ROOT)
POLICY_CORES <- as.integer(arg("--policy-cores", "14"))
CAL <- "output/calibration"; OUT <- file.path(CAL, "npe_vs_mcmc")
dir.create(file.path(OUT, "figures"), recursive = TRUE, showWarnings = FALSE)
TH <- paste0("theta", 1:7); probs <- c(.025, .05, .25, .5, .75, .95, .975)
TOL <- list(smd = 0.1, endpoint = 0.2, width = c(0.9, 1.1))

# ---------------- posteriors ----------------
ref <- as.matrix(read.csv(file.path(CAL, "mcmc_reference/mcmc_draws.csv"))[, TH])
ref_diag <- read.csv(file.path(CAL, "mcmc_reference/mcmc_diagnostics.csv"))
ref_diag <- ref_diag[ref_diag$subset == "all_retained", ]
# NPE draws get the MCMC's support: 4,000 random draws, keeping those within the bounds that pass
# the 2015 equilibrium check (the MCMC gives the others -Inf).
npe_all <- as.matrix(read.csv(file.path(CAL, "npe_lhs/posterior_samples_npe.csv"))[, TH])
old <- setwd("scripts/model"); setup_env <- new.env(parent = globalenv()); sys.source("setup.R", setup_env); setwd(old)
Rcpp::sourceCpp("scripts/model/sim.cpp")
targets <- jsonlite::fromJSON("data/calibration_targets.json")
for (f in c("model_metrics.R", "equilibrium.R", "likelihood.R")) sys.source(file.path("scripts/calibration", f), envir = environment())
eq_data <- setup_env$data; eq_data$t_start <- 0; eq_data$t_end <- targets$equilibrium_horizon_years
eq_data$save_every <- round(1 / eq_data$dt)
set.seed(20261015); npe_try <- npe_all[sample.int(nrow(npe_all), min(4000, nrow(npe_all))), ]
npe_ok <- unlist(parallel::mclapply(seq_len(nrow(npe_try)), function(i) within_bounds(npe_try[i, ]) &&
  isTRUE(tryCatch(check_equilibrium(run_sim(build_params(npe_try[i, ], setup_env$params), eq_data))$pass,
                  error = function(e) FALSE)), mc.cores = POLICY_CORES))
writeLines(c(paste("npe_draws_checked", length(npe_ok)), paste("npe_draws_equilibrium_pass", sum(npe_ok))),
           file.path(OUT, "npe_support.txt"))
samples <- list(npe_lhs = npe_try[npe_ok, , drop = FALSE])
label <- function(k) "Single-round NPE, 50,000 LHS prior simulations"

stats_of <- function(m) { s <- t(apply(m, 2, function(v) c(mean(v), sd(v), quantile(v, probs, names = FALSE))))
  colnames(s) <- c("mean", "sd", paste0("q", c("025", "05", "25", "50", "75", "95", "975"))); s }
rs <- stats_of(ref)
pc <- do.call(rbind, lapply(names(samples), function(k) {
  s <- stats_of(samples[[k]])
  data.frame(run = k, label = label(k), parameter = TH, n_draws = nrow(samples[[k]]), s,
             ref_mean = rs[, "mean"], ref_sd = rs[, "sd"], ref_q025 = rs[, "q025"], ref_q975 = rs[, "q975"],
             smd = (s[, "mean"] - rs[, "mean"]) / rs[, "sd"],
             lo95_diff = (s[, "q025"] - rs[, "q025"]) / rs[, "sd"],
             hi95_diff = (s[, "q975"] - rs[, "q975"]) / rs[, "sd"],
             sd_ratio = s[, "sd"] / rs[, "sd"],
             width95_ratio = (s[, "q975"] - s[, "q025"]) / (rs[, "q975"] - rs[, "q025"]),
             width90_ratio = (s[, "q95"] - s[, "q05"]) / (rs[, "q95"] - rs[, "q05"]),
             width50_ratio = (s[, "q75"] - s[, "q25"]) / (rs[, "q75"] - rs[, "q25"]),
             finite_draw_se_sd_units = 1 / sqrt(nrow(samples[[k]])),
             ref_mcse_mean_sd_units = ref_diag$mcse_mean / rs[, "sd"], row.names = NULL)
}))
pc$pass_mean <- abs(pc$smd) <= TOL$smd
pc$pass_endpoints <- pmax(abs(pc$lo95_diff), abs(pc$hi95_diff)) <= TOL$endpoint
pc$pass_width <- pc$width95_ratio >= TOL$width[1] & pc$width95_ratio <= TOL$width[2]
pc$pass_all <- pc$pass_mean & pc$pass_endpoints & pc$pass_width
write.csv(pc, file.path(OUT, "parameter_comparison.csv"), row.names = FALSE)

rc <- cor(ref)
joint <- do.call(rbind, lapply(names(samples), function(k) { dc <- cor(samples[[k]]) - rc
  data.frame(run = k, max_abs_corr_diff = max(abs(dc[upper.tri(dc)])),
             frobenius_corr_diff = sqrt(sum(dc[upper.tri(dc)]^2))) }))
summary_tab <- do.call(rbind, lapply(names(samples), function(k) { x <- pc[pc$run == k, ]
  data.frame(run = k, label = label(k), max_abs_smd = max(abs(x$smd)),
             max_abs_endpoint = max(pmax(abs(x$lo95_diff), abs(x$hi95_diff))),
             mean_width95_ratio = mean(x$width95_ratio), min_width95_ratio = min(x$width95_ratio),
             max_width95_ratio = max(x$width95_ratio), share_pass_mean = mean(x$pass_mean),
             share_pass_endpoints = mean(x$pass_endpoints), share_pass_width = mean(x$pass_width),
             share_pass_all = mean(x$pass_all)) }))
summary_tab <- merge(summary_tab, joint, by = "run", sort = FALSE)

# ---------------- costs ----------------
read_json <- function(f) if (file.exists(f)) jsonlite::fromJSON(f) else list()
cost <- function(k) {
  t <- read.csv(file.path(CAL, "npe_lhs/r_timings.csv")); p <- read_json(file.path(CAL, "npe_lhs/npe_timing.json"))
  ss <- read_json(file.path(CAL, "npe_lhs/simulation_summary.json"))
  data.frame(run = k, forward_runs = ss$n_attempted, valid_simulations = ss$n_valid,
             simulation_sec = sum(t$seconds[t$stage == "prior_simulation"]), simulation_workers = t$workers[1],
             training_sec = p$train_sec %||% NA)
}
summary_tab <- merge(summary_tab, do.call(rbind, lapply(names(samples), cost)), by = "run", sort = FALSE)
write.csv(summary_tab, file.path(OUT, "accuracy_summary.csv"), row.names = FALSE)

# ---------------- policy outcomes ----------------
pol_in <- file.path(OUT, "policy", "npe_lhs")
pol_draws <- file.path(pol_in, "scenario_2030_draws.csv")
# Rerun when the cache is missing or older than the NPE posterior or the production scenarios.
deps <- file.path(ROOT, c("output/calibration/npe_lhs/posterior_samples_npe.csv",
                          "output/scenarios_2026/scenario_2030_draws.csv"))
if (!file.exists(pol_draws) || any(file.mtime(deps) > file.mtime(pol_draws))) {
  dir.create(pol_in, recursive = TRUE, showWarnings = FALSE)
  set.seed(20261014); f_in <- file.path(pol_in, "policy_input_states.csv")
  write.csv(samples$npe_lhs[sample.int(nrow(samples$npe_lhs), min(1000, nrow(samples$npe_lhs))), ], f_in, row.names = FALSE)
  st <- system2("Rscript", c("scripts/scenarios_2026.R", "--root", ".", "--posterior-files", f_in,
                             "--draws-per-seed", 1000, "--cores", POLICY_CORES, "--out-dir", pol_in,
                             "--sample-seed", 20261010, "--min-valid-share", 0),
                stdout = file.path(pol_in, "run.log"), stderr = file.path(pol_in, "run.log"))
  if (st != 0) stop("NPE-LHS policy run failed; see ", file.path(pol_in, "run.log"))
  # Keep only what this comparison reads; the production scenario set is output/scenarios_2026/.
  unlink(file.path(pol_in, setdiff(list.files(pol_in), c("scenario_2030_draws.csv", "draw_validity.txt",
                                                         "policy_input_states.csv", "calibration_source.txt", "run.log"))))
}
key <- data.frame(scenario = c("J_only", "J_only", "S5000_F25", "S5000_F25", "S8000_F25", "S8000_F75", "S8000_F25"),
                  outcome = c("incidence_reduction", "mortality_reduction", "mortality_reduction", "joint_pass_probability",
                              "mortality_reduction", "mortality_reduction", "total_discounted_cost_2026_2030_sgd"))
pol_dirs <- c(npe_lhs = pol_in)
read_pol <- function(d) list(draws = read.csv(file.path(d, "scenario_2030_draws.csv")),
                             valid = read.delim(file.path(d, "draw_validity.txt"), sep = " ", header = FALSE)$V2)
summ_out <- function(dr, sc, oc, idx = NULL) {
  x <- dr[dr$scenario_id == sc, ]; if (!is.null(idx)) x <- x[idx, ]
  if (oc == "joint_pass_probability") return(c(median = mean(as.logical(x$who_joint_pass)), lo = NA, hi = NA))
  v <- x[[oc]]; c(median = median(v), lo = unname(quantile(v, .025)), hi = unname(quantile(v, .975)))
}
boot_se <- function(dr, sc, oc, B = 200) { n <- sum(dr$scenario_id == sc); set.seed(1)
  apply(replicate(B, summ_out(dr, sc, oc, sample.int(n, n, TRUE))), 1, sd, na.rm = TRUE) }
R0 <- read_pol("output/scenarios_2026")
pol <- do.call(rbind, lapply(names(pol_dirs), function(k) { P <- read_pol(pol_dirs[[k]])
  do.call(rbind, lapply(seq_len(nrow(key)), function(i) {
    sc <- key$scenario[i]; oc <- key$outcome[i]; a <- summ_out(P$draws, sc, oc); b <- summ_out(R0$draws, sc, oc)
    is_p <- oc == "joint_pass_probability"
    err <- if (is_p) 100 * abs(a[["median"]] - b[["median"]]) else abs(a[["median"]] - b[["median"]]) / abs(b[["median"]])
    data.frame(run = k, label = label(k), scenario = sc, outcome = oc, n_valid = P$valid[2], n_requested = P$valid[1],
               median = a[["median"]], lo = a[["lo"]], hi = a[["hi"]], ref_median = b[["median"]], ref_lo = b[["lo"]],
               ref_hi = b[["hi"]], error = err, error_unit = if (is_p) "percentage points" else "relative",
               mc_se_median = boot_se(P$draws, sc, oc)[["median"]], ref_mc_se_median = boot_se(R0$draws, sc, oc)[["median"]],
               lo_rel_error = if (is_p) NA else abs(a[["lo"]] - b[["lo"]]) / abs(b[["lo"]]),
               hi_rel_error = if (is_p) NA else abs(a[["hi"]] - b[["hi"]]) / abs(b[["hi"]]),
               within_1 = err <= (if (is_p) 1 else 0.01), within_5 = err <= (if (is_p) 5 else 0.05),
               within_10 = if (is_p) NA else err <= 0.10, row.names = NULL)
  })) }))
write.csv(pol, file.path(OUT, "policy_comparison.csv"), row.names = FALSE)

# ---------------- figures ----------------
long <- function(m, lab) data.frame(value = as.vector(m), parameter = factor(rep(TH, each = nrow(m)), TH), source = lab)
set.seed(1); ref_long <- long(ref[sample.int(nrow(ref), min(nrow(ref), 20000)), ], "Reference MCMC")
dd <- long(samples$npe_lhs, "NPE, 50,000 LHS")
p <- ggplot() + geom_density(data = ref_long, aes(value), fill = "grey80", colour = NA) +
  geom_density(data = dd, aes(value, colour = source), linewidth = .5) +
  facet_wrap(~parameter, scales = "free", ncol = 4) +
  labs(title = "Single-round NPE versus reference MCMC posterior", x = "theta (log scale)", y = "Density", colour = NULL) +
  theme_classic(base_size = 9) + theme(legend.position = "bottom")
ggsave(file.path(OUT, "figures", "marginals.png"), p, width = 10, height = 7, dpi = 200)
acc <- pc; acc$parameter <- factor(acc$parameter, rev(TH))
acc$source <- "NPE, 50,000 LHS"
p1 <- ggplot(acc, aes(parameter, smd, colour = source)) + annotate("rect", xmin = -Inf, xmax = Inf, ymin = -.1, ymax = .1, alpha = .15) +
  geom_point(position = position_dodge(width = .6)) + coord_flip() +
  labs(y = "Standardised mean difference vs reference", x = NULL, colour = NULL) + theme_classic(base_size = 9)
p2 <- ggplot(acc, aes(parameter, width95_ratio, colour = source)) + annotate("rect", xmin = -Inf, xmax = Inf, ymin = .9, ymax = 1.1, alpha = .15) +
  geom_point(position = position_dodge(width = .6)) + coord_flip() +
  labs(y = "95% interval width ratio vs reference", x = NULL, colour = NULL) + theme_classic(base_size = 9)
ggsave(file.path(OUT, "figures", "accuracy.png"), (p1 | p2) + plot_layout(guides = "collect") & theme(legend.position = "bottom"),
       width = 11, height = 6, dpi = 200)
pairs_of <- list(c("theta1", "theta7"), c("theta6", "theta7"), c("theta2", "theta3"), c("theta4", "theta5"))
pp <- lapply(pairs_of, function(pr) {
  set.seed(2); mk <- function(m, lab) data.frame(x = m[, pr[1]], y = m[, pr[2]], source = lab)
  d2 <- rbind(mk(ref[sample.int(nrow(ref), 5000), ], "Reference MCMC"),
              mk(samples$npe_lhs[sample.int(nrow(samples$npe_lhs), min(5000, nrow(samples$npe_lhs))), ], "NPE, 50,000 LHS"))
  ggplot(d2, aes(x, y)) + geom_density_2d(aes(colour = source), bins = 6) + labs(x = pr[1], y = pr[2], colour = NULL) +
    theme_classic(base_size = 9) })
ggsave(file.path(OUT, "figures", "pairs.png"), wrap_plots(pp, ncol = 2, guides = "collect") & theme(legend.position = "bottom"),
       width = 9, height = 8, dpi = 200)
print(summary_tab[, c("run", "max_abs_smd", "max_abs_endpoint", "mean_width95_ratio", "share_pass_all", "forward_runs")])
print(pol[, c("run", "scenario", "outcome", "median", "ref_median", "error", "within_1")])
cat("comparison written to", OUT, "\n")
