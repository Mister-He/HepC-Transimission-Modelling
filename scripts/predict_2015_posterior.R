#!/usr/bin/env Rscript
# =============================================================================
# predict_2015_posterior.R — posterior predictions of 2015 population and active HCV by
# stratum (D, J, F, X), from the reference MCMC draws unless --posterior-files is given.
# =============================================================================
# Active HCV cases are a+c+t (viremic states); recovered r is not active.

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) {
  i <- match(name, args); if (is.na(i)) default else args[i + 1]
}
ROOT <- normalizePath(arg("--root", getwd()), mustWork = TRUE)
OUT <- file.path(ROOT, arg("--out-dir", "output/calibration/mcmc_reference"))
FIT <- normalizePath(arg("--fit", file.path(ROOT,
  "output/calibration/r_20261013_agepop/fit.rds")), mustWork = TRUE)
CORES <- as.integer(arg("--cores", "4"))
MAX_DRAWS <- as.integer(arg("--max-draws", "4000"))   # 0 = use all draws
TARGET_YEAR <- as.integer(arg("--year", "2015"))
if (!TARGET_YEAR %in% c(2015L, 2016L, 2017L))
  stop("--year must be 2015, 2016, or 2017")
MODEL_YEAR <- 150L + TARGET_YEAR - 2015L
FILE_TAG <- as.character(TARGET_YEAR)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

source(file.path(ROOT, "scripts/calibration/targets.R"))
targets <- read_targets(file.path(ROOT, "data/calibration_targets.json"))
old <- setwd(file.path(ROOT, "scripts/model")); setup_env <- new.env(parent = globalenv())
sys.source("setup.R", setup_env); setwd(old)
Rcpp::sourceCpp(file.path(ROOT, "scripts/model/sim.cpp"))
for (f in c("model_metrics.R", "equilibrium.R", "likelihood.R"))
  sys.source(file.path(ROOT, "scripts/calibration", f), envir = environment())
fit <- readRDS(FIT)

sample_files <- strsplit(arg("--posterior-files", "output/calibration/mcmc_reference/mcmc_draws.csv"), ",")[[1]]
sample_files <- ifelse(grepl("^/", sample_files), sample_files, file.path(ROOT, sample_files))
names(sample_files) <- basename(dirname(sample_files))
if (length(sample_files) == 0) stop("No posterior samples found")
draws <- do.call(rbind, lapply(names(sample_files), function(source) {
  x <- read.csv(sample_files[[source]]); x$posterior_source <- source; x
}))
theta_names <- paste0("theta", 1:12)
if (MAX_DRAWS > 0 && nrow(draws) > MAX_DRAWS) {
  set.seed(20261011); draws <- draws[sample.int(nrow(draws), MAX_DRAWS), , drop = FALSE]
}

idx0 <- function(s, k, h, i) s * 120 + (k - 1) * 30 + h * 6 + i + 2L
stratum_idx <- function(s, states) as.vector(sapply(1:4, function(k)
  sapply(states, function(h) sapply(0:5, function(i) idx0(s, k, h, i) - 1L))))
population_idx <- lapply(0:3, stratum_idx, states = 0:4)
active_idx <- lapply(0:3, stratum_idx, states = 1:3)

one_prediction <- function(theta) {
  tryCatch({
    pm <- build_params(theta, fit$base_params)
    data <- setup_env$data; data$t_start <- 0
    data$t_end <- MODEL_YEAR
    out <- run_sim(pm, data)
    eq_out <- out[out[, 1] <= targets$equilibrium_horizon_years + 1e-8, , drop = FALSE]
    eq <- check_equilibrium(eq_out)
    state <- out[which.min(abs(out[, 1] - MODEL_YEAR)), -1]
    total <- vapply(population_idx, function(ii) sum(state[ii]), numeric(1))
    active <- vapply(active_idx, function(ii) sum(state[ii]), numeric(1))
    values <- c(total[1:4], total[1] + total[3], sum(total),
                active[1:4], active[1] + active[3], sum(active))
    if (any(!is.finite(values))) return(NULL)
    names(values) <- c(paste0("population_", c("D", "J", "F", "X", "current_PWID", "ever_PWID")),
                       paste0("active_hcv_", c("D", "J", "F", "X", "current_PWID", "ever_PWID")))
    list(values = values, equilibrium_pass = isTRUE(eq$pass))
  }, error = function(e) NULL)
}

theta <- as.matrix(draws[, theta_names, drop = FALSE])
cat("Evaluating", nrow(theta), "posterior draws with", CORES, "cores...\n")
if (CORES > 1 && .Platform$OS.type != "windows") {
  res <- parallel::mclapply(seq_len(nrow(theta)), function(i) one_prediction(theta[i, ]),
    mc.cores = CORES)
} else {
  res <- lapply(seq_len(nrow(theta)), function(i) one_prediction(theta[i, ]))
}
keep <- vapply(res, function(x) is.list(x) && is.numeric(x$values) &&
  length(x$values) == 12, logical(1))
if (!any(keep)) stop("No finite posterior draws available")
valid_res <- res[keep]
all_mat <- do.call(rbind, lapply(valid_res, `[[`, "values"))
equilibrium_flag <- vapply(valid_res, `[[`, logical(1), "equilibrium_pass")
n_equilibrium <- sum(equilibrium_flag)
if (n_equilibrium == 0) stop("No finite posterior draws passed the equilibrium gate")
mat <- all_mat[equilibrium_flag, , drop = FALSE]
write.csv(as.data.frame(mat), file.path(OUT, paste0("model_prediction_", FILE_TAG, "_draws.csv")), row.names = FALSE)

summary_one <- function(x) c(median = median(x), lower_95 = unname(quantile(x, .025)),
                              upper_95 = unname(quantile(x, .975)))
summary <- do.call(rbind, lapply(seq_len(ncol(mat)), function(j) summary_one(mat[, j])))
summary <- data.frame(measure = colnames(mat), summary, row.names = NULL, check.names = FALSE)
summary$n_draws_used <- nrow(mat)
summary$n_finite_draws <- nrow(all_mat)
summary$n_equilibrium_pass <- n_equilibrium
summary$posterior_sources <- paste(names(sample_files), collapse = "+")

data_best <- setup_env$data; data_best$t_start <- 0
data_best$t_end <- MODEL_YEAR
best_out <- run_sim(build_params(fit$best$theta, fit$base_params), data_best)
best_state <- best_out[which.min(abs(best_out[, 1] - MODEL_YEAR)), -1]
best_total <- vapply(population_idx, function(ii) sum(best_state[ii]), numeric(1))
best_active <- vapply(active_idx, function(ii) sum(best_state[ii]), numeric(1))
point <- c(best_total[1:4], best_total[1] + best_total[3], sum(best_total),
           best_active[1:4], best_active[1] + best_active[3], sum(best_active))
names(point) <- c(paste0("population_", c("D", "J", "F", "X", "current_PWID", "ever_PWID")),
                  paste0("active_hcv_", c("D", "J", "F", "X", "current_PWID", "ever_PWID")))
summary$point_estimate <- unname(point[match(summary$measure, names(point))])
summary$definition <- ifelse(grepl("active_hcv", summary$measure),
  "a+c+t; viremic/active HCV states", "u+a+c+t+r; all model states")
write.csv(summary, file.path(OUT, paste0("model_prediction_", FILE_TAG, "_strata.csv")), row.names = FALSE)

pretty <- data.frame(
  group = c("D", "J", "F", "X", "Current PWID (D+F)", "Ever-PWID (D+J+F+X)"),
  total_population = sprintf("%.0f (%.0f–%.0f)", summary$median[1:6],
    summary$lower_95[1:6], summary$upper_95[1:6]),
  active_hcv_cases = sprintf("%.0f (%.0f–%.0f)", summary$median[7:12],
    summary$lower_95[7:12], summary$upper_95[7:12]), stringsAsFactors = FALSE)
write.csv(pretty, file.path(OUT, paste0("model_prediction_", FILE_TAG, "_strata_pretty.csv")), row.names = FALSE)
writeLines(c(
  paste0("# Model prediction at ", TARGET_YEAR), "",
  paste0("Posterior medians and equal-tailed 95% Bayesian credible intervals from ", nrow(draws), " posterior draws (", paste(sample_files, collapse = "; "), ")."),
  "Active HCV cases are a+c+t (viremic states); recovered r is antibody-reactive but not active.", "",
  "| Group | Total population | Active HCV cases |", "|---|---:|---:|",
  apply(pretty, 1, function(z) paste0("| ", z[1], " | ", z[2], " | ", z[3], " |")), "",
  paste0("Finite posterior draws: ", nrow(all_mat), " / ", nrow(draws), "; equilibrium-pass draws used: ", n_equilibrium, "."),
  paste0("Posterior sources: ", paste(names(sample_files), collapse = "+"), ".")
), file.path(OUT, paste0("model_prediction_", FILE_TAG, "_strata.md")))

png(file.path(OUT, paste0("model_prediction_", FILE_TAG, "_strata.png")), width = 2000, height = 1300, res = 180)
par(mar = c(13, 5, 3, 1), las = 2)
groups <- pretty$group; xpos <- seq_along(groups)
pop <- summary[1:6, ]; active <- summary[7:12, ]
ylim <- range(c(pop$lower_95, pop$upper_95, active$lower_95, active$upper_95))
plot(xpos, pop$median, pch = 19, col = "#0072B2", ylim = ylim, xaxt = "n",
  xlab = "", ylab = "People", main = paste0(TARGET_YEAR, " posterior prediction by stratum"))
segments(xpos - .12, pop$lower_95, xpos - .12, pop$upper_95, col = "#0072B2", lwd = 2)
points(xpos - .12, pop$median, pch = 19, col = "#0072B2")
segments(xpos + .12, active$lower_95, xpos + .12, active$upper_95, col = "#D55E00", lwd = 2)
points(xpos + .12, active$median, pch = 17, col = "#D55E00")
axis(1, at = xpos, labels = groups)
legend("topright", c("Total population", "Active HCV cases (a+c+t)"),
  col = c("#0072B2", "#D55E00"), pch = c(19, 17), bty = "n")
dev.off()
cat("Wrote", TARGET_YEAR, "posterior predictions using", nrow(mat), "draws to", OUT, "\n")
