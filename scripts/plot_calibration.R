#!/usr/bin/env Rscript
# =============================================================================
# plot_calibration.R — compares the calibration target intervals (seroprevalence and first
# arrests by age, community PWID) with posterior predictive intervals of the model values from
# the reference MCMC posterior and single-round NPE (calibration figure).
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) {
  i <- match(name, args)
  if (is.na(i)) default else args[i + 1]
}
root <- normalizePath(arg("--root", getwd()), mustWork = TRUE)
fit_path <- normalizePath(arg("--fit", file.path(root,
  "output/calibration/r_20261027_agepop/fit.rds")), mustWork = TRUE)
mcmc_file <- normalizePath(arg("--mcmc-file", file.path(root,
  "output/calibration/mcmc_reference/mcmc_draws.csv")), mustWork = TRUE)
npe_file <- normalizePath(arg("--npe-file", file.path(root,
  "output/calibration/npe_lhs/posterior_samples_npe.csv")), mustWork = TRUE)
MCMC_LABEL <- arg("--mcmc-label", "Reference MCMC (production)")
NPE_LABEL <- arg("--npe-label", "Single-round NPE (50,000 LHS)")
out_file <- file.path(root, arg("--out", "figures/calibration_agepop.png"))
cores <- as.integer(arg("--cores", "4"))
max_draws <- as.integer(arg("--max-draws", "2000"))
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)

source(file.path(root, "scripts/calibration/targets.R"))
targets <- read_targets(file.path(root, "data/calibration_targets.json"))
old <- setwd(file.path(root, "scripts/model"))
setup_env <- new.env(parent = globalenv())
sys.source("setup.R", setup_env)
setwd(old)
Rcpp::sourceCpp(file.path(root, "scripts/model/sim.cpp"))
for (f in c("model_metrics.R", "equilibrium.R", "likelihood.R"))
  sys.source(file.path(root, "scripts/calibration", f), envir = environment())
fit <- readRDS(fit_path)
check_fit_current(fit, setup_env$params, setup_env$data)

age_names <- targets$age_groups
theta_names <- paste0("theta", 1:N_THETA)

evaluate_draw <- function(theta) {
  tryCatch({
    pm <- build_params(theta, fit$base_params)
    data <- setup_env$data
    data <- equilibrium_data(data)
    out <- run_sim(pm, data)
    eq <- check_equilibrium(out)
    cs <- calib_summary(out, pm)
    v <- c(setNames(cs$p_sero, paste0("p", seq_along(age_names))),
           setNames(cs$first_arrests, paste0("fa", seq_along(age_names))),
           community = cs$community_pwid)
    if (any(!is.finite(v))) return(NULL)
    c(v, equilibrium_pass = isTRUE(eq$pass))
  }, error = function(e) NULL)
}

load_draws <- function(files) {
  x <- do.call(rbind, lapply(files, function(f) read.csv(f, check.names = FALSE)))
  x <- x[, theta_names, drop = FALSE]
  if (nrow(x) > max_draws) {
    set.seed(20260929)
    x <- x[sample.int(nrow(x), max_draws), , drop = FALSE]
  }
  as.matrix(x)
}

simulate_or_load <- function(theta, cache_file, label, source_files) {
  cache_fresh <- file.exists(cache_file) &&
    file.info(cache_file)$mtime >= max(file.info(source_files)$mtime)
  if (isTRUE(cache_fresh)) {
    cached <- read.csv(cache_file, check.names = FALSE)
    if (identical(cached$n_requested[1], nrow(theta))) return(cached)
  }
  message("Evaluating ", nrow(theta), " ", label, " draws with ", cores, " cores...")
  if (cores > 1 && .Platform$OS.type != "windows") {
    result <- parallel::mclapply(seq_len(nrow(theta)),
      function(i) evaluate_draw(theta[i, ]), mc.cores = cores)
  } else {
    result <- lapply(seq_len(nrow(theta)), function(i) evaluate_draw(theta[i, ]))
  }
  keep <- vapply(result, function(x) is.numeric(x) && length(x) == 2L * length(age_names) + 2L, logical(1))
  if (!any(keep)) stop("No finite ", label, " predictive draws available")
  out <- as.data.frame(do.call(rbind, result[keep])); out$n_requested <- nrow(theta)
  write.csv(out, cache_file, row.names = FALSE)
  out
}

# The cache is fresh only if newer than the draws and every model input and made from the same
# number of draws (n_requested column), so a different --max-draws never reuses it.
model_inputs <- c(file.path(root, c("scripts/model/setup.R", "scripts/model/sim.cpp",
  "data/calibration_targets.json", "scripts/calibration/likelihood.R",
  "scripts/calibration/model_metrics.R", "scripts/calibration/equilibrium.R")), fit_path)
cache_for <- function(f) file.path(dirname(f), "model_prediction_2015_age_draws.csv")
npe <- simulate_or_load(load_draws(npe_file), cache_for(npe_file), "NPE posterior", c(npe_file, model_inputs))
mcmc <- simulate_or_load(load_draws(mcmc_file), cache_for(mcmc_file), "MCMC posterior", c(mcmc_file, model_inputs))

# One row per target: measure, label, column of the predictive draws, and the target interval.
# Target intervals are the 95% range of each likelihood around the observation: Beta for
# seroprevalence, negative binomial for first arrests, log-normal for community PWID.
target_table <- function() {
  y <- targets$seroprevalence; phi <- targets$beta_precision
  fa <- targets$first_arrests; size <- targets$first_arrests_nb_size
  m <- targets$community_pwid; sdl <- targets$community_pwid_log_sd
  data.frame(
    measure = c(rep("seroprevalence", 6), rep("first_arrests", 6), "community_pwid"),
    group = c(age_names, age_names, "All ages (D + F)"),
    column = c(paste0("p", 1:6), paste0("fa", 1:6), "community"),
    target = c(y, fa, m),
    target_lower = c(qbeta(0.025, y * phi, (1 - y) * phi), qnbinom(0.025, size = size, mu = fa),
                     m * exp(-1.959964 * sdl)),
    target_upper = c(qbeta(0.975, y * phi, (1 - y) * phi), qnbinom(0.975, size = size, mu = fa),
                     m * exp(1.959964 * sdl)),
    stringsAsFactors = FALSE)
}

summarise_draws <- function(x, source) {
  x <- x[as.logical(x$equilibrium_pass), , drop = FALSE]
  if (!nrow(x)) stop("No equilibrium-pass ", source, " draws")
  tt <- target_table()
  data.frame(tt[c("measure", "group")], source = source,
             estimate = vapply(tt$column, function(k) median(x[[k]]), numeric(1)),
             lower = vapply(tt$column, function(k) unname(quantile(x[[k]], 0.025)), numeric(1)),
             upper = vapply(tt$column, function(k) unname(quantile(x[[k]], 0.975)), numeric(1)),
             n_draws = nrow(x), row.names = NULL, stringsAsFactors = FALSE)
}

target_intervals <- function() {
  tt <- target_table()
  data.frame(tt[c("measure", "group")], source = "Target", estimate = tt$target,
             lower = tt$target_lower, upper = tt$target_upper, n_draws = NA_integer_)
}

source_levels <- c("Target", MCMC_LABEL, NPE_LABEL)
source_shapes <- setNames(c(17, 15, 16), source_levels)
source_cols <- setNames(c("black", "#0072B2", "#D55E00"), source_levels)

panel <- function(d, which, y_label, title, formatter, groups) {
  d <- d[d$measure == which, ]
  d$group <- factor(d$group, levels = groups)
  d$source <- factor(d$source, levels = source_levels)
  ggplot(d, aes(x = group, y = estimate, colour = source, shape = source)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.12, linewidth = 0.6,
                  position = position_dodge(width = 0.45)) +
    geom_point(size = 2.2, position = position_dodge(width = 0.45)) +
    scale_colour_manual(values = source_cols, name = NULL) +
    scale_shape_manual(values = source_shapes, name = NULL) +
    scale_y_continuous(labels = formatter, limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = y_label, title = title) +
    theme_classic(base_size = 11) + theme(legend.position = "bottom")
}

d_all <- rbind(target_intervals(), summarise_draws(mcmc, MCMC_LABEL), summarise_draws(npe, NPE_LABEL))
figure <- (panel(d_all, "seroprevalence", "HCV seroprevalence in prison", "Seroprevalence, 2015",
                 function(x) paste0(round(100 * x), "%"), age_names) /
           panel(d_all, "first_arrests", "People per year", "First arrests (CNB new drug abusers), 2015",
                 scales::comma, age_names) |
           panel(d_all, "community_pwid", "People", "Community PWID (D + F)", scales::comma,
                 "All ages (D + F)")) +
  plot_layout(widths = c(3, 1), guides = "collect") & theme(legend.position = "bottom")
ggsave(out_file, figure, width = 10, height = 7, dpi = 300, bg = "white")

tt <- target_table(); sm <- summarise_draws(mcmc, "MCMC"); sn <- summarise_draws(npe, "NPE")
check <- data.frame(measure = tt$measure, age_group = tt$group, target = tt$target,
  target_lower = tt$target_lower, target_upper = tt$target_upper,
  NPE_lower = sn$lower, NPE_upper = sn$upper, MCMC_lower = sm$lower, MCMC_upper = sm$upper)
check$NPE_overlap <- check$NPE_upper >= check$target_lower & check$NPE_lower <= check$target_upper
check$MCMC_overlap <- check$MCMC_upper >= check$target_lower & check$MCMC_lower <= check$target_upper
write.csv(check, file.path(dirname(out_file), "calibration_interval_overlap.csv"), row.names = FALSE)
message("Wrote ", out_file, "; NPE draws used: ",
        sum(as.logical(npe$equilibrium_pass)), "; MCMC draws used: ",
        sum(as.logical(mcmc$equilibrium_pass)), "; MCMC overlaps: ", sum(check$MCMC_overlap), "/", nrow(check))
