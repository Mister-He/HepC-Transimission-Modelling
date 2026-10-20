#!/usr/bin/env Rscript
# =============================================================================
# plot_calibration.R — compares the 2015 target intervals with posterior predictive
# intervals from the reference MCMC posterior and single-round NPE (calibration figure).
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
  "output/calibration/r_20261020_agepop/fit.rds")), mustWork = TRUE)
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

age_names <- targets$age_groups
theta_names <- paste0("theta", 1:12)

evaluate_draw <- function(theta) {
  tryCatch({
    pm <- build_params(theta, fit$base_params)
    data <- setup_env$data
    data$t_start <- 0
    data$t_end <- targets$equilibrium_horizon_years
    out <- run_sim(pm, data)
    eq <- check_equilibrium(out)
    s <- J_summary_final(out)
    v <- c(setNames(s$N_hat, paste0("N", seq_along(age_names))),
           setNames(s$p_sero, paste0("p", seq_along(age_names))))
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
  if (isTRUE(cache_fresh)) return(read.csv(cache_file, check.names = FALSE))
  message("Evaluating ", nrow(theta), " ", label, " draws with ", cores, " cores...")
  if (cores > 1 && .Platform$OS.type != "windows") {
    result <- parallel::mclapply(seq_len(nrow(theta)),
      function(i) evaluate_draw(theta[i, ]), mc.cores = cores)
  } else {
    result <- lapply(seq_len(nrow(theta)), function(i) evaluate_draw(theta[i, ]))
  }
  keep <- vapply(result, function(x) is.numeric(x) && length(x) == 13, logical(1))
  if (!any(keep)) stop("No finite ", label, " predictive draws available")
  out <- as.data.frame(do.call(rbind, result[keep]))
  write.csv(out, cache_file, row.names = FALSE)
  out
}

npe <- simulate_or_load(load_draws(npe_file),
  file.path(dirname(npe_file), "model_prediction_2015_age_draws.csv"), "NPE posterior", npe_file)
mcmc <- simulate_or_load(load_draws(mcmc_file),
  file.path(dirname(mcmc_file), "model_prediction_2015_age_draws.csv"), "MCMC posterior", mcmc_file)

summarise_draws <- function(x, source, prefix) {
  x <- x[as.logical(x$equilibrium_pass), , drop = FALSE]
  if (!nrow(x)) stop("No equilibrium-pass ", source, " draws")
  data.frame(
    age_group = age_names,
    source = source,
    estimate = vapply(seq_along(age_names), function(i)
      median(x[[paste0(prefix, i)]]), numeric(1)),
    lower = vapply(seq_along(age_names), function(i)
      unname(quantile(x[[paste0(prefix, i)]], 0.025)), numeric(1)),
    upper = vapply(seq_along(age_names), function(i)
      unname(quantile(x[[paste0(prefix, i)]], 0.975)), numeric(1)),
    n_draws = nrow(x), stringsAsFactors = FALSE
  )
}

target_intervals <- function(measure) {
  if (measure == "p") {
    lower <- qbeta(0.025, targets$x_prev + 0.5,
                   targets$population - targets$x_prev + 0.5)
    upper <- qbeta(0.975, targets$x_prev + 0.5,
                   targets$population - targets$x_prev + 0.5)
    estimate <- targets$p_prev_binom
  } else {
    lower <- exp(log(targets$population) - qnorm(0.975) * targets$population_sigma_log)
    upper <- exp(log(targets$population) + qnorm(0.975) * targets$population_sigma_log)
    estimate <- targets$population
  }
  data.frame(age_group = age_names, source = "Target", estimate = estimate,
             lower = lower, upper = upper, n_draws = NA_integer_)
}

age_cols <- c("#FDE725", "#B4DE2C", "#35B779", "#31688E", "#443983", "#440154")
source_levels <- c("Target", MCMC_LABEL, NPE_LABEL)
source_shapes <- setNames(c(17, 15, 16), source_levels)

plot_calibration <- function(npe, mcmc, measure, y_label, title, formatter) {
  d <- rbind(target_intervals(measure),
             summarise_draws(mcmc, MCMC_LABEL, measure),
             summarise_draws(npe, NPE_LABEL, measure))
  d$age_group <- factor(d$age_group, levels = age_names)
  d$source <- factor(d$source, levels = source_levels)
  ggplot(d, aes(x = age_group, y = estimate, colour = age_group, shape = source)) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  width = 0.12, linewidth = 0.75,
                  position = position_dodge(width = 0.18)) +
    geom_point(size = 2.8, position = position_dodge(width = 0.18)) +
    scale_colour_manual(values = setNames(age_cols, age_names), name = "Age Group") +
    scale_shape_manual(values = source_shapes, name = "Source") +
    scale_y_continuous(labels = formatter, limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(x = "Age group", y = y_label, title = title) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, colour = "grey60", linewidth = 0.4),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      axis.ticks = element_line(colour = "grey60"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom", legend.box = "horizontal",
      legend.box.just = "center"
    ) +
    guides(
      shape = guide_legend(order = 1, override.aes = list(colour = "grey25")),
      colour = guide_legend(order = 2, nrow = 1)
    )
}

p_sero <- plot_calibration(npe, mcmc, "p", "HCV seroprevalence",
  "Age-specific HCV seroprevalence", function(x) paste0(round(100 * x), "%"))
p_population <- plot_calibration(npe, mcmc, "N", "Prison population",
  "Age-specific prison population", scales::label_comma(accuracy = 1))
figure <- (p_sero + p_population) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal",
        legend.box.just = "center")
ggsave(out_file, figure, width = 15, height = 8, dpi = 300, bg = "white")

overlap <- function(x, lo, hi) x$upper >= lo & x$lower <= hi
target_p <- target_intervals("p")
target_N <- target_intervals("N")
check <- rbind(
  data.frame(measure = "seroprevalence", age_group = age_names,
    target_lower = target_p$lower, target_upper = target_p$upper,
    NPE_lower = summarise_draws(npe, "NPE posterior", "p")$lower,
    NPE_upper = summarise_draws(npe, "NPE posterior", "p")$upper,
    MCMC_lower = summarise_draws(mcmc, "MCMC posterior", "p")$lower,
    MCMC_upper = summarise_draws(mcmc, "MCMC posterior", "p")$upper),
  data.frame(measure = "population", age_group = age_names,
    target_lower = target_N$lower, target_upper = target_N$upper,
    NPE_lower = summarise_draws(npe, "NPE posterior", "N")$lower,
    NPE_upper = summarise_draws(npe, "NPE posterior", "N")$upper,
    MCMC_lower = summarise_draws(mcmc, "MCMC posterior", "N")$lower,
    MCMC_upper = summarise_draws(mcmc, "MCMC posterior", "N")$upper)
)
check$NPE_overlap <- check$NPE_upper >= check$target_lower & check$NPE_lower <= check$target_upper
check$MCMC_overlap <- check$MCMC_upper >= check$target_lower & check$MCMC_lower <= check$target_upper
write.csv(check, file.path(dirname(out_file), "calibration_interval_overlap.csv"), row.names = FALSE)
message("Wrote ", out_file, "; NPE draws used: ",
        sum(as.logical(npe$equilibrium_pass)), "; MCMC draws used: ",
        sum(as.logical(mcmc$equilibrium_pass)))
