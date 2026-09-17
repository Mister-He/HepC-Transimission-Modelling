# =============================================================================
# run_analysis_ver2.R — stratum-specific HCV elimination scenarios
#
# Run from the project root:
#   Rscript run_analysis_ver2.R
#
# This pipeline is independent of scripts/run_analysis.R.  It uses the
# calibrated model and posterior draws, but expands treatment allocation to a
# 4 x 4 stratum-by-stage hazard matrix through the optional tau_by_stratum
# parameter in sim.cpp.
# =============================================================================

suppressPackageStartupMessages({
  library(Rcpp)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(name, default = NULL) {
  pos <- match(name, args)
  if (is.na(pos)) return(default)
  args[pos + 1]
}

ROOT <- normalizePath(arg_val("--root", getwd()), mustWork = TRUE)
FIT_PATH <- normalizePath(
  arg_val("--fit", file.path(ROOT, "output", "calibration", "run1_4strata", "fit.rds")),
  mustWork = TRUE
)
POST_PATH <- normalizePath(
  arg_val("--posterior", file.path(ROOT, "output", "calibration", "npe_bayes",
                                    "posterior_samples_mcmc.csv")),
  mustWork = TRUE
)
N_DRAWS <- as.integer(arg_val("--n-draws", "50"))
N_CORES <- as.integer(arg_val("--n-cores", "4"))
if (!is.finite(N_DRAWS) || N_DRAWS < 2) stop("--n-draws must be at least 2")
if (!is.finite(N_CORES) || N_CORES < 1) stop("--n-cores must be at least 1")

FIG_DIR <- file.path(ROOT, "figures", "analysis_ver2")
OUT_DIR <- file.path(ROOT, "output", "analysis_ver2")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
PLOT_ONLY <- identical(tolower(arg_val("--plot-only", "false")), "true")

# -----------------------------------------------------------------------------
# 1. Setup / load model and posterior draws
# -----------------------------------------------------------------------------

if (!PLOT_ONLY) {
old_wd <- setwd(file.path(ROOT, "src"))
setup_env <- new.env(parent = globalenv())
sys.source("setup.R", envir = setup_env)
setwd(old_wd)
if (!exists("run_sim", envir = globalenv(), inherits = FALSE)) {
  Rcpp::sourceCpp(file.path(ROOT, "src", "sim.cpp"))
}
run_sim <- get("run_sim", envir = globalenv(), inherits = FALSE)
for (f in c("targets.R", "model_metrics.R", "equilibrium.R",
            "likelihood.R", "calibrate_nm.R")) {
  sys.source(file.path(ROOT, "src", "calibration", f), envir = environment())
}

fit <- readRDS(FIT_PATH)
post <- read.csv(POST_PATH, check.names = FALSE)
theta_names <- paste0("theta", 1:12)
if (!all(theta_names %in% names(post))) stop("Posterior file lacks theta1:theta12")
theta_all <- as.matrix(post[, theta_names])
if (any(!is.finite(theta_all))) stop("Posterior contains NA/NaN/Inf values")
set.seed(20260908)
draw_index <- sample(seq_len(nrow(theta_all)), min(N_DRAWS, nrow(theta_all)))
theta_draws <- theta_all[draw_index, , drop = FALSE]

eq_data <- setup_env$data
eq_data$t_start <- fit$t_start
eq_data$t_end <- fit$t_end

comp_cols <- function(s, k, h) s * 96 + (k - 1) * 24 + h * 6 + (0:5) + 2L
hcv_cols_by_stratum <- lapply(0:3, function(s) {
  unlist(lapply(1:4, function(k) {
    unlist(lapply(1:3, function(h) comp_cols(s, k, h)))
  }))
})
hcv_cols <- unlist(hcv_cols_by_stratum)

equilibrium_one <- function(theta) {
  pm <- tryCatch(build_params(theta, setup_env$params), error = function(e) NULL)
  if (is.null(pm)) return(NULL)
  out <- tryCatch(run_sim(pm, eq_data), error = function(e) NULL)
  if (is.null(out) || nrow(out) < 2 || any(!is.finite(out))) return(NULL)
  list(pm = pm, y0 = out[nrow(out), 2:385])
}

message("Equilibrating ", nrow(theta_draws), " posterior draws...")
eq_states <- if (N_CORES > 1 && requireNamespace("parallel", quietly = TRUE)) {
  parallel::mclapply(seq_len(nrow(theta_draws)), function(j) {
    equilibrium_one(theta_draws[j, ])
  }, mc.cores = N_CORES)
} else {
  lapply(seq_len(nrow(theta_draws)), function(j) equilibrium_one(theta_draws[j, ]))
}
valid_eq <- !vapply(eq_states, is.null, logical(1))
if (sum(valid_eq) < 2) stop("Fewer than two posterior draws reached equilibrium")
if (any(!valid_eq)) warning(sum(!valid_eq), " equilibrium draws excluded")

# -----------------------------------------------------------------------------
# 2. Scenario registry
# -----------------------------------------------------------------------------

stage_names <- c("NC", "CC", "DC", "HCC")
stratum_names <- c("D", "J", "F", "X")
outcome_names <- c("D", "J", "F", "X", "Total")
age_group_labels <- c("<20", "20-29", "30-39", "40-49", "50-59", "60+")

# Every treatment design uses the same 12-point C x J grid.  Community rates
# are deliberately lower, while prison uptake is capped at 100%/y.
community_levels <- c(0.01, 0.05, 0.10, 0.15)
prison_levels <- c(0.40, 0.70, 1.00)
rate_grid <- expand.grid(C_rate = community_levels, J_rate = prison_levels,
                         KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE) %>%
  arrange(match(C_rate, community_levels), match(J_rate, prison_levels))
n_rate_levels <- nrow(rate_grid)

rate_label <- function(x) paste0(trimws(formatC(100 * x, format = "fg", digits = 3)), "%/y")
rate_palette <- function(c_rate, j_rate) {
  c_key <- trimws(formatC(c_rate, format = "fg", digits = 3))
  endpoints <- switch(c_key,
    "0.01" = c("#e5f5e0", "#238b45"),
    "0.05" = c("#efedf5", "#54278f"),
    "0.1" = c("#deebf7", "#08519c"),
    "0.15" = c("#fee0d2", "#a50f15"),
    stop("Unexpected community uptake: ", c_rate)
  )
  grDevices::colorRampPalette(endpoints)(length(prison_levels))[match(j_rate, prison_levels)]
}
grid_palette <- vapply(seq_len(nrow(rate_grid)), function(i) {
  rate_palette(rate_grid$C_rate[i], rate_grid$J_rate[i])
}, character(1))
age_label <- function(min_age_group) age_group_labels[min_age_group]

# Matrix rows are D/J/F/X and columns are NC/CC/DC/HCC.
matrix_from_rates <- function(D = rep(0, 4), J = rep(0, 4),
                              F = rep(0, 4), X = rep(0, 4)) {
  out <- rbind(D = D, J = J, F = F, X = X)
  if (!all(dim(out) == c(4, 4)) || any(!is.finite(out)) || any(out < 0)) {
    stop("Invalid stratum-by-stage treatment matrix")
  }
  out
}
stage_rate_vector <- function(base, multiplier) {
  base * multiplier / max(multiplier)
}

make_scenario <- function(
    scenario_id, scenario_name, scenario_family, scenario_design,
    description, tau_matrix, coverage_rank, coverage_label,
    treatment_regimen = "SOF/VEL+RBV", min_age_group = 1L,
    alpha_DC_mode = "pos", phi_mode = "baseline",
    target_panel = "all outcomes", combination_strategy = "none",
    age_specific_design = NULL, disease_stage_specific_design = NULL) {
  stopifnot(is.matrix(tau_matrix), all(dim(tau_matrix) == c(4, 4)))
  if (any(!is.finite(tau_matrix)) || any(tau_matrix < 0)) {
    stop("Invalid treatment matrix in ", scenario_id)
  }
  if (min_age_group < 1 || min_age_group > 6) stop("Invalid minimum age")
  if (coverage_rank < 1 || coverage_rank > n_rate_levels) stop("Invalid coverage rank")
  values <- as.numeric(t(tau_matrix))
  names(values) <- as.vector(t(outer(stratum_names, stage_names, paste, sep = "_")))
  if (is.null(age_specific_design)) age_specific_design <- paste0(age_label(min_age_group), " and older")
  if (is.null(disease_stage_specific_design)) disease_stage_specific_design <- paste(stage_names, collapse = "=")
  row <- data.frame(
    scenario_id = scenario_id,
    scenario_name = scenario_name,
    scenario_family = scenario_family,
    scenario_design = scenario_design,
    description = description,
    treatment_regimen = treatment_regimen,
    overall_treatment_coverage = mean(tau_matrix),
    targeted_treatment_coverage = ifelse(any(tau_matrix > 0),
                                          mean(tau_matrix[tau_matrix > 0]), 0),
    coverage_rank = coverage_rank,
    coverage_label = coverage_label,
    age_specific_design = age_specific_design,
    disease_stage_specific_design = disease_stage_specific_design,
    other_stratum_specific_design = target_panel,
    target_panel = target_panel,
    combination_strategy = combination_strategy,
    min_age_group = as.integer(min_age_group),
    min_age_label = age_label(min_age_group),
    alpha_DC_mode = alpha_DC_mode,
    phi_mode = phi_mode,
    stringsAsFactors = FALSE
  )
  for (nm in names(values)) row[[paste0("tau_", nm)]] <- values[[nm]]
  row
}

rows <- list()
add <- function(...) rows[[length(rows) + 1L]] <<- make_scenario(...)
zero_matrix <- matrix(0, nrow = 4, ncol = 4,
                      dimnames = list(stratum_names, stage_names))

# Reference scenario.
add("status_quo", "Status quo", "reference", "Status quo",
    "No DAA treatment.", zero_matrix, 1, "Status quo",
    treatment_regimen = "None", target_panel = "all outcomes")

# Family A: all 12 C x J combinations. Each column represents a community
# target supported by prison treatment, while J represents prison treatment
# supported by community treatment in D/F/X.
for (i in seq_len(n_rate_levels)) {
  c_rate <- rate_grid$C_rate[i]
  j_rate <- rate_grid$J_rate[i]
  label <- paste0("C=", rate_label(c_rate), "; J=", rate_label(j_rate))
  add(paste0("target_D_", sprintf("%02d", i)),
      paste0("D-focused, ", label), "target_coverage", "D",
      "Treat D in the community and J in prison.",
      matrix_from_rates(D = rep(c_rate, 4), J = rep(j_rate, 4)), i, label,
      target_panel = "D")
  add(paste0("target_J_", sprintf("%02d", i)),
      paste0("J-focused, ", label), "target_coverage", "J",
      "Treat J in prison and D/F/X in the community.",
      matrix_from_rates(D = rep(c_rate, 4), J = rep(j_rate, 4),
                        F = rep(c_rate, 4), X = rep(c_rate, 4)), i, label,
      target_panel = "J")
  add(paste0("target_F_", sprintf("%02d", i)),
      paste0("F-focused, ", label), "target_coverage", "F",
      "Treat F in the community and J in prison.",
      matrix_from_rates(F = rep(c_rate, 4), J = rep(j_rate, 4)), i, label,
      target_panel = "F")
  add(paste0("target_X_", sprintf("%02d", i)),
      paste0("X-focused, ", label), "target_coverage", "X",
      "Treat X in the community and J in prison.",
      matrix_from_rates(X = rep(c_rate, 4), J = rep(j_rate, 4)), i, label,
      target_panel = "X")
  add(paste0("target_total_", sprintf("%02d", i)),
      paste0("Total mixed, ", label),
      "target_coverage", "Total mixed",
      "Treat all strata with the common C x J uptake combination.",
      matrix_from_rates(D = rep(c_rate, 4), J = rep(j_rate, 4),
                        F = rep(c_rate, 4), X = rep(c_rate, 4)), i,
      label, target_panel = "Total")
}

# Family B: age targeting, using the model's supported minimum-age rule.
age_specs <- list(
  age_40plus = list(label = "Age 40+", min_age = 4L),
  age_50plus = list(label = "Age 50+", min_age = 5L),
  age_60plus = list(label = "Age 60+", min_age = 6L)
)
for (design in names(age_specs)) {
  spec <- age_specs[[design]]
  for (i in seq_len(n_rate_levels)) {
    c_rate <- rate_grid$C_rate[i]; j_rate <- rate_grid$J_rate[i]
    add(paste0(design, "_", sprintf("%02d", i)),
        paste0(spec$label, ", C=", rate_label(c_rate), ", J=", rate_label(j_rate)),
        "age_targeting", spec$label,
        paste0("Mixed community/prison treatment restricted to ", spec$label, "."),
        matrix_from_rates(D = rep(c_rate, 4), J = rep(j_rate, 4),
                          F = rep(c_rate, 4), X = rep(c_rate, 4)), i,
        paste0("C=", rate_label(c_rate), "; J=", rate_label(j_rate)),
        min_age_group = spec$min_age, target_panel = "all outcomes")
  }
}

# Family C: ten levels for three disease-stage allocation designs.
stage_specs <- list(
  stage_equal = list(label = "Equal stage coverage", mult = c(1, 1, 1, 1)),
  stage_transmission = list(label = "Transmission priority", mult = c(1.4, 1.4, 0.6, 0.4)),
  stage_advanced = list(label = "Advanced-disease priority", mult = c(0.5, 0.7, 1.3, 1.5))
)
for (design in names(stage_specs)) {
  spec <- stage_specs[[design]]
  for (i in seq_len(n_rate_levels)) {
    c_base <- rate_grid$C_rate[i]
    j_base <- rate_grid$J_rate[i]
    c_rate <- stage_rate_vector(c_base, spec$mult)
    j_rate <- stage_rate_vector(j_base, spec$mult)
    add(paste0(design, "_", sprintf("%02d", i)),
        paste0(spec$label, ", level ", i), "stage_targeting", spec$label,
        paste0(spec$label, ": stage multipliers NC/CC/DC/HCC = ",
               paste(spec$mult, collapse = "/"), "."),
        matrix_from_rates(D = c_rate, J = j_rate, F = c_rate, X = c_rate), i,
        paste0("C=", rate_label(c_base), "; J=", rate_label(j_base)), target_panel = "all outcomes",
        disease_stage_specific_design = paste(stage_names, collapse = "="))
  }
}

# Family D: ten matched levels for two regimen settings already implemented in
# the model. Only the DC alpha switch differs; alpha_HCC remains fixed.
regimen_specs <- list(
  regimen_sofvel_rbv = list(label = "SOF/VEL+RBV", alpha = "pos"),
  regimen_sofvel = list(label = "SOF/VEL", alpha = "neg")
)
for (design in names(regimen_specs)) {
  spec <- regimen_specs[[design]]
  for (i in seq_len(n_rate_levels)) {
    c_rate <- rate_grid$C_rate[i]; j_rate <- rate_grid$J_rate[i]
    add(paste0(design, "_", sprintf("%02d", i)),
        paste0(spec$label, ", C=", rate_label(c_rate), ", J=", rate_label(j_rate)),
        "regimen", spec$label,
        paste0(spec$label, " with stage-matched setting-specific uptake."),
        matrix_from_rates(D = rep(c_rate, 4), J = rep(j_rate, 4),
                          F = rep(c_rate, 4), X = rep(c_rate, 4)), i,
        paste0("C=", rate_label(c_rate), "; J=", rate_label(j_rate)),
        treatment_regimen = spec$label, alpha_DC_mode = spec$alpha,
        target_panel = "all outcomes")
  }
}

# Family E: six selected combinations, each with all 12 C x J levels. These combine
# strong prison delivery with realistic weak community delivery, then add
# stage- or age-prioritisation in only a few interpretable policies.
combo_specs <- list(
  combo_J_DF = list(
    label = "J + D/F", combination = "strong prison J plus weak community D/F",
    make = function(c_rate, j_rate) matrix_from_rates(
      D = rep(c_rate, 4), J = rep(j_rate, 4), F = rep(c_rate, 4))),
  combo_J_F = list(
    label = "J + F", combination = "strong prison J plus post-release/community F",
    make = function(c_rate, j_rate) matrix_from_rates(
      J = rep(j_rate, 4), F = rep(c_rate, 4))),
  combo_J_DFX = list(
    label = "J + D/F/X", combination = "strong prison J plus weak community D/F/X",
    make = function(c_rate, j_rate) matrix_from_rates(
      D = rep(c_rate, 4), J = rep(j_rate, 4), F = rep(c_rate, 4), X = rep(c_rate, 4))),
  combo_J_DFX_early = list(
    label = "J + D/F/X early-stage", combination = "mixed setting plus transmission-focused stage allocation",
    make = function(c_rate, j_rate) matrix_from_rates(
      D = stage_rate_vector(c_rate, c(1.4, 1.4, 0.6, 0.4)),
      J = stage_rate_vector(j_rate, c(1.4, 1.4, 0.6, 0.4)),
      F = stage_rate_vector(c_rate, c(1.4, 1.4, 0.6, 0.4)),
      X = stage_rate_vector(c_rate, c(1.4, 1.4, 0.6, 0.4)))),
  combo_J_DFX_age40 = list(
    label = "J + D/F/X, age 40+", combination = "mixed setting plus age-40+ prioritisation",
    min_age = 4L,
    make = function(c_rate, j_rate) matrix_from_rates(
      D = rep(c_rate, 4), J = rep(j_rate, 4), F = rep(c_rate, 4), X = rep(c_rate, 4))),
  combo_J_DFX_age20_early = list(
    label = "J + D/F/X, age 20+ early-stage", combination = "mixed setting plus age-20+ and transmission-focused stage allocation",
    min_age = 2L,
    make = function(c_rate, j_rate) matrix_from_rates(
      D = stage_rate_vector(c_rate, c(1.4, 1.4, 0.6, 0.4)),
      J = stage_rate_vector(j_rate, c(1.4, 1.4, 0.6, 0.4)),
      F = stage_rate_vector(c_rate, c(1.4, 1.4, 0.6, 0.4)),
      X = stage_rate_vector(c_rate, c(1.4, 1.4, 0.6, 0.4))))
)
for (design in names(combo_specs)) {
  spec <- combo_specs[[design]]
  for (i in seq_len(n_rate_levels)) {
    c_rate <- rate_grid$C_rate[i]
    j_rate <- rate_grid$J_rate[i]
    min_age <- if (is.null(spec$min_age)) 1L else spec$min_age
    add(paste0(design, "_", sprintf("%02d", i)),
        paste0(spec$label, ", level ", i), "combination", spec$label,
        spec$combination, spec$make(c_rate, j_rate), i,
        paste0("C=", rate_label(c_rate), "; J=", rate_label(j_rate)),
        min_age_group = min_age, target_panel = "all outcomes",
        combination_strategy = spec$combination)
  }
}

registry <- bind_rows(rows)
if (anyDuplicated(registry$scenario_id)) stop("Duplicate scenario IDs")
if (any(!nzchar(registry$scenario_id))) stop("Empty scenario ID")

        level_palette <- grid_palette
registry$figure_colour <- NA_character_
registry$figure_colour[registry$scenario_id == "status_quo"] <- "#222222"
registry$figure_colour[registry$scenario_id != "status_quo"] <-
  level_palette[registry$coverage_rank[registry$scenario_id != "status_quo"]]
metadata_path <- file.path(ROOT, "output", "strategy_ver2.csv")
write.csv(registry, metadata_path, row.names = FALSE, na = "")

scenario_from_row <- function(row) {
  tau_cols <- as.vector(t(outer(stratum_names, stage_names, paste, sep = "_")))
  tau_matrix <- matrix(as.numeric(row[paste0("tau_", tau_cols)]),
                       nrow = 4, byrow = TRUE,
                       dimnames = list(stratum_names, stage_names))
  list(
    id = row$scenario_id,
    tau = colMeans(tau_matrix),
    tau_by_stratum = as.vector(t(tau_matrix)),
    tau_stratum = rep(1, 4),
    tau_min_age = as.integer(row$min_age_group),
    alpha_DC = if (row$alpha_DC_mode == "pos") setup_env$params$alpha_DC_pos
      else setup_env$params$alpha_DC_neg,
    rho = setup_env$params$rho,
    phi = c(setup_env$params$phi_CC_DC, setup_env$params$phi_CC_HCC,
            setup_env$params$phi_DC_HCC)
  )
}
scenario_objects <- lapply(seq_len(nrow(registry)), function(i) {
  scenario_from_row(registry[i, , drop = FALSE])
})

# -----------------------------------------------------------------------------
# 3. Run scenarios and calculate posterior medians
# -----------------------------------------------------------------------------

project_one <- function(eq_state, scenario) {
  if (is.null(eq_state)) return(NULL)
  pm <- eq_state$pm
  pm$tau <- scenario$tau
  pm$tau_by_stratum <- scenario$tau_by_stratum
  pm$tau_stratum <- scenario$tau_stratum
  pm$tau_min_age <- scenario$tau_min_age
  pm$rho <- scenario$rho
  pm$alpha_DC_pos <- scenario$alpha_DC
  pm$phi_CC_DC <- scenario$phi[1]
  pm$phi_CC_HCC <- scenario$phi[2]
  pm$phi_DC_HCC <- scenario$phi[3]
  proj_data <- eq_data
  proj_data$t_start <- 47
  proj_data$t_end <- 97
  proj_data$y0 <- eq_state$y0
  out <- tryCatch(run_sim(pm, proj_data), error = function(e) NULL)
  if (is.null(out) || any(!is.finite(out))) return(NULL)
  years <- 0:50
  row_index <- vapply(years, function(yr) which.min(abs(out[, 1] - (47 + yr))), integer(1))
  result <- data.frame(year_since_2017 = years, year = 2017 + years)
  for (s in seq_along(stratum_names)) {
    result[[paste0("hcv_", stratum_names[s])]] <-
      rowSums(out[row_index, hcv_cols_by_stratum[[s]], drop = FALSE])
  }
  result$hcv_Total <- rowSums(out[row_index, hcv_cols, drop = FALSE])
  result
}

run_one_scenario <- function(sc) {
  one_draw <- function(j) {
    out <- project_one(eq_states[[j]], sc)
    if (is.null(out)) return(NULL)
    out$draw <- j
    out
  }
  pieces <- if (N_CORES > 1 && requireNamespace("parallel", quietly = TRUE)) {
    parallel::mclapply(seq_along(eq_states), one_draw, mc.cores = N_CORES)
  } else {
    lapply(seq_along(eq_states), one_draw)
  }
  pieces <- pieces[!vapply(pieces, is.null, logical(1))]
  if (length(pieces) < 2) stop("Fewer than two valid runs for ", sc$id)
  bind_rows(pieces) %>% mutate(scenario_id = sc$id)
}

message("Running ", nrow(registry), " scenarios x ", sum(valid_eq),
        " valid posterior draws...")
raw_results <- bind_rows(lapply(seq_along(scenario_objects), function(i) {
  sc <- scenario_objects[[i]]
  message("  ", i, "/", length(scenario_objects), ": ", sc$id)
  run_one_scenario(sc)
}))
if (any(!is.finite(raw_results$hcv_Total)) || any(raw_results$hcv_Total < 0)) {
  stop("Invalid total HCV trajectory")
}

median_results <- raw_results %>%
  group_by(scenario_id, year_since_2017, year) %>%
  summarise(across(starts_with("hcv_"), median), .groups = "drop") %>%
  left_join(registry %>% select(scenario_id, scenario_name, scenario_family,
                                scenario_design, target_panel, coverage_rank,
                                coverage_label), by = "scenario_id") %>%
  arrange(scenario_family, scenario_design, coverage_rank, year_since_2017)
write.csv(median_results, file.path(OUT_DIR, "hcv_median_trajectories.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 4. Elimination assessment
# -----------------------------------------------------------------------------

ELIMINATION_THRESHOLD <- 1
first_year_at_threshold <- function(years, values) {
  hit <- which(values <= ELIMINATION_THRESHOLD)
  if (!length(hit)) return(NA_integer_)
  years[hit[1]]
}
baseline <- median_results %>% filter(scenario_id == "status_quo")
baseline_10 <- baseline$hcv_Total[baseline$year_since_2017 == 10]
baseline_20 <- baseline$hcv_Total[baseline$year_since_2017 == 20]

elimination <- median_results %>%
  group_by(scenario_id, scenario_name, scenario_family, scenario_design) %>%
  summarise(
    hcv_2017 = hcv_Total[year_since_2017 == 0],
    hcv_10y = hcv_Total[year_since_2017 == 10],
    hcv_20y = hcv_Total[year_since_2017 == 20],
    hcv_50y = hcv_Total[year_since_2017 == 50],
    reduction_10y_vs_status_quo = 1 - hcv_10y / baseline_10,
    reduction_20y_vs_status_quo = 1 - hcv_20y / baseline_20,
    minimum_hcv = min(hcv_Total),
    minimum_year = year[which.min(hcv_Total)],
    first_year_hcv_le_1_since_2017 = first_year_at_threshold(year_since_2017, hcv_Total),
    approximately_zero_by_10y = hcv_10y <= ELIMINATION_THRESHOLD,
    approximately_zero_by_20y = hcv_20y <= ELIMINATION_THRESHOLD,
    .groups = "drop"
  ) %>%
  left_join(registry %>% select(-scenario_name, -scenario_family, -scenario_design),
            by = "scenario_id") %>%
  arrange(scenario_family, scenario_design, coverage_rank)
write.csv(elimination, file.path(ROOT, "output", "elimination_assessment_ver2.csv"),
          row.names = FALSE, na = "")
} else {
  # Reuse the completed posterior-median results when only the layout or
  # visual encoding is being revised; this avoids rerunning the ODE model.
  stage_names <- c("NC", "CC", "DC", "HCC")
  stratum_names <- c("D", "J", "F", "X")
  outcome_names <- c("D", "J", "F", "X", "Total")
  metadata_path <- file.path(ROOT, "output", "strategy_ver2.csv")
  trajectories_path <- file.path(OUT_DIR, "hcv_median_trajectories.csv")
  if (!file.exists(metadata_path) || !file.exists(trajectories_path)) {
    stop("--plot-only requires existing strategy_ver2.csv and hcv_median_trajectories.csv")
  }
  registry <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
  median_results <- read.csv(trajectories_path, stringsAsFactors = FALSE,
                             check.names = FALSE)
  if (!nrow(registry) || !nrow(median_results)) stop("Plot-only inputs are empty")
}

# -----------------------------------------------------------------------------
# 5. F2-style combined figures: columns are treatment designs and rows are the
# same outcome measure (total infections in D, J, F, X, and the total population).
# Each column is drawn as an independent ggplot. All columns use the same
# ordered treatment palette and their legends list all 12 C x J combinations.
# -----------------------------------------------------------------------------

theme_ver2 <- function(base_size = 28, show_y = TRUE, show_strips = TRUE) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(fill = NA, colour = "grey55", linewidth = 0.55),
      axis.ticks = element_line(colour = "grey35", linewidth = 0.45),
      axis.title = element_text(size = 34, face = "bold"),
      axis.text = element_text(size = 29, face = "bold", colour = "black"),
      axis.text.x = element_text(size = 29, angle = 0, hjust = 0.5,
                                 face = "bold", colour = "black"),
      plot.title = element_text(size = 34, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 25, hjust = 0),
      strip.background = element_rect(fill = "white", colour = NA),
      strip.placement = "outside",
      legend.position = "bottom",
      legend.justification = "left",
      legend.box.just = "left",
      legend.title = element_text(size = 25, face = "bold"),
      legend.text = element_text(size = 21, face = "bold"),
      legend.margin = margin(0, 0, 0, 0),
      legend.spacing.y = grid::unit(0.02, "lines"),
      plot.margin = margin(10, 10, 8, 10),
      panel.spacing = grid::unit(1.2, "lines"),
      axis.title.y = if (show_y) element_text(size = 34, face = "bold",
                                               margin = margin(r = 5)) else element_blank(),
      axis.text.y = if (show_y) element_text(size = 29, face = "bold",
                                              colour = "black", margin = margin(r = 2)) else element_blank(),
      axis.ticks.y = if (show_y) element_line(colour = "grey35", linewidth = 0.45) else element_blank(),
      strip.text.y.left = if (show_strips) element_text(size = 32, face = "bold", angle = 0,
                                                        colour = "black")
        else element_blank()
    )
}

long_results <- median_results %>%
  pivot_longer(starts_with("hcv_"), names_to = "outcome", values_to = "hcv") %>%
  mutate(outcome = sub("^hcv_", "", outcome),
         outcome = factor(outcome, levels = outcome_names),
         line_label = ifelse(scenario_id == "status_quo", "Status quo", coverage_label))

if (!exists("n_rate_levels")) {
  n_rate_levels <- max(registry$coverage_rank[registry$scenario_id != "status_quo"],
                       na.rm = TRUE)
}
shared_treatment_palette <- registry %>%
  filter(scenario_id != "status_quo") %>%
  distinct(coverage_rank, figure_colour) %>%
  arrange(coverage_rank) %>%
  pull(figure_colour)
column_palette <- function(column_key = NULL, n = n_rate_levels) {
  if (n == n_rate_levels) return(shared_treatment_palette)
  shared_treatment_palette[seq_len(n)]
}

column_legend_title <- function(family, column_key) {
  if (family == "target_coverage") {
    if (column_key == "Total") return("C (D/F/X); J\n(%/year)")
    if (column_key == "J") return("Prison\n(%/year)")
    return("Community\n(%/year)")
  }
  if (family %in% c("age_targeting", "stage_targeting", "regimen", "combination")) {
    return("C; J uptake\n(%/year)")
  }
  "Treatment level"
}

make_column_plot <- function(family, column_key, column_index, target_family = FALSE,
                             title_override = NULL, scenario_ids = NULL,
                             line_colours = NULL) {
  fam_meta <- registry %>% filter(scenario_family == family)
  col_meta <- if (is.null(scenario_ids)) {
    if (target_family) fam_meta %>% filter(target_panel == column_key)
    else fam_meta %>% filter(scenario_design == column_key)
  } else {
    registry %>% filter(scenario_id %in% scenario_ids)
  }
  if (!nrow(col_meta)) stop("No scenarios found for ", family, ": ", column_key)
  if (!is.null(scenario_ids)) {
    col_meta <- col_meta %>% mutate(.order = match(scenario_id, scenario_ids)) %>%
      arrange(.order) %>% select(-.order)
  } else {
    col_meta <- col_meta %>% arrange(coverage_rank)
  }
  intervention_ids <- col_meta$scenario_id
  plot_data <- bind_rows(
    long_results %>% filter(scenario_id == "status_quo") %>% mutate(line_label = "Status quo"),
    long_results %>% filter(scenario_id %in% intervention_ids)
  ) %>% mutate(outcome = factor(outcome, levels = outcome_names))

  if (is.null(scenario_ids)) {
    intervention_labels <- unique(col_meta$coverage_label)
    colour_values <- c("Status quo" = "#202020",
                       setNames(column_palette(column_key, length(intervention_labels)),
                                intervention_labels))
    colour_breaks <- c("Status quo", intervention_labels)
  } else {
    intervention_labels <- col_meta$scenario_id
    label_map <- setNames(paste0(col_meta$scenario_design, ": ",
                                 format_combined_rate_label(col_meta$coverage_label)),
                          col_meta$scenario_id)
    plot_data$line_label <- ifelse(plot_data$scenario_id == "status_quo",
                                   "Status quo", unname(label_map[plot_data$scenario_id]))
    if (is.null(line_colours)) {
      line_colours <- grDevices::hcl.colors(length(intervention_ids), "Dark 3")
    }
    colour_values <- c("Status quo" = "#202020",
                       setNames(line_colours, unname(label_map[intervention_ids])))
    colour_breaks <- c("Status quo", unname(label_map[intervention_ids]))
  }
  show_y <- column_index == 1L
  show_strips <- column_index == 1L

  plot_object <- ggplot(plot_data, aes(x = year_since_2017, y = hcv / 1000,
                                       group = scenario_id, colour = line_label)) +
    geom_vline(xintercept = c(10, 20), colour = "grey72", linewidth = 0.55,
               linetype = "dashed") +
    geom_hline(yintercept = 0, colour = "grey35", linewidth = 0.45) +
    geom_line(linewidth = 1.6) +
    facet_grid(outcome ~ ., scales = "free_y", switch = "y",
               labeller = labeller(outcome = function(x) x)) +
    scale_colour_manual(values = colour_values, breaks = colour_breaks,
                        name = column_legend_title(family, column_key)) +
    scale_x_continuous(breaks = seq(0, 50, 10), limits = c(0, 50),
                       expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1),
                       expand = expansion(mult = c(0.02, 0.08))) +
    labs(title = if (is.null(title_override)) column_key else title_override, x = NULL,
         y = if (show_y) "Total infections (000s)" else NULL) +
    guides(colour = "none") +
    theme_ver2(base_size = 28, show_y = show_y, show_strips = show_strips) +
    theme(legend.position = "none")
  # With strip.placement = "outside", ggplot2 places the row labels first and
  # the y-axis coordinates immediately before the panel.  Keep that native
  # order so the coordinates remain adjacent to the leftmost subfigures.
  plot_object
}

format_combined_rate_label <- function(x) {
  gsub("=", ": ", gsub("; ", " + ", x, fixed = TRUE), fixed = TRUE)
}

make_shared_legend <- function() {
  legend_meta <- registry %>%
    filter(scenario_id != "status_quo") %>%
    distinct(coverage_rank, coverage_label) %>%
    arrange(coverage_rank)
  labels <- format_combined_rate_label(legend_meta$coverage_label)
  legend_spacing <- 1 / 6
  legend_center <- 3.6
  row_y <- legend_center +
    (seq(length(labels), 1) - mean(seq(length(labels), 1))) * legend_spacing
  entries <- data.frame(
    x = rep(0.15, length(labels)), xend = rep(0.95, length(labels)),
    y = row_y, yend = row_y, colour = shared_treatment_palette,
    label = labels, stringsAsFactors = FALSE
  )
  ggplot() +
    geom_text(data = data.frame(x = 0, y = max(row_y) + 0.6,
                                label = "Treatment scenarios"),
              aes(x = x, y = y, label = label), hjust = 0,
              size = 11.4, fontface = "bold") +
    geom_segment(data = entries,
                 aes(x = x, xend = xend, y = y, yend = y, colour = colour),
                 linewidth = 3.2, lineend = "butt") +
    geom_text(data = entries, aes(x = 1.25, y = y, label = label),
              hjust = 0, size = 9.84, fontface = "bold") +
    scale_colour_identity() +
    coord_cartesian(xlim = c(0, 8.5), ylim = c(0, 7.2),
                    expand = FALSE, clip = "off") +
    theme_void() +
    theme(plot.margin = margin(18, 0, 18, 2))
}

# The best-treatment figure has one selected scenario per outcome column, so
# its legend explains the small set of selected designs instead of all 12 rates.
make_best_legend <- function(best_table, best_colours) {
  legend_meta <- best_table %>%
    distinct(best_scenario_id, scenario_design, coverage_label, description) %>%
    mutate(description = gsub("strong prison J plus weak community",
                              "prison J + community", description),
           description = gsub("strong prison J plus post-release/community",
                              "prison J + community", description),
           label = paste0(scenario_design, " | ",
                          format_combined_rate_label(coverage_label),
                          " | ", description)) %>%
    arrange(match(best_scenario_id, names(best_colours)))
  legend_spacing <- 1 / 6
  legend_center <- 3.6
  row_y <- legend_center +
    (seq(nrow(legend_meta), 1) - mean(seq(nrow(legend_meta), 1))) * legend_spacing
  entries <- data.frame(
    x = rep(0.15, nrow(legend_meta)), xend = rep(0.95, nrow(legend_meta)),
    y = row_y, yend = row_y,
    colour = unname(best_colours[legend_meta$best_scenario_id]),
    label = legend_meta$label, stringsAsFactors = FALSE
  )
  ggplot() +
    geom_text(data = data.frame(x = 0, y = max(row_y) + 0.6,
                                label = "Selected treatment designs"),
              aes(x = x, y = y, label = label), hjust = 0,
              size = 11.4, fontface = "bold") +
    geom_segment(data = entries,
                 aes(x = x, xend = xend, y = y, yend = y, colour = colour),
                 linewidth = 3.2, lineend = "butt") +
    geom_text(data = entries, aes(x = 1.25, y = y, label = label),
              hjust = 0, size = 9.84, fontface = "bold") +
    scale_colour_identity() +
    coord_cartesian(xlim = c(0, 18), ylim = c(0, 7.2),
                    expand = FALSE, clip = "off") +
    theme_void() +
    theme(plot.margin = margin(18, 0, 18, 2))
}

# Select one treatment scenario for each outcome using the lowest median total
# infection count at year 20, with year 10 as a tie-breaker. The resulting
# figure shows each selected scenario across all five outcome rows.
design_performance <- long_results %>%
  # target_coverage is an exploratory allocation family and is excluded from
  # the best-treatment comparison at the user's request.
  filter(scenario_id != "status_quo",
         scenario_family != "target_coverage",
         year_since_2017 %in% c(10, 20)) %>%
  group_by(scenario_id, scenario_family, scenario_design, outcome) %>%
  summarise(
    hcv_10 = hcv[year_since_2017 == 10][1],
    hcv_20 = hcv[year_since_2017 == 20][1],
    .groups = "drop"
  )
design_summary <- design_performance %>%
  group_by(scenario_family, scenario_design, outcome) %>%
  summarise(
    best_hcv_10 = min(hcv_10, na.rm = TRUE),
    best_hcv_20 = min(hcv_20, na.rm = TRUE),
    best_scenario_id = scenario_id[which.min(hcv_20)],
    .groups = "drop"
  )
best_designs <- design_summary %>%
  group_by(outcome) %>%
  arrange(best_hcv_20, best_hcv_10, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(factor(outcome, levels = outcome_names))
if (nrow(best_designs) != length(outcome_names)) {
  stop("Best-design selection did not produce one design per outcome")
}
best_design_registry <- best_designs %>%
  left_join(registry %>% select(scenario_id, scenario_name, coverage_label,
                                treatment_regimen, age_specific_design,
                                disease_stage_specific_design,
                                combination_strategy),
            by = c("best_scenario_id" = "scenario_id"))
write.csv(best_design_registry, file.path(OUT_DIR, "best_treatment_by_stratum.csv"),
          row.names = FALSE, na = "")

combine_family_plots <- function(family, columns, target_family = FALSE) {
  plots <- Map(function(column_key, column_index) {
    make_column_plot(family, column_key, column_index, target_family)
  }, columns, seq_along(columns))
  legend <- make_shared_legend()
  left_row <- wrap_plots(plots, nrow = 1, widths = rep(1, length(plots)),
                         guides = "keep")
  shared_x <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "Time (years)",
             size = 10.5, fontface = "bold") +
    theme_void()
  left_with_x <- wrap_plots(left_row, shared_x, ncol = 1,
                            heights = c(1, 0.075))
  wrap_plots(list(left_with_x, plot_spacer(), legend), nrow = 1,
             widths = c(length(plots), 0.12, 2.25), guides = "keep")
}

figure_specs <- c(
  target_coverage = "fig_ver2_target_coverage_by_stratum.png",
  age_targeting = "fig_ver2_age_targeting_by_stratum.png",
  stage_targeting = "fig_ver2_stage_targeting_by_stratum.png",
  regimen = "fig_ver2_regimen_by_stratum.png",
  combination = "fig_ver2_best_treatment_by_stratum.png"
)

target_columns <- c("D", "J", "F", "X", "Total")
ggsave(file.path(FIG_DIR, figure_specs[["target_coverage"]]),
       combine_family_plots("target_coverage", target_columns, target_family = TRUE),
       width = 42, height = 29, dpi = 300, bg = "white", limitsize = FALSE)

family_columns <- list(
  age_targeting = c("Age 40+", "Age 50+", "Age 60+"),
  stage_targeting = c("Equal stage coverage", "Transmission priority", "Advanced-disease priority"),
  regimen = c("SOF/VEL+RBV", "SOF/VEL")
)
family_widths <- c(age_targeting = 31, stage_targeting = 31, regimen = 25,
                   combination = 42)
family_heights <- c(age_targeting = 29, stage_targeting = 29, regimen = 29,
                    combination = 31)
for (family in names(family_columns)) {
  ggsave(file.path(FIG_DIR, figure_specs[[family]]),
         combine_family_plots(family, family_columns[[family]]),
         width = family_widths[[family]], height = family_heights[[family]],
         dpi = 300, bg = "white", limitsize = FALSE)
}
best_columns <- outcome_names
best_scenario_ids <- best_designs$best_scenario_id
best_scenario_ids_unique <- unique(best_scenario_ids)
best_selected_colours <- setNames(
  c("#006d2c", "#54278f", "#08519c", "#a50f15", "#8c510a")[seq_along(best_scenario_ids_unique)],
  best_scenario_ids_unique
)
best_plots <- Map(function(column_key, column_index) {
  choice <- best_designs %>% filter(outcome == column_key)
  make_column_plot(choice$scenario_family, choice$scenario_design,
                   column_index, target_family = FALSE,
                   title_override = choice$scenario_design,
                   scenario_ids = choice$best_scenario_id,
                   line_colours = unname(best_selected_colours[choice$best_scenario_id]))
}, best_columns, seq_along(best_columns))
best_legend <- make_best_legend(
  best_designs %>%
    left_join(registry %>% select(scenario_id, coverage_label, description),
              by = c("best_scenario_id" = "scenario_id")),
  best_selected_colours
)
best_left <- wrap_plots(best_plots, nrow = 1,
                        widths = rep(1, length(best_plots)), guides = "keep")
best_shared_x <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "Time (years)",
           size = 10.5, fontface = "bold") +
  theme_void()
best_left_with_x <- wrap_plots(best_left, best_shared_x, ncol = 1,
                               heights = c(1, 0.075))
ggsave(file.path(FIG_DIR, figure_specs[["combination"]]),
       wrap_plots(list(best_left_with_x, plot_spacer(), best_legend), nrow = 1,
                  widths = c(length(best_plots), 0.12, 2.25), guides = "keep"),
       width = 42, height = 29, dpi = 300, bg = "white", limitsize = FALSE)
write.csv(data.frame(
  figure_file = unname(figure_specs),
  scenario_family = names(figure_specs),
  stringsAsFactors = FALSE
), file.path(OUT_DIR, "figure_manifest.csv"), row.names = FALSE)

message("Scenarios evaluated: ", nrow(registry))
message("Valid posterior draws: ", if (PLOT_ONLY) "reused from existing results" else sum(valid_eq))
message("Scenario registry: ", metadata_path)
message("Figures written to: ", FIG_DIR)
message("Elimination assessment: ", file.path(ROOT, "output", "elimination_assessment_ver2.csv"))
