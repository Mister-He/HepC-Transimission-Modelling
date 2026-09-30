#!/usr/bin/env Rscript
# 2026 treatment scenarios. Treatment starts at model time 160 and are
# evaluated through 2030.

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) {
  i <- match(name, args); if (is.na(i)) default else args[i + 1]
}
ROOT <- normalizePath(arg("--root", getwd()), mustWork = TRUE)
FIT <- normalizePath(arg("--fit", file.path(ROOT,
  "output/calibration/r_20260929_agepop/fit.rds")), mustWork = TRUE)
POSTERIOR_DIR <- arg("--posterior-dir", NULL)
OUT <- file.path(ROOT, arg("--out-dir", "output/scenarios_2026"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

old <- setwd(file.path(ROOT, "scripts/model")); setup_env <- new.env(parent = globalenv())
sys.source("setup.R", setup_env); setwd(old)
Rcpp::sourceCpp(file.path(ROOT, "scripts/model/sim.cpp"))
source(file.path(ROOT, "scripts/calibration/targets.R"))
source(file.path(ROOT, "scripts/calibration/model_metrics.R"))
source(file.path(ROOT, "scripts/calibration/equilibrium.R"))
source(file.path(ROOT, "scripts/calibration/likelihood.R"))
targets <- read_targets(file.path(ROOT, "data/calibration_targets.json"))
fit <- readRDS(FIT)
if (!is.null(POSTERIOR_DIR)) {
  posterior_dir <- normalizePath(file.path(ROOT, POSTERIOR_DIR), mustWork = TRUE)
  posterior_files <- c(file.path(posterior_dir, "posterior_samples_npe.csv"),
                       file.path(posterior_dir, "posterior_samples_npe_seed2.csv"))
  posterior_files <- posterior_files[file.exists(posterior_files)]
  if (!length(posterior_files)) stop("No NPE posterior samples found in ", posterior_dir)
  posterior <- do.call(rbind, lapply(posterior_files, function(f)
    read.csv(f)[, paste0("theta", 1:12), drop = FALSE]))
  calibration_theta <- apply(posterior, 2, median)
  calibration_source <- paste0("NPE posterior median: ", paste(posterior_files, collapse = "; "))
} else {
  calibration_theta <- fit$best$theta
  calibration_source <- "deterministic fit best theta"
}
base <- build_params(calibration_theta, fit$base_params)
snapshot_data <- setup_env$data; snapshot_data$t_start <- 0; snapshot_data$t_end <- 150
snapshot <- run_sim(base, snapshot_data)
if (!check_equilibrium(snapshot)$pass) stop("2015 snapshot failed equilibrium check")

# Model age groups are <20, 20-29, 30-39, 40-49, 50-59 and 60+.
# The requested 20-50 and 30-60 settings therefore mean 20-49 and 30-59.
age_settings <- data.frame(
  age_group_setting = c("20_50", "30_60", "all_ages"),
  age_label = c("20-50", "30-60", "All age groups"),
  age_min_group = c(2L, 3L, 1L), age_max_group = c(4L, 5L, 6L),
  stringsAsFactors = FALSE)
probabilities <- c(0, .25, .50, .75, 1)
scenario_table <- merge(age_settings,
  expand.grid(initiation_prob_J = probabilities, initiation_prob_FX = probabilities,
              KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE),
  by = NULL, sort = FALSE)
scenario_table$initiation_prob_D <- 0
scenario_table$initiation_prob_F <- scenario_table$initiation_prob_FX
scenario_table$initiation_prob_X <- scenario_table$initiation_prob_FX
fmt_prob <- function(x) sprintf("%g", 100 * x)
scenario_table$scenario_id <- sprintf("%s_J%s_FX%s", scenario_table$age_group_setting,
  fmt_prob(scenario_table$initiation_prob_J), fmt_prob(scenario_table$initiation_prob_FX))
scenario_table$scenario <- ifelse(scenario_table$initiation_prob_J == 0 &
  scenario_table$initiation_prob_FX == 0, "status quo", scenario_table$scenario_id)
scenario_table$description <- sprintf("%s: J %s%%, F/X %s%% initiation probability",
  scenario_table$age_label, fmt_prob(scenario_table$initiation_prob_J),
  fmt_prob(scenario_table$initiation_prob_FX))
scenario_table$plain_language <- sprintf(
  "Among eligible chronic people aged %s, each year approximately %s%% in J and %s%% in F/X start a treatment course; D is not eligible.",
  scenario_table$age_label, fmt_prob(scenario_table$initiation_prob_J),
  fmt_prob(scenario_table$initiation_prob_FX))
scenario_table$scenario_type <- "age_probability_grid"
scenario_table <- scenario_table[order(scenario_table$age_group_setting,
  scenario_table$initiation_prob_J, scenario_table$initiation_prob_FX), ]
rownames(scenario_table) <- NULL

idx0 <- function(s, k, h, i) s * 120 + (k - 1) * 30 + h * 6 + i + 2L
all_idx <- function(strata, states = 0:4) as.vector(sapply(strata, function(s)
  sapply(1:4, function(k) sapply(states, function(h) sapply(0:5, function(i) idx0(s, k, h, i))))))
viremic_idx <- all_idx(0:3, 1:3); current_idx <- all_idx(c(0, 2), 1:3)
sero_idx <- all_idx(0:3, 1:4)

make_params <- function(row) {
  pm <- base
  probability_to_tau <- function(x) -log1p(-pmin(as.numeric(x), 0.99995))
  pm$tau <- rep(1, 4)
  pm$tau_stratum <- probability_to_tau(row[paste0("initiation_prob_", c("D", "J", "F", "X"))])
  pm$tau_min_age <- as.integer(row$age_min_group); pm$tau_max_age <- as.integer(row$age_max_group)
  # C++ indexes this vector as stratum-major; transpose before flattening.
  pm$tau_by_stratum <- as.vector(t(outer(pm$tau_stratum, pm$tau)))
  pm
}

rate_components <- function(state, time, pm) {
  y <- as.numeric(state); idxs <- function(s, k, h, age) idx0(s, k, h, age) - 1L
  incidence <- 0
  for (age in 0:5) {
    lambda <- 0
    for (source_age in 0:5) {
      infectious <- 0; active <- 0
      for (s in c(0, 2)) for (k in 1:4) {
        infectious <- infectious + y[idxs(s, k, 1, source_age)] + y[idxs(s, k, 2, source_age)]
        active <- active + sum(y[sapply(0:4, function(h) idxs(s, k, h, source_age))])
      }
      if (active > 0) lambda <- lambda + pm$C_contact[age + 1, source_age + 1] * infectious / active
    }
    m <- pm$m_min + (pm$m_max - pm$m_min) / (1 + exp(-(time - pm$m_t0) / pm$m_tau))
    gam <- pm$q * m * lambda
    for (s in c(0, 2)) for (k in 1:4)
      incidence <- incidence + gam * (y[idxs(s, k, 0, age)] + y[idxs(s, k, 4, age)])
  }
  mortality <- 0
  for (s in 0:3) for (k in 1:4) for (h in 0:4) for (age in 0:5) {
    disease_rate <- if (k == 3) pm$mu_DC * ifelse(h %in% c(3, 4), pm$psi_DC, 1) else
      if (k == 4) pm$mu_HCC * ifelse(h %in% c(3, 4), pm$psi_HCC, 1) else 0
    mortality <- mortality + disease_rate * y[idxs(s, k, h, age)]
  }
  c(incidence = incidence, mortality = mortality)
}

annual <- function(z, row, pm) {
  keep <- which(abs(z[, 1] - round(z[, 1])) < 1e-8); z <- z[keep, , drop = FALSE]
  total <- rowSums(z[, -1, drop = FALSE]); v <- rowSums(z[, viremic_idx, drop = FALSE])
  cv <- rowSums(z[, current_idx, drop = FALSE]); sero <- rowSums(z[, sero_idx, drop = FALSE])
  current_pop <- rowSums(z[, all_idx(c(0, 2)), drop = FALSE])
  rates <- do.call(rbind, lapply(seq_len(nrow(z)), function(i) rate_components(z[i, -1], z[i, 1], pm)))
  data.frame(year = round(2015 + z[, 1] - 150), scenario_id = row$scenario_id,
    scenario = row$scenario,
    total_viremic_prev = v / total, current_viremic_prev = cv / current_pop,
    seroprevalence = sero / total, current_population = current_pop,
    total_population = total, hcv_incidence_rate = rates[, "incidence"],
    hcv_mortality_rate = rates[, "mortality"])
}

flow_by_stratum <- function(z, pm) {
  flow <- matrix(0, nrow(z), 4, dimnames = list(NULL, c("D", "J", "F", "X")))
  for (s in 0:3) for (k in 1:4) for (age in 0:5)
    if (age + 1 >= pm$tau_min_age && age + 1 <= pm$tau_max_age)
      flow[, s + 1] <- flow[, s + 1] + pm$tau_by_stratum[s * 4 + k] * z[, idx0(s, k, 2, age)]
  flow
}
integrate_columns <- function(z, x) colSums((head(x, -1) + tail(x, -1)) * diff(z[, 1]) / 2)

treatment_summary <- function(z, pm, start_state) {
  starts <- integrate_columns(z, flow_by_stratum(z, pm))
  pop <- do.call(c, lapply(0:3, function(s) sum(start_state[all_idx(s) - 1L])))
  person_years <- integrate_columns(z, do.call(cbind, lapply(0:3, function(s)
    rowSums(z[, all_idx(s), drop = FALSE]))))
  names(pop) <- names(starts); names(person_years) <- names(starts)
  data.frame(treatment_starts_D = starts["D"], treatment_starts_J = starts["J"],
    treatment_starts_F = starts["F"], treatment_starts_X = starts["X"],
    treatment_starts_JFX = sum(starts[c("J", "F", "X")]), treatment_starts_total = sum(starts),
    population_2026_D = pop["D"], population_2026_J = pop["J"], population_2026_F = pop["F"],
    population_2026_X = pop["X"], population_2026_total = sum(pop),
    treatment_courses_per_2026_person_D = starts["D"] / pop["D"],
    treatment_courses_per_2026_person_J = starts["J"] / pop["J"],
    treatment_courses_per_2026_person_F = starts["F"] / pop["F"],
    treatment_courses_per_2026_person_X = starts["X"] / pop["X"],
    treatment_courses_per_2026_person_JFX = sum(starts[c("J", "F", "X")]) / sum(pop[c("J", "F", "X")]),
    treatment_courses_per_2026_person_total = sum(starts) / sum(pop),
    treatment_courses_per_person_year_total = sum(starts) / sum(person_years))
}
burden_summary <- function(z) {
  v <- rowSums(z[, viremic_idx, drop = FALSE])
  data.frame(cumulative_viremic_person_years = integrate_columns(z, matrix(v, ncol = 1))[[1]],
    mean_annual_viremic_population = mean(v[-1]))
}

run_one <- function(row) {
  pm <- make_params(row); pm0 <- pm
  pm0$tau <- rep(0, 4); pm0$tau_stratum <- rep(0, 4); pm0$tau_by_stratum <- rep(0, 16)
  d0 <- setup_env$data; d0$y0 <- snapshot[nrow(snapshot), -1]; d0$t_start <- 150; d0$t_end <- 160
  pre <- run_sim(pm0, d0)
  d1 <- d0; d1$y0 <- pre[nrow(pre), -1]; d1$t_start <- 160; d1$t_end <- 165
  treated <- run_sim(pm, d1)
  trajectory <- rbind(annual(pre, row, pm0), annual(treated, row, pm))
  b <- trajectory[trajectory$year == 2015, , drop = FALSE][1, ]; e <- trajectory[trajectory$year == 2030, , drop = FALSE][1, ]
  list(trajectory = trajectory, treatment = treatment_summary(treated, pm, pre[nrow(pre), -1]),
    burden = burden_summary(treated), rates = data.frame(incidence_2015 = b$hcv_incidence_rate,
      mortality_2015 = b$hcv_mortality_rate, incidence_2030 = e$hcv_incidence_rate,
      mortality_2030 = e$hcv_mortality_rate))
}

cat("Running", nrow(scenario_table), "age/probability combinations using", calibration_source, "...\n")
scenario_runs <- lapply(seq_len(nrow(scenario_table)), function(i) run_one(scenario_table[i, , drop = FALSE]))
trajectories <- do.call(rbind, lapply(seq_along(scenario_runs), function(i) {
  x <- scenario_runs[[i]]$trajectory; x$scenario_type <- scenario_table$scenario_type[i]; x
}))
scenario_meta <- scenario_table[, c("scenario_id", "scenario", "scenario_type", "age_group_setting", "age_label",
  "age_min_group", "age_max_group", "initiation_prob_D", "initiation_prob_J", "initiation_prob_F",
  "initiation_prob_X", "initiation_prob_FX", "description", "plain_language")]
trajectories <- dplyr::left_join(trajectories, scenario_meta,
  by = c("scenario_id", "scenario", "scenario_type"))
endpoint <- trajectories[trajectories$year == 2030, ]
treatment <- do.call(rbind, lapply(seq_along(scenario_runs), function(i) {
  x <- scenario_runs[[i]]$treatment; x$scenario_id <- scenario_table$scenario_id[i]
  x$scenario <- scenario_table$scenario[i]
  x$scenario_type <- scenario_table$scenario_type[i]; x
}))
burden <- do.call(rbind, lapply(seq_along(scenario_runs), function(i) {
  x <- scenario_runs[[i]]$burden; x$scenario_id <- scenario_table$scenario_id[i]
  x$scenario <- scenario_table$scenario[i]
  x$scenario_type <- scenario_table$scenario_type[i]; x
}))
rates <- do.call(rbind, lapply(seq_along(scenario_runs), function(i) {
  x <- scenario_runs[[i]]$rates; x$scenario_id <- scenario_table$scenario_id[i]
  x$scenario <- scenario_table$scenario[i]
  x$scenario_type <- scenario_table$scenario_type[i]; x
}))
trajectories <- dplyr::left_join(trajectories, rates,
  by = c("scenario_id", "scenario", "scenario_type"))
endpoint <- dplyr::left_join(endpoint, burden,
  by = c("scenario_id", "scenario", "scenario_type"))
endpoint <- dplyr::left_join(endpoint, treatment,
  by = c("scenario_id", "scenario", "scenario_type"))
endpoint <- dplyr::left_join(endpoint, rates,
  by = c("scenario_id", "scenario", "scenario_type"))
endpoint$viremic_prevalence_screen_pass <- endpoint$total_viremic_prev <= .01
endpoint$who_incidence_target_ratio <- .20; endpoint$who_mortality_target_ratio <- .35
endpoint$incidence_ratio_2030_vs_2015 <- endpoint$incidence_2030 / endpoint$incidence_2015
endpoint$mortality_ratio_2030_vs_2015 <- endpoint$mortality_2030 / endpoint$mortality_2015
endpoint$incidence_reduction <- 1 - endpoint$incidence_ratio_2030_vs_2015
endpoint$mortality_reduction <- 1 - endpoint$mortality_ratio_2030_vs_2015
endpoint$who_target_attainment <- pmin(endpoint$incidence_reduction / .80, endpoint$mortality_reduction / .65)
endpoint$who_target_pass <- endpoint$incidence_ratio_2030_vs_2015 <= .20 & endpoint$mortality_ratio_2030_vs_2015 <= .35
endpoint$treatment_start_year <- 2026
endpoint$interpretation <- "WHO impact target: incidence reduction >=80% and HCV-related mortality reduction >=65% versus 2015; prevalence is secondary."

baseline <- endpoint[endpoint$initiation_prob_J == 0 & endpoint$initiation_prob_FX == 0, ][1, ]
endpoint$cumulative_viremic_person_years_averted <- baseline$cumulative_viremic_person_years - endpoint$cumulative_viremic_person_years
endpoint$endpoint_viremic_reduction_pp <- 100 * (baseline$total_viremic_prev - endpoint$total_viremic_prev)
endpoint$treatment_starts_per_1000_viremic_person_years_averted <- ifelse(endpoint$cumulative_viremic_person_years_averted > 0,
  1000 * endpoint$treatment_starts_total / endpoint$cumulative_viremic_person_years_averted, NA_real_)
tradeoff <- endpoint
tradeoff$pareto_efficient <- vapply(seq_len(nrow(tradeoff)), function(i) {
  other <- tradeoff[-i, , drop = FALSE]
  !any(other$treatment_starts_total <= tradeoff$treatment_starts_total[i] &
    other$cumulative_viremic_person_years <= tradeoff$cumulative_viremic_person_years[i] &
    (other$treatment_starts_total < tradeoff$treatment_starts_total[i] |
      other$cumulative_viremic_person_years < tradeoff$cumulative_viremic_person_years[i]))
}, logical(1))

write.csv(scenario_table, file.path(OUT, "scenario_definitions.csv"), row.names = FALSE)
writeLines(calibration_source, file.path(OUT, "calibration_source.txt"))
write.csv(trajectories, file.path(OUT, "scenario_trajectories.csv"), row.names = FALSE)
write.csv(endpoint, file.path(OUT, "scenario_2030_summary.csv"), row.names = FALSE)
write.csv(treatment, file.path(OUT, "scenario_treatment_summary.csv"), row.names = FALSE)
write.csv(tradeoff, file.path(OUT, "scenario_tradeoff.csv"), row.names = FALSE)
write.csv(subset(tradeoff, who_target_pass), file.path(OUT, "target_scenario_tradeoff.csv"), row.names = FALSE)
if (requireNamespace("writexl", quietly = TRUE))
  writexl::write_xlsx(list(Definitions = scenario_table, Endpoint = endpoint,
    Treatment = treatment, Tradeoff = tradeoff), file.path(OUT, "scenario_results.xlsx"))

library(ggplot2); library(patchwork); library(tidyr)
prob_labels <- paste0(fmt_prob(probabilities), "%")
shape_values <- c(`20_50` = 21, `30_60` = 24, all_ages = 22)
age_files <- setNames(age_settings$age_group_setting, age_settings$age_group_setting)

# Nature-style figures are separate for each age setting. Every 2015 point is
# exactly 100% because its denominator is its own 2015 reference rate.
for (ag in names(age_files)) {
  tr <- subset(trajectories, age_group_setting == ag)
  tr$Incidence <- tr$hcv_incidence_rate / tr$incidence_2015; tr$Mortality <- tr$hcv_mortality_rate / tr$mortality_2015
  long <- pivot_longer(tr, cols = c(Incidence, Mortality), names_to = "metric", values_to = "ratio")
  targets <- data.frame(metric = c("Incidence", "Mortality"), target = c(.20, .35))
  p1 <- ggplot(long, aes(year, ratio, colour = scenario)) +
    geom_line(linewidth = .7) + geom_point(data = subset(long, year == 2015), size = 1.4) +
    geom_hline(yintercept = 1, linetype = 3, colour = "grey35") +
    geom_hline(data = targets, aes(yintercept = target), linetype = 2, colour = "black") +
    scale_colour_viridis_d(option = "D", end = .95, guide = guide_legend(nrow = 5, byrow = TRUE)) +
    facet_wrap(~ metric) + scale_x_continuous(breaks = seq(2015, 2030, 5)) + scale_y_continuous(labels = scales::percent) +
    labs(x = NULL, y = "Rate relative to 2015 baseline", colour = "Scenarios", title = "A  WHO impact indicators from 2026") +
    theme_classic(base_size = 10)
  ep <- subset(endpoint, age_group_setting == ag)
  ep_long <- pivot_longer(ep, cols = c(incidence_reduction, mortality_reduction), names_to = "metric", values_to = "reduction")
  ep_long$metric <- ifelse(ep_long$metric == "incidence_reduction", "Incidence", "Mortality")
  ep_long$bar_colour <- ifelse(ep_long$who_target_pass, "#009E73", "#D55E00")
  p2 <- ggplot(ep_long, aes(reorder(scenario, reduction), reduction)) +
    geom_col(aes(fill = bar_colour), width = .72, show.legend = FALSE) + scale_fill_identity() +
    geom_hline(data = data.frame(metric = c("Incidence", "Mortality"), target = c(.80, .65)), aes(yintercept = target), linetype = 2) +
    coord_flip() + facet_wrap(~ metric) +
    scale_y_continuous(labels = scales::percent) + labs(x = NULL, y = "Reduction versus 2015", title = "B  2030 WHO target attainment") +
    theme_classic(base_size = 9)
  combined <- p1 + p2 + plot_layout(guides = "collect") & theme(legend.position = "bottom",
    legend.title.position = "left", legend.justification = "center")
  ggsave(file.path(OUT, paste0("scenario_2025_nature_style_", ag, ".png")), combined, width = 13, height = 8, dpi = 300)
}

# Treatment demand: integrated expected starts in J/F/X during 2026--2030.
# These are courses, not unique people; D is absent because it is prohibited.
pass <- subset(endpoint, who_target_pass)
pass_treatment <- pass[, c("scenario_id", "scenario", "scenario_type", "age_group_setting", "age_label", "total_viremic_prev",
  "incidence_reduction", "mortality_reduction", "initiation_prob_J", "initiation_prob_FX", "description", "plain_language")]
pass_treatment <- dplyr::left_join(pass_treatment, treatment,
  by = c("scenario_id", "scenario", "scenario_type"))
pass_treatment <- pass_treatment[order(pass_treatment$age_group_setting, pass_treatment$treatment_starts_total), ]
write.csv(pass_treatment, file.path(OUT, "target_scenario_treatment_summary.csv"), row.names = FALSE)
long_starts <- pivot_longer(pass_treatment, cols = c(treatment_starts_J, treatment_starts_F, treatment_starts_X), names_to = "stratum", values_to = "treatment_starts")
long_starts$stratum <- sub("treatment_starts_", "", long_starts$stratum)
for (ag in unique(pass_treatment$age_group_setting)) {
  x <- subset(long_starts, age_group_setting == ag); x$scenario <- factor(x$scenario, levels = rev(unique(x$scenario)))
  totals <- unique(x[c("scenario", "treatment_starts_total")])
  totals$label <- format(round(totals$treatment_starts_total), big.mark = ",", trim = TRUE)
  totals$scenario <- factor(totals$scenario, levels = levels(x$scenario))
  p3 <- ggplot(x, aes(scenario, treatment_starts, fill = stratum)) + geom_col() + coord_flip() +
    geom_text(data = totals, aes(x = scenario, y = treatment_starts_total, label = label),
      inherit.aes = FALSE, hjust = -.15, size = 3) +
    scale_y_continuous(expand = expansion(mult = c(0, .14))) +
    labs(x = NULL, y = "Treatment starts, 2026-2030", fill = "Stratum", title = paste0("A  Treatment demand: ", unique(x$age_label))) + theme_classic(base_size = 10)
  if (nrow(x) > 0) p3 <- p3 + scale_fill_manual(values = c(J = "#0072B2", F = "#D55E00", X = "#009E73"))
  q <- subset(pass_treatment, age_group_setting == ag)
  p4 <- ggplot(q, aes(reorder(scenario, treatment_courses_per_2026_person_total), treatment_courses_per_2026_person_total)) + geom_col(fill = "#6A3D9A") +
    scale_y_continuous(expand = expansion(mult = c(0, .08))) + coord_flip(ylim = c(.2, NA)) +
    labs(x = NULL, y = "Starts / population present in 2026", title = "B  Descriptive treatment intensity") + theme_classic(base_size = 10)
  ggsave(file.path(OUT, paste0("scenario_treatment_demand_", ag, ".png")), p3 + p4, width = 13, height = 8, dpi = 300)
}

# Heatmaps use discrete axes so all cells touch. D is absent because treatment
# in D is prohibited. The three age panels share one fixed gradient and legend.
heatmaps <- lapply(names(age_files), function(ag) {
  x <- subset(endpoint, age_group_setting == ag)
  x$J_label <- factor(paste0(fmt_prob(x$initiation_prob_J), "%"), levels = prob_labels)
  x$FX_label <- factor(paste0(fmt_prob(x$initiation_prob_FX), "%"), levels = prob_labels)
  p <- ggplot(x, aes(J_label, FX_label, fill = who_target_attainment)) + geom_tile(colour = "white", linewidth = .3, width = 1, height = 1) +
    geom_text(aes(label = ifelse(who_target_pass, "PASS", "")), size = 3) + scale_x_discrete(drop = FALSE) + scale_y_discrete(drop = FALSE) +
    scale_fill_viridis_c(limits = c(0, 1), breaks = seq(0, 1, .25), labels = scales::percent,
      oob = scales::squish) + labs(x = "J treatment initiation probability", y = "F/X treatment initiation probability",
      fill = "Minimum WHO\ntarget attainment", title = unique(x$age_label)) +
    theme_classic(base_size = 10) + theme(plot.title = element_text(hjust = .5))
  p
})
heatmap_combined <- patchwork::wrap_plots(heatmaps, nrow = 1, guides = "collect") &
  theme(legend.position = "right")
ggsave(file.path(OUT, "scenario_initiation_probability_heatmap.png"), heatmap_combined,
  width = 18, height = 7, dpi = 300)

# One combined trade-off plot contains all policies. Dominated policies are
# hollow, non-dominated policies are filled by the viremic-prevalence scale,
# and the three duplicate no-treatment runs are represented by one status quo.
status_quo <- tradeoff[tradeoff$initiation_prob_J == 0 & tradeoff$initiation_prob_FX == 0, ][1, ]
tradeoff_plot <- subset(tradeoff, !(initiation_prob_J == 0 & initiation_prob_FX == 0))
p6 <- ggplot() +
  geom_point(data = subset(tradeoff_plot, !pareto_efficient),
    aes(treatment_starts_total, cumulative_viremic_person_years_averted,
      colour = total_viremic_prev, shape = age_group_setting),
    fill = NA, size = 3.1, alpha = .9, stroke = 1) +
  geom_point(data = subset(tradeoff_plot, pareto_efficient),
    aes(treatment_starts_total, cumulative_viremic_person_years_averted,
      colour = total_viremic_prev, fill = total_viremic_prev, shape = age_group_setting),
    size = 3.1, alpha = .9, stroke = 1) +
  geom_point(data = status_quo, aes(treatment_starts_total, cumulative_viremic_person_years_averted),
    shape = 16, colour = "black", fill = "black", size = 3.2, inherit.aes = FALSE) +
  geom_text(data = subset(tradeoff_plot, who_target_pass),
    aes(x = treatment_starts_total, y = cumulative_viremic_person_years_averted, label = scenario), colour = "black",
    hjust = -.05, vjust = -.35, size = 2.2, show.legend = FALSE) +
  scale_shape_manual(values = shape_values, labels = c(`20_50` = "20-50", `30_60` = "30-60", all_ages = "All age groups"), name = "Age setting") +
  scale_colour_viridis_c(trans = "sqrt", labels = scales::percent, name = "2030\nviremic prevalence") +
  scale_fill_viridis_c(trans = "sqrt", labels = scales::percent, name = "2030\nviremic prevalence",
    guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(.02, .18))) + labs(x = "Treatment starts, 2026-2030",
    y = "Cumulative viremic person-years averted vs status quo", title = "Scenario trade-off across all age settings") + coord_cartesian(clip = "off") +
  theme_classic(base_size = 10) + theme(plot.margin = margin(5.5, 70, 5.5, 5.5))
ggsave(file.path(OUT, "scenario_tradeoff.png"), p6, width = 12, height = 7, dpi = 300)
cat("Wrote", nrow(scenario_table), "scenarios to", OUT, "\n")
