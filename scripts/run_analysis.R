# =============================================================================
# run_analysis.R — 50-year sensitivity projection with posterior CrIs
#
# Usage:
#   Rscript scripts/run_analysis.R \
#     --root . \
#     --fit output/calibration/<run>/fit.rds \
#     --posterior output/calibration/npe_bayes/posterior_samples_mcmc.csv \
#     --out-dir output/analysis --n-draws 300 --n-cores 4
# =============================================================================

suppressMessages({
  library(Rcpp)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(name, default = NULL) {
  pos <- match(name, args)
  if (is.na(pos)) return(default)
  args[pos + 1]
}

ROOT <- normalizePath(arg_val("--root", getwd()), mustWork = TRUE)
FIT_PATH <- normalizePath(arg_val("--fit", file.path(ROOT, "output", "calibration", "run1_4strata", "fit.rds")), mustWork = TRUE)
POST_PATH <- normalizePath(arg_val("--posterior", file.path(ROOT, "output", "calibration", "npe_bayes", "posterior_samples_mcmc.csv")), mustWork = TRUE)
OUT_DIR <- normalizePath(arg_val("--out-dir", file.path(ROOT, "output", "analysis")), mustWork = FALSE)
N_DRAWS <- as.integer(arg_val("--n-draws", "300"))
N_CORES <- as.integer(arg_val("--n-cores", "4"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

old_wd <- setwd(file.path(ROOT, "src"))
setup_env <- new.env(parent = globalenv()); sys.source("setup.R", envir = setup_env); setwd(old_wd)
if (!exists("run_sim", envir = globalenv(), inherits = FALSE)) Rcpp::sourceCpp(file.path(ROOT, "src", "sim.cpp"))
run_sim <- get("run_sim", envir = globalenv(), inherits = FALSE)
for (f in c("targets.R", "model_metrics.R", "equilibrium.R", "likelihood.R", "calibrate_nm.R")) {
  sys.source(file.path(ROOT, "src", "calibration", f), envir = environment())
}
fit <- readRDS(FIT_PATH)
base_params <- setup_env$params
post <- read.csv(POST_PATH)
th_all <- as.matrix(post[, paste0("theta", 1:12)])
set.seed(2026)
idx <- sample(seq_len(nrow(th_all)), min(N_DRAWS, nrow(th_all)))
th_draws <- th_all[idx, , drop = FALSE]

# Equilibrium horizon from the calibration run
eq_data <- setup_env$data
eq_data$t_start <- fit$t_start
eq_data$t_end <- fit$t_end

scenario_path <- file.path(ROOT, "src", "calibration", "scenarios.csv")
scenario_df <- read.csv(scenario_path, stringsAsFactors = FALSE)

load_scenarios <- function(sc_df, base_params) {
  lapply(seq_len(nrow(sc_df)), function(j) {
    sc <- sc_df[j, ]
    tau_stratum <- c(sc$elig_D, sc$elig_J, sc$elig_F, sc$elig_X)
    phi <- if (sc$phi_mode == "baseline") c(0.07, 0.23, 1) else c(1, 1, 1)
    alpha_DC <- if (sc$alpha_DC_mode == "pos") base_params$alpha_DC_pos else
      base_params$alpha_DC_neg
    list(
      id = sc$scenario,
      description = sc$description,
      tau = as.numeric(sc[c("tau_NC", "tau_CC", "tau_DC", "tau_HCC")]),
      rho = sc$rho,
      alpha_DC = alpha_DC,
      phi = phi,
      tau_stratum = as.numeric(tau_stratum),
      tau_min_age = as.integer(sc$min_age_group),
      figure_family = NA_character_,
      coverage = NA_real_,
      coverage_label = NA_character_
    )
  })
}

scenarios <- load_scenarios(scenario_df, base_params)

### F2 & F3 scenarios
# Additional scenario families used for figures analogous to F2/F3.  For D,
# F, and X we vary the annual treatment-initiation rate from 1% to 15% of
# chronic HCV cases in the targeted stratum.  J has four 100%-uptake,
# infection-stage-specific scenarios.  Combined varies community treatment
# from 1% to 15%/y while prison treatment is fixed at 100%/y.  The panels report total Singapore outcomes for each targeted
# stratum and do not apply interactions between targeted strata.
figure_rate_sets <- list(
  D = seq(0.01, 0.15, by = 0.01),
  # J = rep(1, 4),
  J = seq(0.1, 1.0, by = 0.1),
  F = seq(0.01, 0.15, by = 0.01),
  X = seq(0.01, 0.15, by = 0.01),
  combined = seq(0.01, 0.15, by = 0.01)
)

stratum_labels <- c(
  D = "D: never-arrested active PWID",
  J = "J: currently detained PWID",
  F = "F: ever-arrested active PWID",
  X = "X: former non-injecting PWID",
  combined = "Combined: D + J + F + X"
)

# Short labels used in the multi-panel figures.
plot_stratum_labels <- c(D = "D", J = "J", F = "F", X = "X",
                         combined = "Total")

family_palette <- list(
  D = c("#d9f0a3", "#006837"),
  J = c("#e0ecf4", "#8856a7"),
  F = c("#c6dbef", "#084594"),
  X = c("#fee8c8", "#e34a33"),
  combined = c("#fddbc7", "#b2182b")
)

make_figure_family <- function(prefix, description, elig, rates,
                               labels = paste0(round(100 * rates), "%/y"),
                               tau_builder = function(rate, i) rep(rate, 4),
                               tau_by_stratum_builder = NULL) {
  colours <- grDevices::colorRampPalette(family_palette[[prefix]])(length(rates))
  lapply(seq_along(rates), function(i) {
    rate <- rates[i]
    list(
      id = paste0("figure_", prefix, "_", sprintf("%02d", i)),
      description = paste(description, labels[i]),
      tau = tau_builder(rate, i),
      tau_by_stratum = if (is.null(tau_by_stratum_builder)) NULL else
        as.vector(t(tau_by_stratum_builder(rate, i))),
      rho = base_params$rho,
      alpha_DC = base_params$alpha_DC_pos,
      phi = c(base_params$phi_CC_DC, base_params$phi_CC_HCC,
              base_params$phi_DC_HCC),
      tau_stratum = as.numeric(elig),
      tau_min_age = 1L,
      figure_family = prefix,
      coverage = rate,
      coverage_label = labels[i],
      figure_colour = colours[i]
    )
  })
}

make_j_family <- function() {
  stages <- c("NC", "CC", "DC", "HCC")
  colours <- grDevices::colorRampPalette(family_palette$J)(length(stages))
  lapply(seq_along(stages), function(i) {
    tau <- rep(0, 4); tau[i] <- 1
    list(
      id = paste0("figure_J_", tolower(stages[i])),
      description = paste("J: 100% uptake, treatment at", stages[i]),
      tau = tau,
      rho = base_params$rho,
      alpha_DC = base_params$alpha_DC_pos,
      phi = c(base_params$phi_CC_DC, base_params$phi_CC_HCC,
              base_params$phi_DC_HCC),
      tau_stratum = c(0, 1, 0, 0),
      tau_min_age = 1L,
      figure_family = "J",
      coverage = 1,
      coverage_label = paste0(stages[i], " (100% uptake)"),
      figure_colour = colours[i]
    )
  })
}

figure_baseline <- list(
  id = "figure_status_quo",
  description = "Status quo: no treatment",
  tau = c(0, 0, 0, 0),
  rho = base_params$rho,
  alpha_DC = base_params$alpha_DC_pos,
  phi = c(base_params$phi_CC_DC, base_params$phi_CC_HCC,
          base_params$phi_DC_HCC),
  tau_stratum = c(0, 0, 0, 0),
  tau_min_age = 1L,
  figure_family = "status_quo",
  coverage = 0,
  coverage_label = "Status quo",
  figure_colour = "#222222"
)

figure_scenarios <- c(
  list(figure_baseline),
  make_figure_family("D", "D treatment", c(1, 0, 0, 0), figure_rate_sets$D),
  # make_j_family(),
  make_figure_family("J", "J treatment", c(0, 1, 0, 0), figure_rate_sets$J),
  make_figure_family("F", "F treatment", c(0, 0, 1, 0), figure_rate_sets$F),
  make_figure_family("X", "X treatment", c(0, 0, 0, 1), figure_rate_sets$X),
  make_figure_family("combined", "Combined treatment", c(1, 1, 1, 1),
                     figure_rate_sets$combined,
                     tau_by_stratum_builder = function(rate, i) {
                       rbind(D = rep(rate, 4), J = rep(1, 4),
                             F = rep(rate, 4), X = rep(rate, 4))
                     })
)

analysis_scenarios <- c(scenarios, figure_scenarios)

# Summary indices
idx_hcv <- function(strata = 0:3, stages = 1:4, states = 1:3, ages = 0:5) {
  as.vector(sapply(strata, function(s) sapply(stages, function(k)
    sapply(states, function(h) sapply(ages, function(i)
      s * 96 + (k - 1) * 24 + h * 6 + i + 2L)))))
}
idx_dc <- function(strata = 0:3, states = 0:3, ages = 0:5) {
  as.vector(sapply(strata, function(s) sapply(states, function(h)
    sapply(ages, function(i) s * 96 + (3 - 1) * 24 + h * 6 + i + 2L))))
}
idx_hcc <- function(strata = 0:3, states = 0:3, ages = 0:5) {
  as.vector(sapply(strata, function(s) sapply(states, function(h)
    sapply(ages, function(i) s * 96 + (4 - 1) * 24 + h * 6 + i + 2L))))
}

comp_cols <- function(s, k, h) {
  s * 96 + (k - 1) * 24 + h * 6 + (0:5) + 2L
}

state_cols <- function(s, k, h) comp_cols(s, k, h) - 1L

sum_comp <- function(out, strata, stages = 1:4, states = 0:3) {
  cols <- unlist(lapply(strata, function(s)
    unlist(lapply(stages, function(k)
      unlist(lapply(states, function(h) comp_cols(s, k, h)))))))
  rowSums(out[, cols, drop = FALSE])
}

stratum_snapshot <- function(state, strata) {
  cols_total <- unlist(lapply(strata, function(s)
    unlist(lapply(1:4, function(k)
      unlist(lapply(0:3, function(h) state_cols(s, k, h)))))))
  cols_hcv <- unlist(lapply(strata, function(s)
    unlist(lapply(1:4, function(k)
      unlist(lapply(1:3, function(h) state_cols(s, k, h)))))))
  c(total = sum(state[cols_total]), hcv = sum(state[cols_hcv]))
}

equilibrium_one <- function(theta) {
  pm <- tryCatch(build_params(theta, base_params), error = function(e) NULL)
  if (is.null(pm)) return(NULL)
  out_eq <- tryCatch(run_sim(pm, eq_data), error = function(e) NULL)
  if (is.null(out_eq)) return(NULL)
  list(pm = pm, y0 = out_eq[nrow(out_eq), 2:385])
}

# Calculate the flows represented by the ODE.  Cumulative DC/HCC include the
# advanced-stage stock present in 2017 plus incident transitions thereafter.
# Treatment and death counters start at zero in 2017.
cumulative_flows <- function(out, pm, scenario) {
  n <- nrow(out)
  flow_dc <- matrix(0, nrow = n, ncol = 4)
  flow_hcc <- matrix(0, nrow = n, ncol = 4)
  flow_treatment <- matrix(0, nrow = n, ncol = 4)
  flow_death <- matrix(0, nrow = n, ncol = 4)

  ptc_CC_DC <- (pm$rho * pm$r3_CC_DC + (1 - pm$rho)) * pm$p_CC_DC
  ptc_CC_HCC <- (pm$rho * pm$r3_CC_HCC + (1 - pm$rho)) * pm$p_CC_HCC
  ptc_DC_HCC <- (pm$rho * pm$r3_DC_HCC + (1 - pm$rho)) * pm$p_DC_HCC

  for (s in 0:3) {
    for (i in 0:5) {
      cc <- out[, comp_cols(s, 2, 2)]
      dc <- out[, comp_cols(s, 3, 2)]
      flow_dc[, s + 1] <- flow_dc[, s + 1] + ptc_CC_DC * cc[, i + 1]
      flow_hcc[, s + 1] <- flow_hcc[, s + 1] +
        ptc_CC_HCC * cc[, i + 1] + ptc_DC_HCC * dc[, i + 1]

      eligible <- scenario$tau_stratum[s + 1] *
        as.numeric((i + 1) >= scenario$tau_min_age)
      for (k in 1:4) {
        chronic <- out[, comp_cols(s, k, 2)[i + 1]]
        treatment_rate <- if (is.null(scenario$tau_by_stratum)) {
          scenario$tau[k] * eligible
        } else {
          scenario$tau_by_stratum[s * 4 + k] *
            as.numeric((i + 1) >= scenario$tau_min_age)
        }
        flow_treatment[, s + 1] <- flow_treatment[, s + 1] +
          treatment_rate * chronic
      }

      base_mu <- pm$mu[i + 1] * pm$omega
      # Death flow is calculated separately for each HCV state because t has
      # the post-SVR mortality modifier.
      for (h in 0:3) {
        dc_rate <- if (h == 3) pm$psi_DC * (base_mu + pm$mu_DC) else
          base_mu + pm$mu_DC
        hcc_rate <- if (h == 3) pm$psi_HCC * (base_mu + pm$mu_HCC) else
          base_mu + pm$mu_HCC
        flow_death[, s + 1] <- flow_death[, s + 1] +
          dc_rate * out[, comp_cols(s, 3, h)[i + 1]] +
          hcc_rate * out[, comp_cols(s, 4, h)[i + 1]]
      }
    }
  }

  integrate_flow <- function(flow) {
    result <- matrix(0, nrow = n, ncol = ncol(flow))
    if (n > 1) {
      increments <- flow[-n, , drop = FALSE] * diff(out[, 1])
      result[-1, ] <- apply(increments, 2, cumsum)
    }
    result
  }

  dc_cum <- integrate_flow(flow_dc)
  hcc_cum <- integrate_flow(flow_hcc)
  treatment_cum <- integrate_flow(flow_treatment)
  death_cum <- integrate_flow(flow_death)

  initial_dc <- sapply(0:3, function(s) sum(out[1, idx_dc(strata = s)]))
  initial_hcc <- sapply(0:3, function(s) sum(out[1, idx_hcc(strata = s)]))
  dc_cum <- sweep(dc_cum, 2, initial_dc, "+")
  hcc_cum <- sweep(hcc_cum, 2, initial_hcc, "+")

  list(
    dc = dc_cum,
    hcc = hcc_cum,
    treatment = treatment_cum,
    death = death_cum
  )
}

project_one <- function(eq_state, scenario) {
  if (is.null(eq_state)) return(NULL)
  pm <- eq_state$pm
  y0 <- eq_state$y0
  pm$tau <- scenario$tau
  pm$rho <- scenario$rho
  pm$alpha_DC_pos <- scenario$alpha_DC
  pm$phi_CC_DC <- scenario$phi[1]
  pm$phi_CC_HCC <- scenario$phi[2]
  pm$phi_DC_HCC <- scenario$phi[3]
  pm$tau_stratum <- scenario$tau_stratum
  if (!is.null(scenario$tau_by_stratum)) {
    pm$tau_by_stratum <- scenario$tau_by_stratum
  }
  pm$tau_min_age <- scenario$tau_min_age
  proj_data <- eq_data
  proj_data$t_start <- 47
  proj_data$t_end <- 97
  proj_data$y0 <- y0
  out <- tryCatch(run_sim(pm, proj_data), error = function(e) NULL)
  if (is.null(out)) return(NULL)
  years <- 0:50
  rows <- vapply(years, function(yr) {
    which.min(abs(out[, 1] - (47 + yr)))
  }, integer(1))
  groups <- list(D = 0, J = 1, F = 2, X = 3, combined = 0:3)
  flows <- cumulative_flows(out, pm, scenario)
  result <- data.frame(year = 2017 + years)
  for (group_name in names(groups)) {
    strata <- groups[[group_name]]
    result[[paste0("hcv_", group_name)]] <- sum_comp(out[rows, , drop = FALSE],
                                                       strata, states = 1:3)
    result[[paste0("dc_cumulative_", group_name)]] <-
      rowSums(flows$dc[rows, strata + 1, drop = FALSE])
    result[[paste0("hcc_cumulative_", group_name)]] <-
      rowSums(flows$hcc[rows, strata + 1, drop = FALSE])
    result[[paste0("treatments_", group_name)]] <-
      rowSums(flows$treatment[rows, strata + 1, drop = FALSE])
    result[[paste0("deaths_", group_name)]] <-
      rowSums(flows$death[rows, strata + 1, drop = FALSE])
  }
  result
}

cat("Running", length(analysis_scenarios), "scenarios x", nrow(th_draws), "posterior draws...\n")
cat("Equilibrating posterior draws...\n")
eq_states <- if (N_CORES > 1 && requireNamespace("parallel", quietly = TRUE)) {
  parallel::mclapply(seq_len(nrow(th_draws)), function(j) {
    equilibrium_one(th_draws[j, ])
  }, mc.cores = N_CORES)
} else {
  lapply(seq_len(nrow(th_draws)), function(j) equilibrium_one(th_draws[j, ]))
}

initial_dat <- bind_rows(lapply(seq_along(eq_states), function(j) {
  eq <- eq_states[[j]]
  if (is.null(eq)) return(NULL)
  strata <- list(D = 0, J = 1, F = 2, X = 3, total = 0:3)
  bind_rows(lapply(names(strata), function(nm) {
    z <- stratum_snapshot(eq$y0, strata[[nm]])
    data.frame(draw = j, stratum = nm, total = z["total"], hcv = z["hcv"])
  }))
}))

initial_bci <- initial_dat %>%
  pivot_longer(c(total, hcv), names_to = "measure", values_to = "value") %>%
  group_by(stratum, measure) %>%
  summarise(
    median = median(value),
    BCI_2.5 = as.numeric(quantile(value, 0.025, names = FALSE)),
    BCI_97.5 = as.numeric(quantile(value, 0.975, names = FALSE)),
    .groups = "drop"
  ) %>%
  mutate(stratum = factor(stratum, levels = c("D", "J", "F", "X", "total"))) %>%
  arrange(stratum, measure) %>%
  mutate(stratum = as.character(stratum))
write.csv(initial_bci, file.path(OUT_DIR, "posterior_bci_strata_2017_long.csv"),
          row.names = FALSE)
initial_bci_wide <- initial_bci %>%
  pivot_wider(names_from = measure,
              values_from = c(median, BCI_2.5, BCI_97.5),
              names_glue = "{measure}_{.value}") %>%
  select(stratum, total_median, total_BCI_2.5, total_BCI_97.5,
         hcv_median, hcv_BCI_2.5, hcv_BCI_97.5)
write.csv(initial_bci_wide, file.path(OUT_DIR, "posterior_bci_strata_2017.csv"),
          row.names = FALSE)

all_res <- list()
for (sc in analysis_scenarios) {
  run_one <- function(j) {
    d <- project_one(eq_states[[j]], sc)
    if (is.null(d)) return(NULL)
    d$draw <- j
    d
  }
  res <- if (N_CORES > 1 && requireNamespace("parallel", quietly = TRUE)) {
    parallel::mclapply(seq_len(nrow(th_draws)), run_one, mc.cores = N_CORES)
  } else {
    lapply(seq_len(nrow(th_draws)), run_one)
  }
  ok <- !sapply(res, is.null)
  d <- do.call(rbind, res[ok])
  d$scenario <- sc$id
  all_res[[sc$id]] <- d
  cat("  scenario", sc$id, "done (", sum(ok), "draws )\n")
}

dat <- do.call(rbind, all_res)
measure_cols <- setdiff(names(dat), c("scenario", "draw", "year"))
summary <- dat %>%
  group_by(scenario, year) %>%
  summarise(across(all_of(measure_cols), list(
    median = median,
    lo = ~as.numeric(quantile(.x, 0.025, names = FALSE)),
    hi = ~as.numeric(quantile(.x, 0.975, names = FALSE))
  )), .groups = "drop")
write.csv(summary, file.path(OUT_DIR, "scenario_summary.csv"), row.names = FALSE)

key_years <- c(2017, 2027, 2037, 2047, 2057, 2067)
key_table <- summary %>% filter(year %in% key_years)
write.csv(key_table, file.path(OUT_DIR, "scenario_key_years.csv"), row.names = FALSE)

theme_pub <- function() {
  theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA),
          strip.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour = "grey60", linewidth = 0.4),
          legend.position = "right")
}

scenario_meta <- bind_rows(lapply(seq_along(analysis_scenarios), function(ii) {
  sc <- analysis_scenarios[[ii]]
  data.frame(scenario = sc$id, description = sc$description,
             figure_family = sc$figure_family,
             coverage = sc$coverage,
             coverage_label = sc$coverage_label,
             figure_colour = if (is.null(sc$figure_colour)) NA_character_ else
               sc$figure_colour,
             stringsAsFactors = FALSE)
}))

figure_targets <- c("D", "J", "F", "X", "combined")
figure_meta <- scenario_meta %>%
  filter(!is.na(figure_family)) %>%
  mutate(figure_colour = ifelse(is.na(figure_colour), "#222222", figure_colour))
write.csv(figure_meta %>%
            select(scenario, description, figure_family, coverage,
                   coverage_label, figure_colour),
          file.path(OUT_DIR, "figure_scenario_definitions.csv"), row.names = FALSE)

make_figure_long <- function(metric_prefix, outcome_labels,
                             total_only = FALSE) {
  metric_names <- if (total_only) "combined" else
    c("D", "J", "F", "X", "combined")
  metric_cols <- paste0(metric_prefix, metric_names)
  raw <- dat %>%
    filter(scenario %in% figure_meta$scenario) %>%
    select(scenario, draw, year, all_of(metric_cols)) %>%
    pivot_longer(all_of(metric_cols), names_to = "outcome_key", values_to = "value") %>%
    mutate(outcome = if (total_only) {
      factor(outcome_labels, levels = outcome_labels)
    } else {
      factor(plot_stratum_labels[sub("^.*_", "", outcome_key)],
             levels = unname(plot_stratum_labels))
    }) %>%
    left_join(figure_meta, by = "scenario") %>%
    mutate(panel = ifelse(scenario == "figure_status_quo", NA_character_,
                          figure_family))

  baseline <- raw %>%
    filter(scenario == "figure_status_quo") %>%
    select(-panel) %>%
    tidyr::crossing(panel = figure_targets)
  interventions <- raw %>% filter(scenario != "figure_status_quo")

  bind_rows(baseline, interventions) %>%
    mutate(panel = factor(panel, levels = figure_targets),
           panel = recode(panel,
                          D = plot_stratum_labels[["D"]],
                          J = plot_stratum_labels[["J"]],
                          F = plot_stratum_labels[["F"]],
                          X = plot_stratum_labels[["X"]],
                          combined = plot_stratum_labels[["combined"]]),
           panel = factor(panel, levels = unname(plot_stratum_labels))) %>%
    group_by(scenario, description, figure_family, coverage_label, figure_colour,
             panel, outcome, year) %>%
    summarise(
      median = median(value),
      lo = as.numeric(quantile(value, 0.025, names = FALSE)),
      hi = as.numeric(quantile(value, 0.975, names = FALSE)),
      .groups = "drop"
    )
}

plot_f2 <- function() {
  f2 <- make_figure_long("hcv_", plot_stratum_labels)
  f2_labels <- expand.grid(panel = factor(unname(plot_stratum_labels),
                                          levels = unname(plot_stratum_labels)),
                           outcome = factor(unname(plot_stratum_labels),
                                             levels = unname(plot_stratum_labels))) %>%
    mutate(label = letters[seq_len(n())])
  ggplot(f2, aes(x = year - 2017, y = median / 1000, group = scenario,
                 colour = figure_colour, fill = figure_colour)) +
    geom_line(linewidth = 0.7) +
    geom_text(data = f2_labels,
              aes(x = 0.5, y = Inf, label = paste0("(", label, ")")),
              inherit.aes = FALSE, hjust = 0, vjust = 1.3, size = 5.5) +
    facet_grid(outcome ~ panel, scales = "free_y", switch = "y") +
    scale_colour_identity(guide = "none") +
    scale_fill_identity(guide = "none") +
    scale_x_continuous(breaks = seq(0, 50, 10)) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 5),
                       expand = expansion(mult = c(0.02, 0.10))) +
    labs(x = "Time (years)", y = "Infections (000s)") +
    theme_pub() +
    theme(text = element_text(size = 19), axis.title = element_text(size = 19),
          axis.text = element_text(size = 16), strip.placement = "outside",
          panel.spacing = grid::unit(1.2, "lines"),
          strip.background = element_blank(), strip.text.y.left = element_text(
            angle = 90, size = 20), strip.text.x = element_text(size = 18),
          legend.position = "none")
}

plot_f3 <- function() {
  f3_parts <- list(
    make_figure_long("dc_cumulative_", "Cumulative DC\ncases (000s)", TRUE),
    make_figure_long("hcc_cumulative_", "Cumulative HCC\ncases (000s)", TRUE),
    make_figure_long("treatments_", "Cumulative\ntreatments (000s)", TRUE),
    make_figure_long("deaths_", "Cumulative advanced-stage\ndeaths (000s)", TRUE)
  )
  f3 <- bind_rows(f3_parts) %>%
    mutate(outcome = factor(outcome, levels = c("Cumulative DC\ncases (000s)",
                                                "Cumulative HCC\ncases (000s)",
                                                "Cumulative\ntreatments (000s)",
                                                "Cumulative advanced-stage\ndeaths (000s)")))
  f3_labels <- expand.grid(panel = factor(unname(plot_stratum_labels),
                                          levels = unname(plot_stratum_labels)),
                           outcome = levels(f3$outcome)) %>%
    mutate(outcome = factor(outcome, levels = levels(f3$outcome)),
           label = letters[seq_len(n())])
  ggplot(f3, aes(x = year - 2017, y = median / 1000, group = scenario,
                 colour = figure_colour, fill = figure_colour)) +
    geom_line(linewidth = 0.7) +
    geom_text(data = f3_labels,
              aes(x = 0.5, y = Inf, label = paste0("(", label, ")")),
              inherit.aes = FALSE, hjust = 0, vjust = 1.3, size = 6.5) +
    facet_grid(outcome ~ panel, scales = "free_y", switch = "y") +
    scale_colour_identity(guide = "none") +
    scale_fill_identity(guide = "none") +
    scale_x_continuous(breaks = seq(0, 50, 10)) +
    scale_y_continuous(breaks = scales::breaks_pretty(n = 5),
                       expand = expansion(mult = c(0.02, 0.10))) +
    labs(x = "Time (years)", y = "Number of people (000s)") +
    theme_pub() +
    theme(text = element_text(size = 22), axis.title = element_text(size = 22),
          axis.text = element_text(size = 18), strip.placement = "outside",
          panel.spacing = grid::unit(1.2, "lines"),
          strip.background = element_blank(), strip.text.y.left = element_text(
            angle = 90, size = 22), strip.text.x = element_text(size = 20),
          legend.position = "none")
}

p_f2 <- plot_f2()
p_f3 <- plot_f3()
legend_title_size <- 6.0
legend_text_size <- 5.0
make_right_legend <- function() {
  min_labels <- c(D = "1%/y", J = "10%/y", F = "1%/y", X = "1%/y",
                  combined = "1%/y")
  max_labels <- c(D = "15%/y", J = "100%y", F = "15%/y", X = "15%/y",
                  combined = "15%/y")
  families <- figure_targets
  lines <- bind_rows(lapply(seq_along(families), function(i) {
    family <- families[i]
    meta <- figure_meta %>% filter(figure_family == family)
    edges <- seq(1.7, 3.9, length.out = nrow(meta) + 1)
    y <- 1.0 + (length(families) - i) * 0.05
    data.frame(family = family, figure_colour = meta$figure_colour,
               x = edges[-length(edges)], xend = edges[-1], y = y, yend = y)
  }))
  labels <- data.frame(
    family = families, display_label = unname(plot_stratum_labels[families]),
    min_label = unname(min_labels[families]), max_label = unname(max_labels[families]),
    y = 1.0 + (length(families) - seq_along(families)) * 0.05,
    stringsAsFactors = FALSE
  )
  ggplot() +
    geom_text(data = data.frame(x = 0.05, y = max(labels$y) + 0.20,
                                label = "Treatment scenarios"),
              aes(x = x, y = y, label = label), hjust = 0,
              size = legend_title_size) +
    geom_segment(data = lines,
                 aes(x = x, xend = xend, y = y, yend = y,
                     colour = figure_colour), linewidth = 2.2, lineend = "butt") +
    geom_text(data = labels, aes(x = 0.05, y = y, label = display_label),
              hjust = 0, size = legend_text_size) +
    geom_text(data = labels, aes(x = 1.6, y = y, label = min_label),
              hjust = 1, size = legend_text_size) +
    geom_text(data = labels, aes(x = 4.0, y = y, label = max_label),
              hjust = 0, size = legend_text_size) +
    scale_colour_identity() +
    coord_cartesian(xlim = c(0, 4.6), ylim = c(0.2, 2.0), expand = FALSE,
                    clip = "off") +
    theme_void() +
    theme(plot.margin = margin(10, 10, 10, 35),
          text = element_text(size = 25))
}
p_legend <- make_right_legend()
p_f2_final <- p_f2 + p_legend + plot_layout(widths = c(1, 0.24))
p_f3_final <- p_f3 + p_legend + plot_layout(widths = c(1, 0.24))
FIG_DIR <- file.path(ROOT, "figures", "analysis")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(OUT_DIR, "fig_F2_hcv_scenarios.png"), p_f2_final,
       width = 18, height = 12, dpi = 300)
ggsave(file.path(OUT_DIR, "fig_F3_cumulative_advanced_hcv.png"), p_f3_final,
       width = 23, height = 14, dpi = 300)
file.copy(file.path(OUT_DIR, "fig_F2_hcv_scenarios.png"),
          file.path(FIG_DIR, "fig_F2_hcv_scenarios.png"), overwrite = TRUE)
file.copy(file.path(OUT_DIR, "fig_F3_cumulative_advanced_hcv.png"),
          file.path(FIG_DIR, "fig_F3_cumulative_advanced_hcv.png"), overwrite = TRUE)

scenario_order <- c(scenario_df$scenario,
                    vapply(figure_scenarios, function(x) x$id, character(1)))
dat$scenario <- factor(dat$scenario, levels = scenario_order)

cat("Analysis outputs written to", OUT_DIR, "\n")
