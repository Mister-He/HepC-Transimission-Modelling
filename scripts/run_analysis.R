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
      tau_min_age = as.integer(sc$min_age_group)
    )
  })
}

scenarios <- load_scenarios(scenario_df, base_params)

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

summarise_at_year <- function(out, years, year0 = 47) {
  times <- out[, 1]
  sapply(years, function(yr) {
    row <- which.min(abs(times - (year0 + yr)))
    c(hcv = sum(out[row, idx_hcv()]),
      dc = sum(out[row, idx_dc()]),
      hcc = sum(out[row, idx_hcc()]))
  })
}

project_one <- function(theta, scenario) {
  pm <- tryCatch(build_params(theta, base_params), error = function(e) NULL)
  if (is.null(pm)) return(NULL)
  out_eq <- tryCatch(run_sim(pm, eq_data), error = function(e) NULL)
  if (is.null(out_eq)) return(NULL)
  y0 <- out_eq[nrow(out_eq), 2:385]
  pm$tau <- scenario$tau
  pm$rho <- scenario$rho
  pm$alpha_DC_pos <- scenario$alpha_DC
  pm$phi_CC_DC <- scenario$phi[1]
  pm$phi_CC_HCC <- scenario$phi[2]
  pm$phi_DC_HCC <- scenario$phi[3]
  pm$tau_stratum <- scenario$tau_stratum
  pm$tau_min_age <- scenario$tau_min_age
  proj_data <- eq_data
  proj_data$t_start <- 47
  proj_data$t_end <- 97
  proj_data$y0 <- y0
  out <- tryCatch(run_sim(pm, proj_data), error = function(e) NULL)
  if (is.null(out)) return(NULL)
  years <- 0:50
  mat <- summarise_at_year(out, years)
  data.frame(year = 2017 + years,
             hcv = mat["hcv", ], dc = mat["dc", ], hcc = mat["hcc", ])
}

cat("Running", length(scenarios), "scenarios x", nrow(th_draws), "posterior draws...\n")
all_res <- list()
for (sc in scenarios) {
  run_one <- function(j) {
    d <- project_one(th_draws[j, ], sc)
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
summary <- dat %>%
  group_by(scenario, year) %>%
  summarise(
    hcv_median = median(hcv), hcv_lo = quantile(hcv, 0.025), hcv_hi = quantile(hcv, 0.975),
    dc_median = median(dc), dc_lo = quantile(dc, 0.025), dc_hi = quantile(dc, 0.975),
    hcc_median = median(hcc), hcc_lo = quantile(hcc, 0.025), hcc_hi = quantile(hcc, 0.975),
    .groups = "drop"
  )
write.csv(summary, file.path(OUT_DIR, "scenario_summary.csv"), row.names = FALSE)

key_years <- c(2017, 2027, 2037, 2047, 2057, 2067)
key_table <- summary %>% filter(year %in% key_years)
write.csv(key_table, file.path(OUT_DIR, "scenario_key_years.csv"), row.names = FALSE)

okabe <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
           "#EDC948", "#B07AA1", "#FF9DA7", "#9C755F", "#BAB0AC",
           "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
           "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF")
theme_pub <- function() {
  theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA),
          strip.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour = "grey60", linewidth = 0.4),
          legend.position = "right")
}

scenario_order <- scenario_df$scenario
dat$scenario <- factor(dat$scenario, levels = scenario_order)

p1 <- ggplot(summary, aes(x = year, y = hcv_median, colour = scenario)) +
  geom_ribbon(aes(ymin = hcv_lo, ymax = hcv_hi, fill = scenario), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = okabe) + scale_fill_manual(values = okabe) +
  labs(x = "Year", y = "Total HCV count",
       title = "Total HCV (acute + chronic + treatment) projection, 2017-2067") +
  theme_pub()
ggsave(file.path(OUT_DIR, "fig_hcv_trajectories.png"), p1, width = 11, height = 7, dpi = 300)

long_dc <- summary %>% select(scenario, year, median = dc_median, lo = dc_lo, hi = dc_hi) %>%
  mutate(measure = "DC")
long_hcc <- summary %>% select(scenario, year, median = hcc_median, lo = hcc_lo, hi = hcc_hi) %>%
  mutate(measure = "HCC")
long <- bind_rows(long_dc, long_hcc)
p2 <- ggplot(long, aes(x = year, y = median, colour = scenario)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = scenario), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~measure, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = okabe) + scale_fill_manual(values = okabe) +
  labs(x = "Year", y = "Count",
       title = "Decompensated cirrhosis and HCC projections, 2017-2067") +
  theme_pub()
ggsave(file.path(OUT_DIR, "fig_dc_hcc_trajectories.png"), p2, width = 11, height = 9, dpi = 300)

cat("Analysis outputs written to", OUT_DIR, "\n")
