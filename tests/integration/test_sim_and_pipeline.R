# Integration tests for the 480-compartment simulator and calibration wiring.
args0 <- commandArgs(FALSE)
file_arg <- sub("^--file=", "", args0[grepl("^--file=", args0)])
root <- normalizePath(file.path(dirname(normalizePath(file_arg)), "..", ".."))
source(file.path(root, "tests", "helper.R"))
suppressMessages({ library(Rcpp); library(RcppArmadillo) })

old <- setwd(file.path(root, "scripts", "model"))
e <- new.env(parent = globalenv()); sys.source("setup.R", e); setwd(old)
if (!exists("run_sim", envir = globalenv(), inherits = FALSE))
  Rcpp::sourceCpp(file.path(root, "scripts", "model", "sim.cpp"))
run_sim <- get("run_sim", envir = globalenv(), inherits = FALSE)
for (f in c("targets.R", "model_metrics.R", "equilibrium.R", "likelihood.R"))
  sys.source(file.path(root, "scripts", "calibration", f), envir = environment())
targets <- read_targets(file.path(root, "data", "calibration_targets.json"))

short_data <- e$data; short_data$t_end <- 0.1; short_data$dt <- 1 / 365
out <- run_sim(e$params, short_data)
check_equal(dim(out)[2], 481, "simulator returns time plus 480 states")
check_true(all(is.finite(out)) && all(out[, -1] >= -1e-8),
           "short simulation has finite non-negative compartments")
check_true(all(J_summary_at(out, nrow(out))$p_sero >= 0 &
               J_summary_at(out, nrow(out))$p_sero <= 1),
           "J seroprevalence remains in [0,1]")

simulate_seed <- function(state, scale = 1) {
  d <- short_data; d$y0 <- rep(0, 480)
  d$y0[idx(0, 1, state, 2) - 1L] <- 100 * scale
  run_sim(e$params, d)
}
a_out <- simulate_seed(1)
t_out <- simulate_seed(3)
check_true(a_out[nrow(a_out), idx(0, 1, 4, 2)] > 0,
           "acute infection flows to recovered")
check_true(t_out[nrow(t_out), idx(0, 1, 4, 2)] > 0,
           "treatment flows to recovered through SVR12")

# Age-specific treatment eligibility uses both bounds: with ages 20--39
# enabled, the 20--29 compartment recovers while the 40--49 compartment does
# not receive the treatment flow.
age_data <- short_data; age_data$y0 <- rep(0, 480); age_data$y0[idx(0, 1, 2, 1) - 1L] <- 100
age_data$y0[idx(0, 1, 2, 3) - 1L] <- 100
age_params <- e$params
age_params$q <- 0; age_params$beta <- rep(0, 6); age_params$beta_base <- e$params$beta
age_params$p_NC_CC <- 0; age_params$tau <- rep(1, 4); age_params$tau_stratum <- rep(1, 4)
age_params$tau_min_age <- 2; age_params$tau_max_age <- 3
age_params$tau_by_stratum <- as.vector(t(outer(age_params$tau_stratum, age_params$tau)))
age_out <- run_sim(age_params, age_data)
check_true(age_out[nrow(age_out), idx(0, 1, 4, 1)] > 0,
           "lower eligible age receives treatment")
check_true(age_out[nrow(age_out), idx(0, 1, 4, 3)] < 1e-3,
           "upper age bound excludes older compartment")
recovered_stage <- short_data; recovered_stage$y0 <- rep(0, 480)
recovered_stage$y0[idx(0, 2, 4, 2) - 1L] <- 100
stage_out <- run_sim(e$params, recovered_stage)
check_true(stage_out[nrow(stage_out), idx(0, 3, 4, 2)] > 0 &&
             stage_out[nrow(stage_out), idx(0, 4, 4, 2)] > 0,
           "recovered post-SVR progression enters later stages")

# Add an infectious active PWID and verify recovered people can be reinfected.
reinfection <- short_data; reinfection$y0 <- rep(0, 480)
reinfection$y0[idx(0, 1, 4, 2) - 1L] <- 100
reinfection$y0[idx(0, 1, 1, 2) - 1L] <- 100
reinfection_params <- e$params; reinfection_params$q <- 1
reinfection_params$C_contact <- matrix(100, 6, 6)
reinf_out <- run_sim(reinfection_params, reinfection)
recovered_only <- reinfection; recovered_only$y0[idx(0, 1, 1, 2) - 1L] <- 0
reinf_base <- run_sim(reinfection_params, recovered_only)
check_true(reinf_out[nrow(reinf_out), idx(0, 1, 1, 2)] >
             reinf_base[nrow(reinf_base), idx(0, 1, 1, 2)],
           "recovered people are susceptible to reinfection")

# The post-2015 schedule is active and uses the calibrated beta scale.
scheduled <- e$data; scheduled$t_start <- 150; scheduled$t_end <- 151
scheduled$y0 <- rep(0, 480); scheduled$y0[idx(0, 1, 0, 2) - 1L] <- 1
scheduled_params <- e$params; scheduled_params$beta <- rep(2, 6)
scheduled_params$beta_base <- e$params$beta
scheduled_out <- run_sim(scheduled_params, scheduled)
check_true(sum(scheduled_out[nrow(scheduled_out), -1]) > 1,
           "post-2015 scheduled inflow is active")

# Calendar alignment: model time 150 is the 2015 snapshot, so the interval
# beginning at 150 must use the 2016 row rather than the 2015 row.
calendar_data <- scheduled; calendar_data$inflow_schedule <- rbind(
  c(100, rep(0, 5)), c(10, rep(0, 5)))
calendar_params <- e$params; calendar_params$beta <- e$params$beta
calendar_2016 <- run_sim(calendar_params, calendar_data)
calendar_data$inflow_schedule <- rbind(c(10, rep(0, 5)), c(100, rep(0, 5)))
calendar_2015 <- run_sim(calendar_params, calendar_data)
check_true(sum(calendar_2016[nrow(calendar_2016), -1]) + 5 <
             sum(calendar_2015[nrow(calendar_2015), -1]),
           "post-2015 schedule starts with the 2016 row")

cal_data <- e$data; cal_data$t_end <- 10
obj <- make_objective(e$params, cal_data)
check_true(is.finite(obj(rep(0, N_THETA))), "likelihood is finite at theta zero")

n_fail <- test_report("tests/integration/test_sim_and_pipeline.R")
if (n_fail > 0) quit(status = 1)
