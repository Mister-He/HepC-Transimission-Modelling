# =============================================================================
# tests/integration/test_sim_and_pipeline.R — integration tests for the 480-compartment simulator
# (compilation, flows, reinfection, treatment, inflow) and calibration wiring.
# =============================================================================

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

# First arrest: run_sim moves D to J at c_first * lambda3.
fa_data <- short_data; fa_data$t_end <- 1; fa_data$y0 <- rep(0, 480); fa_data$y0[idx(0, 1, 0, 2) - 1L] <- 100
fa <- modifyList(e$params, list(q = 0, beta = rep(0, 6)))
fa1 <- run_sim(fa, fa_data); fa2 <- run_sim(modifyList(fa, list(c_first = 2 * fa$c_first)), fa_data)
check_true(fa2[nrow(fa2), idx(0, 1, 0, 2)] < fa1[nrow(fa1), idx(0, 1, 0, 2)] &&
             fa2[nrow(fa2), idx(1, 1, 0, 2)] > fa1[nrow(fa1), idx(1, 1, 0, 2)],
           "doubling c_first moves more of D into J")
check_equal(fa1[nrow(fa1), idx(0, 1, 0, 2)], 100 * exp(-(fa$c_first * fa$lambda3[3] + fa$mu[3] * fa$smr_nonviraemic + 0.1)),
            tol = 1e-4, "D leaves at c * lambda3 + mu * SMR + ageing")
check_error(run_sim(modifyList(fa, list(c_first = -1)), fa_data), "a negative c_first is rejected")
# PWID mortality: in J (no injecting, no release, no progression) an uninfected person leaves
# at mu * smr_nonviraemic + ageing and a chronic one at mu * smr_viraemic + ageing.
md <- fa_data; md$y0 <- rep(0, 480); md$y0[idx(1, 1, 0, 2) - 1L] <- 100; md$y0[idx(1, 1, 2, 3) - 1L] <- 100
mp <- modifyList(fa, list(lambda2 = rep(0, 6), p_NC_CC = 0, mu = rep(0.01, 6)))
mo <- run_sim(mp, md); last <- nrow(mo)
check_equal(c(mo[last, idx(1, 1, 0, 2)], mo[last, idx(1, 1, 2, 3)]),
            100 * exp(-0.01 * c(mp$smr_nonviraemic, mp$smr_viraemic) - 0.1),
            tol = 1e-4, "J deaths follow mu * SMR by viraemic status")
check_true(mp$smr_viraemic > mp$smr_nonviraemic, "viraemic people die faster than non-viraemic")

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
age_params$q <- 0; age_params$beta <- rep(0, 6)
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

# The post-2015 schedule is active (CNB counts, unscaled).
scheduled <- e$data; scheduled$t_start <- 150; scheduled$t_end <- 151
scheduled$y0 <- rep(0, 480); scheduled$y0[idx(0, 1, 0, 2) - 1L] <- 1
scheduled_params <- e$params; scheduled_params$beta <- rep(0, 6)
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

# Event-located annual quotas: F and X stop independently; existing courses
# finish after the cap, and restarting the solver resets the annual counter.
cap_data <- e$data; cap_data$t_start <- 160; cap_data$t_end <- 161
cap_data$y0 <- rep(0,480)
for (s in 1:3) cap_data$y0[idx(s,1,2,2)-1L] <- 100
cap_data$treatment_capacity <- c(Inf,Inf,5,10)
cap_pm <- e$params; cap_pm$q <- 0; cap_pm$beta <- rep(0,6)
cap_pm$lambda2 <- cap_pm$lambda3 <- rep(0,6)
cap_pm$tau_by_stratum <- c(rep(0,4),rep(10,12))
cap_pm$tau_min_age <- 1; cap_pm$tau_max_age <- 6
cap_out <- run_sim(cap_pm,cap_data)
cap_events <- cap_out[nrow(cap_out),482:493]
check_true(abs(cap_events[3]-5)<1e-6 && abs(cap_events[4]-10)<1e-6 &&
             cap_events[2]>10,"F/X independent quotas stop exactly; J continues")
check_true(cap_events[7]>0 && cap_events[8]>0 &&
             cap_events[11]>0 && cap_events[12]>0,
           "completion and cure continue after quota exhaustion")
cap_data$y0 <- cap_out[nrow(cap_out),2:481]
cap_data$t_start <- 161; cap_data$t_end <- 162
reset_out <- run_sim(cap_pm,cap_data)
check_true(abs(reset_out[nrow(reset_out),484]-5)<1e-6,
           "annual F quota resets with retained disease state")
cap_data$treatment_capacity <- c(Inf,Inf,0,0)
zero_out <- run_sim(cap_pm,cap_data)
check_true(all(abs(zero_out[nrow(zero_out),484:485])<1e-8) &&
             zero_out[nrow(zero_out),483]>0,
           "zero F/X quota leaves J treatment active")
balance_pm <- cap_pm; balance_pm$mu <- rep(0,6)
balance_pm$mu_DC <- 0; balance_pm$mu_HCC <- 0
balance_pm$p_NC_CC <- 0; balance_pm$p_CC_DC <- 0
balance_pm$p_CC_HCC <- 0; balance_pm$p_DC_HCC <- 0
balance_data <- cap_data; balance_data$t_start <- 160; balance_data$t_end <- 160.2
balance_data$y0 <- rep(0,480); balance_data$y0[idx(2,1,2,2)-1L] <- 100
balance_data$treatment_capacity <- c(Inf,Inf,5,Inf)
balance_out <- run_sim(balance_pm,balance_data)
balance_events <- balance_out[nrow(balance_out),482:493]
f_t <- sum(balance_out[nrow(balance_out),sapply(0:5,function(a) idx(2,1,3,a))])
check_true(abs(f_t-(balance_events[3]-balance_events[7]))<1e-5,
           "c-to-t start flow matches t inflow before treatment completions")

# The equilibrium gate compares the final state with T-10
# (scripts/calibration/equilibrium.R), so a horizon shorter than 10 years
# compares against the t = 0 seed, where the DC/HCC compartments are still zero
# and the penalty becomes non-finite. Calibration always runs model years 0--150
# (year 150 is the 2015 snapshot), so evaluate the objective on that horizon.
cal_data <- e$data
cal_data$t_end <- targets$equilibrium_horizon_years
obj <- make_objective(e$params, cal_data)
check_true(is.finite(obj(rep(0, N_THETA))), "likelihood is finite at theta zero")

cal_out <- run_sim(build_params(rep(0, N_THETA), e$params), cal_data)
check_true(check_equilibrium(cal_out)$pass,
           "base parameters satisfy the 10-year equilibrium gate")

# The MCMC target is the likelihood.R prior + Beta seroprevalence (no other terms).
sys.source(file.path(root, "scripts", "calibration", "mcmc.R"), envir = environment())
lp0 <- log_posterior_exact(rep(0, N_THETA), e$params, cal_data)$lp
s0 <- J_summary_final(cal_out)
check_equal(lp0, -nlp(rep(0, N_THETA)) - nll_prev(s0$p_sero),
            tol = 1e-8, "the MCMC log posterior is prior + prevalence")

n_fail <- test_report("tests/integration/test_sim_and_pipeline.R")
if (n_fail > 0) quit(status = 1)
