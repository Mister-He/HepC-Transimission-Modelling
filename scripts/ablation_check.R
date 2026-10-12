#!/usr/bin/env Rscript
# =============================================================================
# ablation_check.R — reference-MCMC ablation (docs/ablation_check.md): compares proposal
# covariance and initialisation strategies on the shared posterior and checks the
# integration step, thinning, adaptation drift and proposal conditioning.
#
# Phases (--phase, default all): prepare (inputs and predeclared config), hessian,
# wave1 (A_mean, A_disp, A_tail, D_prior), pilot, wave2 (B_pilot, C_map, dt_daily,
# dt_monthly), analyze. Completed arms are skipped. Outputs: output/ablation_check/<run-id>/.
# Usage: nohup Rscript scripts/ablation_check.R --run-id ab_20261012 > log 2>&1 &
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) { i <- match(name, args); if (is.na(i)) default else args[i + 1] }
ROOT <- normalizePath(arg("--root", getwd()), mustWork = TRUE); setwd(ROOT)
RUN_ID <- arg("--run-id", "ab_20261012")
PHASE <- arg("--phase", "all")
SMOKE <- "--smoke" %in% args
OUT <- file.path("output/ablation_check", RUN_ID)
FIT <- "output/calibration/r_20261010_agepop/fit.rds"
NPE <- "output/calibration/npe_lhs/posterior_samples_npe.csv"
REF <- "output/calibration/mcmc_reference"
TH <- paste0("theta", 1:12); D <- 12; C_D <- 2.38^2 / D
B <- if (SMOKE) list(warm = 120, adapt = 60, every = 20, evals = 300, pilot = 150, pilot_warm = 80, pilot_adapt = 30,
                     dt_warm = 60, dt_adapt = 30, dt_evals = 150, n_det = 20) else
  list(warm = 10000, adapt = 5000, every = 50, evals = 30000, pilot = 3000, pilot_warm = 1500, pilot_adapt = 500,
       dt_warm = 2000, dt_adapt = 500, dt_evals = 6000, n_det = 500)
dir.create(file.path(OUT, "inputs"), recursive = TRUE, showWarnings = FALSE)
inp <- function(f) file.path(OUT, "inputs", f)
logmsg <- function(...) cat(format(Sys.time()), ..., "\n")

# ---------------- prepare: inputs and predeclared configuration ----------------
if (PHASE %in% c("all", "prepare") && !file.exists(file.path(OUT, "ablation_config.json"))) {
  npe <- as.matrix(read.csv(NPE)[, TH]); mu <- colMeans(npe); S <- cov(npe)
  write.csv(S, inp("cov_npe.csv"), row.names = FALSE)
  write.csv(as.data.frame(matrix(mu, 4, D, byrow = TRUE, dimnames = list(NULL, TH))), inp("starts_npe_mean.csv"), row.names = FALSE)
  set.seed(20261021); disp <- npe[sample.int(nrow(npe), 3), ]
  over <- mu + 3 * sqrt(diag(S)) * sign(rnorm(D))          # one overdispersed start, 3 NPE SD from the mean
  write.csv(as.data.frame(rbind(disp, over)), inp("starts_npe_dispersed.csv"), row.names = FALSE)
  md <- mahalanobis(npe, mu, S); tail_ix <- which(md >= quantile(md, 0.995))
  dirs <- sign(sweep(npe[tail_ix, ], 2, mu)); pick <- tail_ix[1]
  for (k in 2:4) {                                          # tail draws in different directions
    cand <- setdiff(tail_ix, pick); sim <- sapply(cand, function(i) max(abs(dirs[match(i, tail_ix), ] %*% t(dirs[match(pick, tail_ix), , drop = FALSE]))))
    pick <- c(pick, cand[which.min(sim)])
  }
  write.csv(as.data.frame(npe[pick, ]), inp("starts_npe_tail.csv"), row.names = FALSE)
  set.seed(20261022); write.csv(as.data.frame(npe[sample.int(nrow(npe), 4), ]), inp("starts_npe_random.csv"), row.names = FALSE)
  set.seed(20261023); write.csv(as.data.frame(npe[sample.int(nrow(npe), 4), ]), inp("starts_pilot.csv"), row.names = FALSE)
  map <- as.numeric(readRDS(FIT)$best$theta)
  write.csv(as.data.frame(matrix(map, 4, D, byrow = TRUE, dimnames = list(NULL, TH))), inp("starts_map.csv"), row.names = FALSE)
  md5 <- function(f) unname(tools::md5sum(f))
  jsonlite::write_json(list(
    run_id = RUN_ID, created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), smoke = SMOKE,
    git_commit = system("git rev-parse HEAD", intern = TRUE),
    git_dirty_scripts = length(system("git status --porcelain scripts", intern = TRUE)) > 0,
    md5 = list(setup.R = md5("scripts/model/setup.R"), sim.cpp = md5("scripts/model/sim.cpp"),
               mcmc.R = md5("scripts/calibration/mcmc.R"), run_mcmc_reference.R = md5("scripts/run_mcmc_reference.R"),
               ablation_check.R = md5("scripts/ablation_check.R")),
    target = "Normal(fit$best$theta, 0.20) prior; Binomial + log-normal likelihood; equilibrium as hard restriction (same as NPE)",
    fit = FIT, npe_posterior = NPE, reference_posterior = REF,
    budget_per_chain = list(total_evaluations = B$evals, warmup = B$warm, adapt_start = B$adapt, adapt_every = B$every,
                            pilot_evaluations = B$pilot, pilot_charged_to_arm_B = TRUE,
                            hessian = "evaluations counted and reported separately for arm C"),
    arms = list(
      A_mean = "cov = (2.38^2/d) cov(NPE); starts = NPE mean (all chains)  [control]",
      A_disp = "cov as A; starts = 3 random NPE draws + 1 overdispersed (mean +/- 3 NPE SD)",
      A_tail = "cov as A; starts = 4 NPE draws from the top 0.5% Mahalanobis distance, different directions",
      B_pilot = "pilot: 4 chains x 3,000 evals (prior-scaled cov, adapt from 500); main cov = (2.38^2/d) pilot cov; starts = pilot end states",
      C_map = "cov = inverse Hessian of -log posterior at fit$best$theta (numDeriv); starts = MAP",
      D_prior = "cov = 0.1^2 diag(prior variance) (production setting); starts = 4 random NPE draws",
      dt_daily = "short chains, dt = 1/365, cov as A, NPE mean starts", dt_monthly = "as dt_daily with dt = 1/12"),
    thinning_rule = "largest thin keeping >= 90% of unthinned bulk and tail ESS on the reference chains; chosen thin = 10 (storage only; diagnostics use all draws)",
    thresholds = list(rhat_max = 1.01, ess_bulk_min = 400, ess_tail_min = 400, acceptance = c(0.15, 0.50),
                      condition_number_flag = 1e6, adaptation_drift_tolerance = 0.01,
                      drift_note = "each update adds 50 draws to a window of >= 5,000, so ~1% change is expected noise; 1e-3 also reported",
                      chain_outlier_sd = 0.5, start_sensitivity = "|mean_arm - mean_ref| <= 3 * sqrt(mcse_arm^2 + mcse_ref^2) for every parameter",
                      decision = "best min(bulk ESS) per CPU-hour among arms passing all checks")
  ), file.path(OUT, "ablation_config.json"), pretty = TRUE, auto_unbox = TRUE, digits = NA)
  logmsg("prepared inputs and predeclared config")
}

# ---------------- helpers ----------------
run_arm <- function(name, extra, wait = TRUE) {
  dir <- file.path(OUT, name)
  if (file.exists(file.path(dir, "MCMC_STATUS.txt"))) { logmsg("skip", name); return(invisible(NULL)) }
  dir.create(dir, showWarnings = FALSE)
  a <- c("scripts/run_mcmc_reference.R", "--fit", FIT, "--out-dir", dir, "--cores", 4, "--n-chains", 4, extra)
  writeLines(paste("Rscript", paste(a, collapse = " ")), file.path(dir, "command.txt"))
  logmsg("start", name)
  system2("Rscript", a, stdout = file.path(dir, "run.log"), stderr = file.path(dir, "run.log"), wait = wait)
}
await <- function(names) {
  repeat {
    done <- file.exists(file.path(OUT, names, "MCMC_STATUS.txt"))
    if (all(done)) break
    Sys.sleep(60)
  }
  logmsg("finished", paste(names, collapse = ", "))
}
main_args <- function(cov = NULL, scale = 1, starts, seed, evals = B$evals) c(
  "--n-warmup", B$warm, "--adapt-start", B$adapt, "--adapt-every", B$every, "--max-evals", evals,
  "--thin", 10, "--seed", seed, if (!is.null(cov)) c("--cov-file", cov, "--cov-scale", scale), "--start-file", starts)

# ---------------- hessian (arm C input) ----------------
if (PHASE %in% c("all", "hessian") && !file.exists(inp("cov_hessian.csv"))) {
  setup_env <- new.env(); old <- setwd("scripts/model"); sys.source("setup.R", setup_env); setwd(old)
  Rcpp::sourceCpp("scripts/model/sim.cpp")
  source("scripts/calibration/targets.R"); targets <- read_targets("data/calibration_targets.json")
  for (f in c("model_metrics.R", "equilibrium.R", "likelihood.R", "mcmc.R")) source(file.path("scripts/calibration", f))
  fit <- readRDS(FIT); map <- as.numeric(fit$best$theta); dat <- setup_env$data; dat$t_start <- 0; dat$t_end <- 150
  n_eval <- 0L; t0 <- proc.time()
  nlp <- function(th) { n_eval <<- n_eval + 1L; -log_posterior_exact(th, setup_env$params, dat, map, rep(0.2, D))$lp }
  H <- numDeriv::hessian(nlp, map, method.args = list(eps = 1e-3, d = 1e-3, r = 2))
  H <- (H + t(H)) / 2; ev <- eigen(H, symmetric = TRUE)$values
  Sig <- if (all(ev > 0)) solve(H) else as.matrix(Matrix::nearPD(solve(H))$mat)
  write.csv(Sig, inp("cov_hessian.csv"), row.names = FALSE)
  jsonlite::write_json(list(map = map, map_source = "fit$best$theta (Nelder-Mead best fit; prior centre)",
                            hessian_method = "numDeriv::hessian (Richardson, r = 2, d = 1e-3) of -log posterior",
                            hessian_eigen = ev, positive_definite = all(ev > 0), evaluations = n_eval,
                            seconds = (proc.time() - t0)[["elapsed"]],
                            note = "local Gaussian approximation at the deterministic best fit, which is not exactly the posterior mode"),
                       inp("hessian_info.json"), pretty = TRUE, auto_unbox = TRUE, digits = NA)
  logmsg("hessian done:", n_eval, "evaluations")
}

# ---------------- wave 1 ----------------
if (PHASE %in% c("all", "wave1")) {
  run_arm("A_mean", main_args(inp("cov_npe.csv"), C_D, inp("starts_npe_mean.csv"), 20261111), wait = FALSE)
  run_arm("A_disp", main_args(inp("cov_npe.csv"), C_D, inp("starts_npe_dispersed.csv"), 20261112), wait = FALSE)
  run_arm("A_tail", main_args(inp("cov_npe.csv"), C_D, inp("starts_npe_tail.csv"), 20261113), wait = FALSE)
  run_arm("D_prior", main_args(NULL, 1, inp("starts_npe_random.csv"), 20261114), wait = FALSE)
  await(c("A_mean", "A_disp", "A_tail", "D_prior"))
}

# ---------------- pilot (arm B input) ----------------
if (PHASE %in% c("all", "pilot")) {
  run_arm("B_pilot_stage1", c("--n-warmup", B$pilot_warm, "--adapt-start", B$pilot_adapt, "--adapt-every", B$every,
                              "--max-evals", B$pilot, "--thin", 1, "--seed", 20261115, "--start-file", inp("starts_pilot.csv")))
  p <- read.csv(file.path(OUT, "B_pilot_stage1", "mcmc_draws.csv"))
  half <- p[p$iteration > max(p$iteration) / 2, TH]
  write.csv(cov(half), inp("cov_pilot.csv"), row.names = FALSE)
  ends <- do.call(rbind, lapply(split(p, p$chain), function(x) x[which.max(x$iteration), TH]))
  write.csv(ends, inp("starts_pilot_end.csv"), row.names = FALSE)
}

# ---------------- wave 2 ----------------
if (PHASE %in% c("all", "wave2")) {
  run_arm("B_pilot", main_args(inp("cov_pilot.csv"), C_D, inp("starts_pilot_end.csv"), 20261116, evals = B$evals - B$pilot), wait = FALSE)
  run_arm("C_map", main_args(inp("cov_hessian.csv"), 1, inp("starts_map.csv"), 20261117), wait = FALSE)
  dt_args <- function(dt, seed) c("--n-warmup", B$dt_warm, "--adapt-start", B$dt_adapt, "--adapt-every", B$every,
    "--max-evals", B$dt_evals, "--thin", 1, "--seed", seed, "--cov-file", inp("cov_npe.csv"), "--cov-scale", C_D,
    "--start-file", inp("starts_npe_mean.csv"), "--dt", dt)
  run_arm("dt_daily", dt_args(1 / 365, 20261118), wait = FALSE)
  run_arm("dt_monthly", dt_args(1 / 12, 20261118), wait = FALSE)
  await(c("B_pilot", "C_map", "dt_daily", "dt_monthly"))
}

# ---------------- analyze ----------------
if (PHASE %in% c("all", "analyze")) {
  cfg <- jsonlite::fromJSON(file.path(OUT, "ablation_config.json")); thr <- cfg$thresholds
  arms <- c("A_mean", "A_disp", "A_tail", "B_pilot", "C_map", "D_prior")
  ref_d <- read.csv(file.path(REF, "mcmc_diagnostics.csv")); ref_d <- ref_d[ref_d$subset == "all_retained", ]
  npe <- as.matrix(read.csv(NPE)[, TH])
  rd <- function(a, f) read.csv(file.path(OUT, a, f))
  # Diagnostics per parameter per arm
  diag <- do.call(rbind, lapply(c(arms, "dt_daily", "dt_monthly"), function(a) {
    d <- rd(a, "mcmc_diagnostics.csv"); cs <- rd(a, "mcmc_chain_summary.csv")
    cbind(arm = a, d[, c("parameter", "mean", "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean", "mcse_q025", "mcse_q975", "iat")],
          acceptance = mean(cs$acceptance))
  }))
  write.csv(diag, file.path(OUT, "diagnostics.csv"), row.names = FALSE)
  # Arm summary: budgets, CPU, ESS per CPU-hour, checks
  hess <- jsonlite::fromJSON(inp("hessian_info.json"))
  arm_sum <- do.call(rbind, lapply(arms, function(a) {
    d <- diag[diag$arm == a, ]; cs <- rd(a, "mcmc_chain_summary.csv"); dr <- rd(a, "covariance_drift.csv")
    extra_cpu <- if (a == "B_pilot") sum(rd("B_pilot_stage1", "mcmc_chain_summary.csv")$cpu_sec) else if (a == "C_map") hess$seconds else 0
    extra_evals <- if (a == "B_pilot") sum(rd("B_pilot_stage1", "mcmc_chain_summary.csv")$total_evals) else if (a == "C_map") hess$evaluations else 0
    last_drift <- do.call(rbind, lapply(split(dr, dr$chain), function(x) x[which.max(x$iteration), ]))
    first_below <- function(tol) max(vapply(split(dr, dr$chain), function(x) {
      ok <- which(!is.na(x$drift) & x$drift < tol); if (length(ok)) x$iteration[min(ok)] else NA_real_ }, numeric(1)))
    pc <- rd(a, "mcmc_per_chain.csv"); eig <- rd(a, "covariance_eigen.csv")
    cpu_h <- (sum(cs$cpu_sec) + extra_cpu) / 3600
    data.frame(arm = a, evaluations_main = sum(cs$total_evals), evaluations_extra = extra_evals,
               cpu_hours = cpu_h, wall_hours = max(cs$total_sec) / 3600,
               max_rhat = max(d$rhat), min_ess_bulk = min(d$ess_bulk), min_ess_tail = min(d$ess_tail),
               ess_bulk_per_cpu_hour = min(d$ess_bulk) / cpu_h,
               acceptance_min = min(cs$acceptance), acceptance_max = max(cs$acceptance),
               final_warmup_drift_max = max(last_drift$drift, na.rm = TRUE),
               iteration_drift_below_0.01 = first_below(0.01), iteration_drift_below_0.001 = first_below(0.001),
               max_chain_median_shift_sd = max(abs(pc$median_shift_sd)),
               cond_initial = eig$condition_number[eig$matrix_label == "initial_proposal"],
               cond_adapted_max = max(eig$condition_number[eig$matrix_label == "adapted_proposal"]),
               max_abs_smd_vs_reference = max(abs((d$mean - ref_d$mean) / ref_d$sd)),
               start_sensitivity_ok = all(abs(d$mean - ref_d$mean) <= 3 * sqrt(d$mcse_mean^2 + ref_d$mcse_mean^2)))
  }))
  arm_sum$pass_rhat <- arm_sum$max_rhat < thr$rhat_max
  arm_sum$pass_ess <- arm_sum$min_ess_bulk > thr$ess_bulk_min & arm_sum$min_ess_tail > thr$ess_tail_min
  arm_sum$pass_acceptance <- arm_sum$acceptance_min >= thr$acceptance[1] & arm_sum$acceptance_max <= thr$acceptance[2]
  arm_sum$pass_drift <- arm_sum$final_warmup_drift_max < thr$adaptation_drift_tolerance
  arm_sum$pass_chains <- arm_sum$max_chain_median_shift_sd <= thr$chain_outlier_sd
  arm_sum$pass_all <- with(arm_sum, pass_rhat & pass_ess & pass_acceptance & pass_drift & pass_chains & start_sensitivity_ok)
  write.csv(arm_sum, file.path(OUT, "arm_summary.csv"), row.names = FALSE)
  # Start sensitivity (Part C)
  ss <- do.call(rbind, lapply(c("A_mean", "A_disp", "A_tail", "D_prior"), function(a) {
    d <- diag[diag$arm == a, ]
    data.frame(arm = a, parameter = d$parameter, mean = d$mean, ref_mean = ref_d$mean,
               diff_sd = (d$mean - ref_d$mean) / ref_d$sd,
               diff_mcse = (d$mean - ref_d$mean) / sqrt(d$mcse_mean^2 + ref_d$mcse_mean^2),
               sd_ratio = d$sd / ref_d$sd, rhat = d$rhat, ess_bulk = d$ess_bulk)
  }))
  write.csv(ss, file.path(OUT, "start_sensitivity.csv"), row.names = FALSE)
  # Covariance eigen-spectra (B4) and drift (B3)
  eigrow <- function(label, S) { e <- eigen(as.matrix(S), symmetric = TRUE, only.values = TRUE)$values
    data.frame(arm = "shared", matrix_label = label, chain = NA, t(setNames(e, paste0("eigen", 1:12))), condition_number = max(e) / min(e), min_eigen = min(e)) }
  eig_all <- rbind(eigrow("cov_npe", read.csv(inp("cov_npe.csv"))), eigrow("cov_hessian_inverse", read.csv(inp("cov_hessian.csv"))),
                   eigrow("cov_pilot", read.csv(inp("cov_pilot.csv"))),
                   eigrow("cov_reference_posterior", cov(as.matrix(read.csv(file.path(REF, "mcmc_draws.csv"))[, TH]))),
                   do.call(rbind, lapply(arms, function(a) cbind(arm = a, rd(a, "covariance_eigen.csv")))))
  write.csv(eig_all, file.path(OUT, "covariance_eigen.csv"), row.names = FALSE)
  write.csv(do.call(rbind, lapply(arms, function(a) cbind(arm = a, rd(a, "covariance_drift.csv")))),
            file.path(OUT, "covariance_drift.csv"), row.names = FALSE)
  # B1: NPE vs MCMC scale
  ref_draws <- as.matrix(read.csv(file.path(REF, "mcmc_draws.csv"))[, TH])
  q <- function(m, p) apply(m, 2, quantile, p)
  b1 <- data.frame(parameter = TH, npe_sd = apply(npe, 2, sd), mcmc_sd = ref_d$sd, ratio_npe_over_mcmc = apply(npe, 2, sd) / ref_d$sd,
                   npe_width95 = q(npe, .975) - q(npe, .025), mcmc_width95 = q(ref_draws, .975) - q(ref_draws, .025),
                   hessian_sd = sqrt(diag(as.matrix(read.csv(inp("cov_hessian.csv"))))),
                   arm_A_mean_ess_bulk = diag$ess_bulk[diag$arm == "A_mean"], arm_A_mean_acceptance = diag$acceptance[diag$arm == "A_mean"])
  b1$width_ratio <- b1$npe_width95 / b1$mcmc_width95
  write.csv(b1, file.path(OUT, "npe_vs_mcmc_scale.csv"), row.names = FALSE)
  # A1: daily vs monthly step
  dd <- diag[diag$arm == "dt_daily", ]; dm <- diag[diag$arm == "dt_monthly", ]
  dtc <- data.frame(parameter = TH, mean_daily = dd$mean, mean_monthly = dm$mean, sd_daily = dd$sd, sd_monthly = dm$sd,
                    diff_sd = (dm$mean - dd$mean) / dd$sd, diff_mcse = (dm$mean - dd$mean) / sqrt(dd$mcse_mean^2 + dm$mcse_mean^2),
                    cpu_hours_daily = sum(rd("dt_daily", "mcmc_chain_summary.csv")$cpu_sec) / 3600,
                    cpu_hours_monthly = sum(rd("dt_monthly", "mcmc_chain_summary.csv")$cpu_sec) / 3600)
  # Deterministic discretisation check at reference posterior draws
  setup_env <- new.env(); old <- setwd("scripts/model"); sys.source("setup.R", setup_env); setwd(old)
  Rcpp::sourceCpp("scripts/model/sim.cpp"); source("scripts/calibration/targets.R"); targets <- read_targets("data/calibration_targets.json")
  for (f in c("model_metrics.R", "equilibrium.R", "likelihood.R", "mcmc.R")) source(file.path("scripts/calibration", f))
  fit <- readRDS(FIT); map <- as.numeric(fit$best$theta)
  set.seed(20261024); det_idx <- sample.int(nrow(ref_draws), B$n_det)
  one <- function(i, dt) { dat <- setup_env$data; dat$t_start <- 0; dat$t_end <- 150; dat$dt <- dt
    out <- run_sim(build_params(ref_draws[i, ], setup_env$params), dat); s <- J_summary_final(out)
    c(s$p_sero, s$N_hat, lp = log_posterior_exact(ref_draws[i, ], setup_env$params, dat, map, rep(0.2, D))$lp) }
  det <- do.call(rbind, parallel::mclapply(det_idx, function(i) c(one(i, 1 / 365), one(i, 1 / 12)), mc.cores = 14))
  k <- 13
  det_sum <- data.frame(quantity = c(paste0("p", 1:6), paste0("N", 1:6), "log_posterior"),
    mean_abs_diff = colMeans(abs(det[, (k + 1):(2 * k)] - det[, 1:k])),
    max_abs_diff = apply(abs(det[, (k + 1):(2 * k)] - det[, 1:k]), 2, max),
    mean_daily = colMeans(det[, 1:k]))
  write.csv(dtc, file.path(OUT, "dt_chain_comparison.csv"), row.names = FALSE)
  write.csv(det_sum, file.path(OUT, "dt_deterministic_comparison.csv"), row.names = FALSE)
  # Thinning evidence on the reference chains (A2)
  rdraw <- read.csv(file.path(REF, "mcmc_draws.csv")); n <- min(table(rdraw$chain)); arr <- array(NA, c(n, 4, D))
  for (c4 in 1:4) arr[, c4, ] <- as.matrix(rdraw[rdraw$chain == c4, TH][1:n, ])
  thin_tab <- do.call(rbind, lapply(c(1, 5, 10, 25, 50, 100, 250, 1000), function(t) { ix <- seq(1, n, by = t)
    data.frame(thin = t, retained = 4 * length(ix), min_ess_bulk = min(sapply(1:D, function(j) posterior::ess_bulk(arr[ix, , j]))),
               min_ess_tail = min(sapply(1:D, function(j) posterior::ess_tail(arr[ix, , j]))),
               max_rhat = max(sapply(1:D, function(j) posterior::rhat(arr[ix, , j])))) }))
  write.csv(thin_tab, file.path(OUT, "thinning_evidence.csv"), row.names = FALSE)
  # Manifest
  jsonlite::write_json(list(config = cfg, arms = lapply(setNames(c(arms, "B_pilot_stage1", "dt_daily", "dt_monthly"), c(arms, "B_pilot_stage1", "dt_daily", "dt_monthly")),
    function(a) list(command = readLines(file.path(OUT, a, "command.txt")), run_config = jsonlite::fromJSON(file.path(OUT, a, "mcmc_run_config.json")),
                     chains = rd(a, "mcmc_chain_summary.csv"))), hessian = hess,
    hardware = list(cpu = system("sysctl -n machdep.cpu.brand_string", intern = TRUE), cores = parallel::detectCores(),
                    concurrency = "4 arms x 4 chains per wave on 16 cores; CPU time per chain is the comparable cost")),
    file.path(OUT, "run_manifest.json"), pretty = TRUE, auto_unbox = TRUE, digits = NA)
  print(arm_sum[, c("arm", "max_rhat", "min_ess_bulk", "min_ess_tail", "ess_bulk_per_cpu_hour", "acceptance_min", "acceptance_max",
                    "final_warmup_drift_max", "max_abs_smd_vs_reference", "start_sensitivity_ok", "pass_all")])
  logmsg("analysis written to", OUT)
}
