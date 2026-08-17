# =============================================================================
# run_npe.R — Neural Posterior Estimation (NPE) + strict MCMC validation
#
# Pipeline (per prompt_mcmc.md):
#   data    generate prior training data + SBC held-out set (R, parallel sims)
#   train   train NPE_C (Python sbi) + posterior sampling + SBC + 2 seeds
#   mcmc    validation adaptive Metropolis (strict: R-hat in [0.99,1.01],
#           ESS > 400 per parameter); chains extended until criteria pass
#   predict posterior predictive credible intervals (NPE + MCMC + Laplace +
#           observed) with equilibrium filter
#   plots   ggplot2 figures (trace, density, predictive, CI compare, SBC)
#   all     run every step
#
# Usage:
#   Rscript src/calibration/run_npe.R --step all \
#     --root . --fit output/calibration/run2_v1_12p_warm/fit.rds \
#     --out-dir output/calibration/npe_bayes \
#     --n-sims 60000 --n-cores 6 --seed 2026 \
#     --python /tmp/bayes-venv/bin/python
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(name, default = NULL) {
  pos <- match(name, args)
  if (is.na(pos)) return(default)
  args[pos + 1]
}

STEP     <- arg_val("--step", "all")
ROOT     <- normalizePath(arg_val("--root", getwd()), mustWork = TRUE)
FIT_PATH <- normalizePath(arg_val("--fit",
              file.path(ROOT, "output", "calibration",
                        "run2_v1_12p_warm", "fit.rds")), mustWork = TRUE)
OUT_DIR  <- normalizePath(arg_val("--out-dir",
              file.path(ROOT, "output", "calibration", "npe_bayes")),
            mustWork = FALSE)
N_SIMS   <- as.integer(arg_val("--n-sims", "60000"))
N_CORES  <- as.integer(arg_val("--n-cores", "6"))
SEED     <- as.integer(arg_val("--seed", "2026"))
N_DRAWS  <- as.integer(arg_val("--n-draws", "60000"))
N_PROPOSAL <- as.integer(arg_val("--n-proposal", "20000"))
N_ROUNDS <- as.integer(arg_val("--n-rounds", "3"))
N_ITER_MCMC <- as.integer(arg_val("--n-iter-mcmc", "30000"))
BURNIN_MCMC <- as.integer(arg_val("--burnin-mcmc", "5000"))
MCMC_MAX_ITER <- as.integer(arg_val("--mcmc-max-iter", "400000"))
THIN_MCMC <- as.integer(arg_val("--thin-mcmc", "20"))
PYTHON   <- arg_val("--python", "/tmp/bayes-venv/bin/python")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("== NPE Bayesian run ==\n")
cat("Step:", STEP, "| Root:", ROOT, "\nOut:", OUT_DIR, "\n")

# ---------------------------------------------------------------------------
# 1. Load model + calibration modules (top level, like the other runners)
# ---------------------------------------------------------------------------
old_wd <- setwd(file.path(ROOT, "src"))
setup_env <- new.env(parent = globalenv())
sys.source("setup.R", envir = setup_env)
setwd(old_wd)
if (!exists("run_sim", envir = globalenv(), inherits = FALSE)) {
  Rcpp::sourceCpp(file.path(ROOT, "src", "sim.cpp"))
}
run_sim <- get("run_sim", envir = globalenv(), inherits = FALSE)
for (f in c("targets.R", "model_metrics.R", "equilibrium.R",
            "likelihood.R", "calibrate_nm.R", "mcmc.R", "plot_mcmc.R")) {
  sys.source(file.path(ROOT, "src", "calibration", f), envir = environment())
}
fit <- readRDS(FIT_PATH)
TARGET_TIME <- 45
base_params <- setup_env$params
data_local  <- setup_env$data
data_local$t_start <- fit$t_start
data_local$t_end   <- fit$t_end
priors <- make_priors()

# ---------------------------------------------------------------------------
# 2. Helpers
# ---------------------------------------------------------------------------
sample_prior <- function(n, priors, seed) {
  set.seed(seed)
  th <- matrix(NA_real_, n, 12)
  for (i in 1:6) {
    th[, i] <- rnorm(n, log(priors$contact_anchor[i]), priors$contact_sd)
  }
  z <- matrix(rt(n * 6, df = priors$beta_df), n, 6)
  th[, 7:12] <- matrix(rep(log(priors$beta_anchor), each = n), nrow = n,
                       ncol = 6) + priors$beta_sd * z
  th
}

sim_summary <- function(theta, base_params, data_local) {
  pm <- tryCatch(build_params(theta, base_params), error = function(e) NULL)
  if (is.null(pm)) return(NULL)
  out <- tryCatch(run_sim(pm, data_local), error = function(e) NULL)
  if (is.null(out)) return(NULL)
  s <- tryCatch(J_summary_final(out), error = function(e) NULL)
  if (is.null(s)) return(NULL)
  if (!all(is.finite(s$p_sero)) || !all(is.finite(s$N_hat)) ||
      any(s$N_hat <= 0)) return(NULL)
  c(p = s$p_sero, N = s$N_hat)
}

sim_parallel <- function(theta_mat, base_params, data_local, n_cores) {
  idx <- seq_len(nrow(theta_mat))
  run_one <- function(j) sim_summary(theta_mat[j, ], base_params, data_local)
  res <- if (n_cores > 1 && requireNamespace("parallel", quietly = TRUE)) {
    parallel::mclapply(idx, run_one, mc.cores = n_cores)
  } else {
    lapply(idx, run_one)
  }
  ok <- !sapply(res, is.null)
  list(x = do.call(rbind, res[ok]), ok = ok)
}

# ---------------------------------------------------------------------------
# STEP data — training + SBC held-out sets
# ---------------------------------------------------------------------------
step_data <- function() {
  cat("Generating training data:", N_SIMS, "prior draws...\n")
  th_train <- sample_prior(N_SIMS, priors, SEED)
  sim <- sim_parallel(th_train, base_params, data_local, N_CORES)
  cat("  finite simulations:", sum(sim$ok), "/", N_SIMS, "\n")
  write.csv(th_train[sim$ok, , drop = FALSE],
            file.path(OUT_DIR, "theta_train.csv"), row.names = FALSE)
  write.csv(sim$x, file.path(OUT_DIR, "x_train.csv"), row.names = FALSE)

  cat("Generating SBC held-out set (300 draws)...\n")
  th_sbc <- sample_prior(300, priors, SEED + 100)
  sim_sbc <- sim_parallel(th_sbc, base_params, data_local, N_CORES)
  cat("  finite SBC simulations:", sum(sim_sbc$ok), "/", 300, "\n")
  write.csv(th_sbc[sim_sbc$ok, , drop = FALSE],
            file.path(OUT_DIR, "theta_sbc.csv"), row.names = FALSE)
  write.csv(sim_sbc$x, file.path(OUT_DIR, "x_sbc.csv"), row.names = FALSE)

  x_obs <- data.frame(
    p1 = cal_targets$prev_binom[1], p2 = cal_targets$prev_binom[2],
    p3 = cal_targets$prev_binom[3], p4 = cal_targets$prev_binom[4],
    p5 = cal_targets$prev_binom[5], p6 = cal_targets$prev_binom[6],
    N1 = cal_targets$prison_total[1], N2 = cal_targets$prison_total[2],
    N3 = cal_targets$prison_total[3], N4 = cal_targets$prison_total[4],
    N5 = cal_targets$prison_total[5], N6 = cal_targets$prison_total[6]
  )
  write.csv(x_obs, file.path(OUT_DIR, "x_obs.csv"), row.names = FALSE)
  cat("Saved theta_train/x_train/theta_sbc/x_sbc/x_obs to", OUT_DIR, "\n")
}

# ---------------------------------------------------------------------------
# STEP train — call Python NPE
# ---------------------------------------------------------------------------
step_train <- function() {
  if (!file.exists(PYTHON)) stop("python not found:", PYTHON)
  py <- file.path(ROOT, "src", "calibration", "npe_train.py")
  base_args <- c("--theta", file.path(OUT_DIR, "theta_train.csv"),
                 "--x", file.path(OUT_DIR, "x_train.csv"),
                 "--x-obs", file.path(OUT_DIR, "x_obs.csv"),
                 "--theta-sbc", file.path(OUT_DIR, "theta_sbc.csv"),
                 "--x-sbc", file.path(OUT_DIR, "x_sbc.csv"),
                 "--out", OUT_DIR)
  run_py <- function(seed, suffix = "", round = 1) {
    cat("NPE round", round, "(seed", seed, ")...\n")
    sname <- if (nzchar(suffix)) paste0("_", suffix) else ""
    args <- c(py, "--theta", file.path(OUT_DIR, "theta_train.csv"),
              "--x", file.path(OUT_DIR, "x_train.csv"),
              "--x-obs", file.path(OUT_DIR, "x_obs.csv"),
              "--theta-sbc", file.path(OUT_DIR, "theta_sbc.csv"),
              "--x-sbc", file.path(OUT_DIR, "x_sbc.csv"),
              "--out", OUT_DIR, "--seed", seed, "--n-draws", N_DRAWS,
              "--round", round, "--n-rounds", N_ROUNDS,
              "--n-proposal", N_PROPOSAL,
              if (nzchar(suffix)) c("--suffix", suffix) else character(0))
    if (round > 1) {
      args <- c(args, "--theta2", file.path(OUT_DIR,
                paste0("proposal_theta_round", round, sname, ".csv")),
                "--x2", file.path(OUT_DIR,
                paste0("x_round", round, sname, ".csv")))
    }
    out <- system2(PYTHON, args, stdout = TRUE, stderr = TRUE,
                   env = c("MPLCONFIGDIR=/tmp/mplcache"))
    status <- attr(out, "status")
    if (!is.null(status) && status != 0) {
      cat(tail(out, 20), sep = "\n")
      stop("NPE round ", round, " failed (seed ", seed, ")")
    }
    cat("  NPE round", round, "seed", seed, "done.\n")
  }
  simulate_round <- function(suffix = "", round = 2) {
    sname <- if (nzchar(suffix)) paste0("_", suffix) else ""
    prop <- read.csv(file.path(OUT_DIR,
                               paste0("proposal_theta_round", round,
                                      sname, ".csv")))
    cat("Simulating round-", round, " proposals (", nrow(prop), "draws)...\n",
        sep = "")
    sim <- sim_parallel(as.matrix(prop), base_params, data_local, N_CORES)
    n_invalid <- sum(!sim$ok)
    x_out <- matrix(NA_real_, nrow(prop), 12)
    x_out[sim$ok, ] <- sim$x
    if (n_invalid > 0) {
      # sbi multi-round does not allow NaN/Inf simulations; replace invalid
      # summaries with sentinels far outside the observed range (p -> 1e-6,
      # N -> 1e12; after the logit/log + standardisation + clip they map to
      # the transform extremes). Documented in bayes_methodology.md.
      x_out[!sim$ok, ] <- matrix(rep(c(rep(1e-6, 6), rep(1e12, 6)), n_invalid),
                                 nrow = n_invalid, byrow = TRUE)
    }
    write.csv(x_out, file.path(OUT_DIR,
                               paste0("x_round", round, sname, ".csv")),
              row.names = FALSE)
    cat("  round-", round, " x written (valid", sum(sim$ok), "/", nrow(prop),
        "; sentinel-filled", n_invalid, ").\n")
  }

  for (s in list(list(seed = SEED, suffix = ""),
                 list(seed = SEED + 1, suffix = "seed2"))) {
    for (r in seq_len(N_ROUNDS)) {
      run_py(s$seed, suffix = s$suffix, round = r)
      if (r < N_ROUNDS) simulate_round(s$suffix, round = r + 1)
    }
  }
  cat("NPE training done.\n")
}

# ---------------------------------------------------------------------------
# STEP mcmc — strict validation MCMC
# ---------------------------------------------------------------------------
step_mcmc <- function() {
  post_npe <- read.csv(file.path(OUT_DIR, "posterior_samples_npe.csv"))
  th_npe <- as.matrix(post_npe[, paste0("theta", 1:12)])
  Sigma0 <- cov(th_npe)
  # ensure positive definite
  Sigma0 <- Sigma0 + 1e-8 * diag(12)

  sol <- fit$solutions
  best_id <- fit$best$start_id
  start_ids <- c(best_id, best_id)
  starts <- list(
    clip_to_bounds(apply(th_npe, 2, median)),
    clip_to_bounds(as.numeric(sol[sol$start_id == start_ids[1], paste0("theta", 1:12)])),
    clip_to_bounds(as.numeric(sol[sol$start_id == start_ids[2], paste0("theta", 1:12)]))
  )
  start_names <- c("npe_median", start_ids)

  logpost <- function(th) log_posterior(th, base_params, data_local, priors,
                                        target_mode = TARGET_MODE)
  block <- N_ITER_MCMC
  burnin <- BURNIN_MCMC
  thin   <- THIN_MCMC
  chains <- NULL
  total_iter <- 0L
  repeat {
    cat("MCMC validation block: 3 chains x", block,
        "iters (total", total_iter + block, ")...\n")
    run_chain <- function(j) {
      init <- if (is.null(chains)) NULL else chains[[j]]$state
      am_mcmc(starts[[j]], logpost, n_iter = block, burnin = burnin,
              Sigma0 = Sigma0, seed = SEED + 10 * j + total_iter,
              thin = thin, init_state = init)
    }
    new <- parallel::mclapply(1:3, run_chain, mc.cores = N_CORES)
    if (is.null(chains)) {
      chains <- new
    } else {
      for (j in 1:3) {
        chains[[j]]$samples <- rbind(chains[[j]]$samples, new[[j]]$samples)
        chains[[j]]$logpost <- c(chains[[j]]$logpost, new[[j]]$logpost)
        chains[[j]]$iteration <- c(chains[[j]]$iteration, new[[j]]$iteration)
        chains[[j]]$state <- new[[j]]$state
        chains[[j]]$acceptance <-
          new[[j]]$state$n_accept / new[[j]]$state$n_iter
      }
    }
    total_iter <- total_iter + block
    burnin <- 0L
    names(chains) <- start_names
    acc <- sapply(chains, `[[`, "acceptance")
    cat("  acceptance:", paste(round(acc, 3), collapse = ", "), "\n")

    mcl <- coda::mcmc.list(lapply(chains, function(ch) {
      m <- ch$samples; colnames(m) <- paste0("theta", 1:12); coda::mcmc(m)
    }))
    rhat <- coda::gelman.diag(mcl, autoburnin = FALSE, multivariate = TRUE)
    ess <- do.call(rbind, lapply(chains, function(ch) {
      m <- ch$samples; colnames(m) <- paste0("theta", 1:12)
      coda::effectiveSize(coda::mcmc(m))
    }))
    ess_pool <- colSums(ess)
    rhat_ok <- all(rhat$psrf[, 1] >= 0.995 & rhat$psrf[, 1] <= 1.005)
    ess_ok  <- all(ess_pool > 1000)
    cat("  R-hat range:", round(range(rhat$psrf[, 1]), 4),
        "| ESS pooled range:", round(range(ess_pool), 1), "\n")
    if (rhat_ok && ess_ok) break
    cat("  Criteria not met; extending chains.\n")
    if (total_iter >= MCMC_MAX_ITER) {
      warning("MCMC max iterations reached without meeting strict criteria")
      break
    }
  }

  post_mcmc <- do.call(rbind, lapply(seq_along(chains), function(j) {
    ch <- chains[[j]]
    smp <- ch$samples; colnames(smp) <- paste0("theta", 1:12)
    data.frame(chain = start_names[j], iteration = ch$iteration,
               smp, logpost = ch$logpost)
  }))
  write.csv(post_mcmc, file.path(OUT_DIR, "posterior_samples_mcmc.csv"),
            row.names = FALSE)

  diag_df <- data.frame(
    parameter = paste0("theta", 1:12),
    rhat = rhat$psrf[, 1],
    ess_chain1 = ess[1, ], ess_chain2 = ess[2, ], ess_chain3 = ess[3, ],
    ess_pooled = ess_pool,
    acceptance_mean = mean(acc),
    n_iter = total_iter, burnin = BURNIN_MCMC, thin = thin,
    posterior_median = apply(do.call(rbind, mcl), 2, median),
    nm_theta = fit$best$theta
  )
  diag_df$mpsrf <- rhat$mpsrf
  write.csv(diag_df, file.path(OUT_DIR, "diagnostics_mcmc.csv"), row.names = FALSE)
  cat("MCMC validation done. Final total iterations:", total_iter, "\n")
  invisible(list(chains = chains, mcl = mcl, diag = diag_df))
}

# ---------------------------------------------------------------------------
# STEP predict — posterior predictive intervals (NPE + MCMC)
# ---------------------------------------------------------------------------
step_predict <- function() {
  lap <- fit$laplace$intervals
  npe <- read.csv(file.path(OUT_DIR, "posterior_samples_npe.csv"))
  mcmc <- read.csv(file.path(OUT_DIR, "posterior_samples_mcmc.csv"))

  # NPE: cap at 30000 draws
  th_npe <- as.matrix(npe[sample(seq_len(nrow(npe)), min(nrow(npe), 30000)),
                          paste0("theta", 1:12)])
  # MCMC: thin by 3 -> ~6000 draws total
  th_mcmc <- as.matrix(mcmc[seq(1, nrow(mcmc), by = 3),
                            paste0("theta", 1:12)])

  cat("Posterior predictive: NPE", nrow(th_npe), "draws, MCMC",
      nrow(th_mcmc), "draws...\n")
  pred_npe <- predictive_intervals(th_npe, base_params, data_local,
                                   target_mode = TARGET_MODE,
                                   n_cores = N_CORES, seed = SEED + 20)
  pred_mcmc <- predictive_intervals(th_mcmc, base_params, data_local,
                                    target_mode = TARGET_MODE,
                                    n_cores = N_CORES, seed = SEED + 21)
  pm_npe <- attr(pred_npe, "predictive_matrix")
  pm_mcmc <- attr(pred_mcmc, "predictive_matrix")

  target_posterior <- data.frame(
    source = "npe",
    setNames(as.data.frame(pm_npe$p), paste0("p", 1:6)),
    setNames(as.data.frame(pm_npe$N), paste0("N", 1:6))
  )
  target_posterior <- rbind(target_posterior, data.frame(
    source = "mcmc",
    setNames(as.data.frame(pm_mcmc$p), paste0("p", 1:6)),
    setNames(as.data.frame(pm_mcmc$N), paste0("N", 1:6))
  ))
  write.csv(target_posterior, file.path(OUT_DIR, "target_posterior.csv"),
            row.names = FALSE)

  ci <- data.frame(
    age_group = pred_npe$age_group,
    p_npe_median = pred_npe$p_median, p_npe_lo = pred_npe$p_lo,
    p_npe_hi = pred_npe$p_hi,
    p_mcmc_median = pred_mcmc$p_median, p_mcmc_lo = pred_mcmc$p_lo,
    p_mcmc_hi = pred_mcmc$p_hi,
    p_la_hat = lap$p_hat, p_la_lo = lap$p_lo, p_la_hi = lap$p_hi,
    p_obs_lo = lap$p_obs_lo, p_obs_hi = lap$p_obs_hi,
    N_npe_median = pred_npe$N_median, N_npe_lo = pred_npe$N_lo,
    N_npe_hi = pred_npe$N_hi,
    N_mcmc_median = pred_mcmc$N_median, N_mcmc_lo = pred_mcmc$N_lo,
    N_mcmc_hi = pred_mcmc$N_hi,
    N_la_hat = lap$N_hat, N_la_lo = lap$N_lo, N_la_hi = lap$N_hi,
    N_obs_lo = lap$N_obs_lo, N_obs_hi = lap$N_obs_hi,
    p_overlap_npe_obs = !(pred_npe$p_hi < lap$p_obs_lo |
                          pred_npe$p_lo > lap$p_obs_hi),
    p_overlap_mcmc_obs = !(pred_mcmc$p_hi < lap$p_obs_lo |
                           pred_mcmc$p_lo > lap$p_obs_hi),
    p_overlap_npe_mcmc = !(pred_npe$p_hi < pred_mcmc$p_lo |
                           pred_npe$p_lo > pred_mcmc$p_hi),
    N_overlap_npe_obs = !(pred_npe$N_hi < lap$N_obs_lo |
                          pred_npe$N_lo > lap$N_obs_hi),
    N_overlap_mcmc_obs = !(pred_mcmc$N_hi < lap$N_obs_lo |
                           pred_mcmc$N_lo > lap$N_obs_hi),
    N_overlap_npe_mcmc = !(pred_npe$N_hi < pred_mcmc$N_lo |
                           pred_npe$N_lo > pred_mcmc$N_hi),
    n_npe_draws = pred_npe$n_draws_used[1],
    n_mcmc_draws = pred_mcmc$n_draws_used[1]
  )
  write.csv(ci, file.path(OUT_DIR, "credible_intervals.csv"), row.names = FALSE)
  cat("Credible intervals written.\n")
}

# ---------------------------------------------------------------------------
# STEP plots — ggplot2
# ---------------------------------------------------------------------------
step_plots <- function() {
  mcmc <- read.csv(file.path(OUT_DIR, "posterior_samples_mcmc.csv"))
  npe  <- read.csv(file.path(OUT_DIR, "posterior_samples_npe.csv"))
  npe2 <- read.csv(file.path(OUT_DIR, "posterior_samples_npe_seed2.csv"))
  diag <- read.csv(file.path(OUT_DIR, "diagnostics_mcmc.csv"))
  sbc  <- read.csv(file.path(OUT_DIR, "sbc_summary.csv"))
  tp   <- read.csv(file.path(OUT_DIR, "target_posterior.csv"))
  ci   <- read.csv(file.path(OUT_DIR, "credible_intervals.csv"))

  chains_list <- lapply(split(mcmc, mcmc$chain), function(d) {
    m <- as.matrix(d[, paste0("theta", 1:12)]); colnames(m) <- paste0("theta", 1:12); m
  })
  plot_trace(chains_list, file.path(OUT_DIR, "fig_trace.png"))

  # density: NPE + MCMC + prior
  plot_density_npe(npe, npe2, mcmc, priors, file.path(OUT_DIR, "fig_density.png"))

  tp_npe <- tp[tp$source == "npe", ]
  plot_predictive_density(tp_npe, fit, file.path(OUT_DIR,
                                                 "fig_predictive_density.png"))
  plot_ci_compare_npe(ci, fit, file.path(OUT_DIR, "fig_ci_compare.png"))
  plot_sbc(sbc, file.path(OUT_DIR, "fig_sbc.png"))

  cfg <- data.frame(
    field = c("run_id", "fit_rds", "git_sha", "setup_R_md5", "sim_cpp_md5",
              "pptx_md5", "R_version", "coda_version",
              "npe_python", "sbi_version", "torch_version",
              "prior", "n_train_sims", "n_rounds", "n_proposal_per_round",
              "npe_draws", "sbc_n", "sbc_draws_per_x",
              "mcmc_chains", "mcmc_total_iter", "mcmc_burnin", "mcmc_thin",
              "mcmc_rhat_range", "mcmc_ess_pooled_range", "mcmc_acceptance",
              "primary_method", "fallback_reason"),
    value = c(
      basename(OUT_DIR), FIT_PATH,
      system("git rev-parse HEAD", intern = TRUE),
      unname(tools::md5sum(file.path(ROOT, "src", "setup.R"))),
      unname(tools::md5sum(file.path(ROOT, "src", "sim.cpp"))),
      unname(tools::md5sum(file.path(ROOT, "Model schematic.pptx"))),
      R.version.string, as.character(packageVersion("coda")),
      PYTHON,
      system(paste(shQuote(PYTHON), "-c",
                   "'import sbi; print(sbi.__version__)'"),
             intern = TRUE),
      system(paste(shQuote(PYTHON), "-c",
                   "'import torch; print(torch.__version__)'"),
             intern = TRUE),
      "log contact ~ Normal(0,2^2); log beta ~ Student-t(3,0,2)",
      as.character(N_SIMS), as.character(N_ROUNDS),
      as.character(N_PROPOSAL), as.character(N_DRAWS),
      "298", "1000",
      "3", as.character(unique(diag$n_iter)), as.character(BURNIN_MCMC),
      "5",
      paste0("[", min(diag$rhat), ", ", max(diag$rhat), "]"),
      paste0("[", min(diag$ess_pooled), ", ", max(diag$ess_pooled), "]"),
      paste0("[", min(diag$acceptance_mean), ", ",
             max(diag$acceptance_mean), "]"),
      "MCMC (adaptive Metropolis) after NPE attempted",
      paste("NPE overconfident in well-identified directions and unstable",
            "across seeds in weakly identified directions",
            "(see bayes_methodology.md)")
    ),
    stringsAsFactors = FALSE
  )
  write.csv(cfg, file.path(OUT_DIR, "run_config.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()),
             file.path(OUT_DIR, "sessionInfo.txt"))
  cat("Figures written.\n")
}

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
if (STEP %in% c("data", "all")) step_data()
if (STEP %in% c("train", "all")) step_train()
if (STEP %in% c("mcmc", "all")) step_mcmc()
if (STEP %in% c("predict", "all")) step_predict()
if (STEP %in% c("plots", "all")) step_plots()
cat("Done.\n")
