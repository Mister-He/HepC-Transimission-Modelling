# =============================================================================
# setup.R  —  Parameter and data lists for the HCV PWID compartmental model
# Compile the C++ ODE solver first:
#   Rcpp::sourceCpp("sim.cpp")
# Then source this file and call run_sim(params, data)
# =============================================================================

library(Rcpp)
library(dplyr)
sourceCpp("sim.cpp")

# =============================================================================
# HELPER: flat compartment index (mirrors C++ idx())
# s in {0=D,1=J,2=F,3=X}, k in {1,2,3,4}, h in {0=u,1=a,2=c,3=t}, i in {0..8}
# =============================================================================
idx <- function(s, k, h, i) s * 4 * 4 * 9 + (k - 1) * 4 * 9 + h * 9 + i + 1L

# =============================================================================
# PARAMETERS
# Sources:
#   q         — Boelen et al. (2014)
#   kappa     — Micallef et al. (2006)
#   iota1/2   — Fasano (2024); WHO (2026)
#   rho       — Zhao et al. (2019)
#   alpha_*   — Pawlotsky (2020 EASL); Xiao et al. (2025)
#   p_*       — Lim et al. (2018)
#   r3_*      — Kanwal et al. (2014); derived (see Section X)
#   phi_*     — van Meer (2012); Morgan et al. (2013)
#   mu_DC/HCC — Lim et al. (2018)
#   psi_*     — Premkumar (2024); Dang et al. (2020)
#   pi_recid  — fitted to SPS data (0.65)
#   lambda*   — CALIBRATED to SPS data (guesses used here)
#   mu        — CALIBRATED to SingStat life tables (guesses used here)
#   c_contact — CALIBRATED to HCV prevalence trajectory (guess used here)
#   beta      — CALIBRATED to PWID population size (guess used here)
# =============================================================================

params <- list(

  # ── Transmission & natural history ────────────────────────────────────────
  q = 0.0057, # transmission prob per needle-sharing event
  kappa = 0.26, # prob spontaneous clearance after acute infection
  iota1 = 180 / 365, # mean acute duration (years); 180 days
  iota2 = 84 / 365, # mean treatment duration (years); 84 days = 12 weeks

  # ── Genotype ──────────────────────────────────────────────────────────────
  rho = 0.78, # proportion with genotype 3

  # ── SVR12 rates (effective, genotype-weighted) ─────────────────────────────
  # Stages 1-2 (NC & CC), RBV-free
  # alpha_NC = rho * alpha3_NC + (1-rho) * alphao_NC
  #          = 0.78*0.976 + 0.22*0.988 = 0.9787
  alpha_NC = 0.9787,

  # Stage 3 (DC) with RBV+
  # alpha_DC_pos = 0.78*0.87 + 0.22*0.92 = 0.8818
  alpha_DC_pos = 0.8818,

  # Stage 3 (DC) without RBV (used in sensitivity)
  # alpha_DC_neg = 0.78*0.45 + 0.22*0.90 = 0.5490
  alpha_DC_neg = 0.5490,

  # Stage 4 (HCC): assumed equal to DC+ (no dedicated trial data; see Section X)
  alpha_HCC = 0.8818,

  # ── Treatment initiation rates (per year, per stage) ──────────────────────
  # GUESS: set to 0 (no treatment) as base scenario; override in scenarios
  tau = c(
    0.0, # tau_1: NC   — SCENARIO-DEFINED; 0 = no treatment
    0.0, # tau_2: CC   — SCENARIO-DEFINED
    0.0, # tau_3: DC   — SCENARIO-DEFINED
    0.0 # tau_4: HCC  — SCENARIO-DEFINED
  ),

  # ── Baseline progression rates (per year, other genotype) ─────────────────
  p_NC_CC = 0.027, # NC  → CC  (Thein et al. 2008)
  p_CC_DC = 0.039, # CC  → DC  (Lim et al. 2018)
  p_CC_HCC = 0.014, # CC  → HCC (Lim et al. 2018)
  p_DC_HCC = 0.014, # DC  → HCC (Lim et al. 2018)

  # ── GT3 relative risks ────────────────────────────────────────────────────
  r3_NC_CC = 1.36, # r3^{NC→CC}  derived; CI: 1.27–1.46
  r3_CC_DC = 1.36, # r3^{CC→DC}  derived; CI: 1.27–1.46
  r3_CC_HCC = 1.93, # r3^{CC→HCC} derived; CI: 1.71–2.17
  r3_DC_HCC = 1.93, # r3^{DC→HCC} derived; CI: 1.71–2.17

  # ── Post-SVR progression modifiers ────────────────────────────────────────
  phi_CC_DC = 0.07, # relative risk CC→DC after SVR (van Meer 2012)
  phi_CC_HCC = 0.23, # relative risk CC→HCC after SVR (Morgan et al. 2013)
  phi_DC_HCC = 1.00, # assumed no reduction DC→HCC after SVR

  # ── Background mortality (per year, age-varying) ──────────────────────────
  # CALIBRATED: placeholder values from SingStat life table
  # Age groups (example boundaries): 15-19, 20-24, 25-29, 30-34,
  #                                   35-39, 40-44, 45-49, 50-54, 55+
  mu = c(
    0.001267, # age group 1 (15-19)
    0.000300, # age group 2 (20-24)
    0.000300, # age group 3 (25-29)
    0.000400, # age group 4 (30-34)
    0.000500, # age group 5 (35-39)
    0.000700, # age group 6 (40-44)
    0.001400, # age group 7 (45-49)
    0.002300, # age group 8 (50-54)
    0.016100 # age group 9 (55+)
  ),

  # ── Standardized Mortality rate of ever-PWIDs  ─────────────────────────────
  omega = 14.68, # SMR for ever-PWIDs (Degenhardt et al. 2011)

  # ── Disease-specific excess mortality ──────────────────────────────────────
  mu_DC = 0.130, # additional mortality in decompensated cirrhosis (Lim 2018)
  mu_HCC = 0.430, # additional mortality in HCC                     (Lim 2018)

  # ── SVR mortality modifiers ────────────────────────────────────────────────
  psi_DC = 0.45, # relative mortality risk DC after SVR  (Premkumar 2024)
  psi_HCC = 0.37, # relative mortality risk HCC after SVR (Dang 2020)

  # ── Incarceration rates (per year, age-varying) ────────────────────────────
  # CALIBRATED: placeholder values — replace with SPS-fitted rates
  lambda1 = c(
    0.9769631, 0.7962602, 0.7842585,
    0.8721055, 0.8675996, 0.9522804,
    1.0377386, 0.9348018, 0.9021672
  ), # first-arrest rate   lambda_i^(1) — GUESS
  lambda2 = c(
    0.489, 0.620, 0.663,
    0.628, 0.533, 0.475,
    0.472, 0.441, 0.451
  ), # release rate        lambda_i^(2) — GUESS (0.5yr avg)
  lambda3 = c(
    0.9769631, 0.7962602, 0.7842585,
    0.8721055, 0.8675996, 0.9522804,
    1.0377386, 0.9348018, 0.9021672
  ), # re-arrest rate      lambda_i^(3) — GUESS
  pi_recid = 0.61610608, # recidivism probability (fitted to SPS; Assumption)

  # ── Needle-sharing contact rate ────────────────────────────────────────────
  # CALIBRATED: scalar homogeneous mixing — replace with 9×9 matrix post-calib.
  C_contact = rbind(
    c(7, 4, 1, 1, 0, 1, 1, 0, 0),
    c(11, 34, 21, 11, 6, 2, 1, 0.5, 0.5),
    c(7, 30, 80, 62, 30, 10, 2, 1.5, 1.5),
    c(2, 10, 60, 121, 65, 38, 15, 2, 2),
    c(1, 11, 22, 67, 107, 41, 18, 2.5, 2.5),
    c(0, 4, 6, 22, 32, 31, 10, 2, 2),
    c(0, 1, 1, 8, 10, 15, 11, 0.5, 0.5),
    c(0, 5, 5, 5.5, 6.5, 6.5, 3.5, 1.5, 1.5),
    c(0, 5, 5, 5.5, 6.5, 6.5, 3.5, 1.5, 1.5)
  ) * 4,

  # ── Population entry rates (per year, age-varying) ────────────────────────
  # CALIBRATED: constant-in-time placeholder — replace with beta_i(t) from calib.
  beta = c(
    235, # age group 1
    565 / 2, # age group 2
    565 / 2, # age group 3
    301 / 2, # age group 4
    301 / 2, # age group 5
    111 / 2, # age group 6
    111 / 2, # age group 7
    33 / 2, # age group 8
    33 / 2 + 4 # age group 9
  )
)

# =============================================================================
# INITIAL CONDITIONS
# 576 compartments initialised to near-zero.
# A small seed of chronic infection is placed in D_c,1,i (never-incarcerated,
# non-cirrhosis, chronic, all age groups) to start the epidemic.
# CALIBRATED: replace with equilibrium-derived or SPS/CNB baseline estimates.
# =============================================================================
y0 <- rep(0.0, 576)

# Susceptible population: put most of PWID in D_u,1,i
# Seed chronic infection: 20 chronic per age group in D_c,1,i
pos <- c(55, 145, 183, 164, 212, 299, 222, 190, 133)
tot <- c(307, 797, 829, 633, 598, 642, 481, 439, 366)
for (i in 0:8) {
  y0[idx(s = 0, k = 1, h = 0, i = i)] <- tot[i + 1] - pos[i + 1] # D_{u,1,i}
  y0[idx(s = 0, k = 1, h = 1, i = i)] <- pos[i + 1] # D_{a,1,i}
}

# =============================================================================
# DATA / SIMULATION SETTINGS
# =============================================================================
# Observation-model dispersion settings.
# These can be tuned or moved into `data` if you want to estimate them.
sigma_N <- 0.10
phi_overdisp <- 50.0

data <- list(
  t_start = 0.0, # start year (0 = model year 0; map to calendar year in R)
  t_end   = 100.0, # simulate 100 years
  dt      = 1 / 365, # daily time steps (1/365 of a year)
  y0      = y0, # initial conditions (length-576 vector)
  sigma_N = sigma_N,
  phi_overdisp = phi_overdisp
)

# =============================================================================
# SCENARIO HELPERS
# =============================================================================

# Scenario 1: no treatment (baseline counterfactual)
params_s1 <- modifyList(params, list(tau = c(0, 0, 0, 0)))

# # Scenario 2: treat NC and CC only (limited access)
# params_s2 <- modifyList(params, list(tau = c(0.5, 0.5, 0, 0)))

# # Scenario 3: treat all stages
# params_s3 <- modifyList(params, list(tau = c(0.5, 0.5, 0.3, 0.2)))

# # Scenario 4: aggressive elimination target
# params_s4 <- modifyList(params, list(tau = c(0.8, 0.8, 0.6, 0.4)))

# # Scenario 5: sensitivity analysis with reduced DC SVR efficacy
# params_s5 <- modifyList(params, list(tau = c(1.0, 1.0, 1.0, 1.0)))

# =============================================================================
# COLUMN NAME HELPER (for labelling the output matrix)
# =============================================================================
strata_names <- c("D", "J", "F", "X")
stage_names <- c("NC", "CC", "DC", "HCC")
state_names <- c("u", "a", "c", "t")
age_names <- paste0("age", 1:9)

expand.grid(age_names, state_names, stage_names, strata_names) %>%
  mutate(col_name = paste(Var4, Var3, Var2, Var1, sep = "_")) %>%
  pull(col_name) %>%
  c("time", .) -> col_names

# =============================================================================
# EXAMPLE RUN
# =============================================================================
start <- Sys.time()
out <- run_sim(params_s1, data)
end <- Sys.time()


# =============================================================================
# hmc_calibration.R  —  HMC Calibration for the HCV PWID Compartmental Model
#
# Append this file after setup.R (which sources sim.cpp and defines params,
# data, idx(), y0, etc.).
#
# Parameters estimated (12 total, all in log-space for positivity):
#   theta[1]       = log(beta_scale)        scalar multiplier on beta[1:9]
#   theta[2]       = log(delta)             shift parameter for shifted Gamma prior
#   theta[3]       = log(alpha)             shape parameter for Gamma hyperprior
#   theta[4]       = log(beta_rate)         rate parameter for Gamma hyperprior
#   theta[5:12]    = log(C_contact_scale[1:8])  8 age-specific row scalings
#
# Scaling factors generated:
#   C_contact_scale[1:8] = delta + Gamma(alpha, beta_rate)
#   C_contact_scale[9]    = 1  (reference row, fixed for identifiability)
#
# Likelihood (three-step observation model, J-stratum only, per age group i = 1..9):
#   N_total_obs            ~ LogNormal(log(N_total_model), sigma_N)
#   obs_tot[1:9] | N_total ~ Multinomial(N_total_obs, p_1:9(theta))
#   obs_pos[a] | obs_tot[a] ~ BetaBinomial(obs_tot[a], q_a(theta), phi)
#
# Prior: LogNormal hyperpriors on (beta_scale, delta, alpha, beta_rate)
#   beta_scale ~ LogNormal(0, 1)
#   delta      ~ LogNormal(log(0.25), 0.5)   [shift/floor for row scalings]
#   alpha      ~ LogNormal(log(2.5), 0.3)    [Gamma shape hyperparameter]
#   beta_rate  ~ LogNormal(log(1.0), 0.3)    [Gamma rate hyperparameter]
#
# The 8 free row scalings follow a shifted Gamma hierarchy, with the 9th row
# fixed to 1 as the reference row for identifiability.
#
# Gradient: central finite differences through the C++ ODE solver
#   Cost per HMC step: (L+1) gradient evals × 2 × N_PARAMS ODE runs
#                    = (L+1) × 24 ODE runs  [plus 1 for current lp]
#   Tip: keep L small (5-10) and use parallel chains.
#
# Outputs:
#   hmc_chains        — raw chain list (all iterations, all chains)
#   post_warmup_list  — post-warmup samples per chain
#   all_samples       — pooled post-warmup matrix (n_samp × 12)
#   diag_table        — R-hat and ESS summary data.frame
#   ppc_out           — posterior predictive replicates
#   plots             — trace, PPC histogram, PPC interval
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# =============================================================================
# 0.  OBSERVATIONS
# =============================================================================

obs_pos <- c(3, 8, 6, 5, 13, 31, 32, 36, 37) * 5
obs_tot <- c(50, 85, 69, 38, 42, 91, 71, 78, 100) * 5
N_total_obs <- sum(obs_tot)

stopifnot(length(obs_pos) == 9L, length(obs_tot) == 9L)

N_PARAMS         <- 12L
param_names_log  <- c(
  "log_beta_scale", "log_delta", "log_alpha", "log_beta_rate",
  paste0("log_C_contact_scale_", 1:8)
)
param_names_orig <- c(
  "beta_scale", "delta", "alpha", "beta_rate",
  paste0("C_contact_scale_", 1:8)
)


# =============================================================================
# 1.  PARAMETER TRANSFORMS
# =============================================================================

#' Unconstrained theta (log-space) -> constrained named list
#' 
#' Generates 8 free contact scaling factors from a shifted Gamma hierarchy.
#' The 9th row is held at 1 as a reference for identifiability.
constrain_theta <- function(theta) {
  beta_scale <- exp(theta[1L])
  delta      <- exp(theta[2L])
  alpha      <- exp(theta[3L])
  beta_rate  <- exp(theta[4L])

  # 8 estimated row scalings; row 9 remains the reference row (scale = 1)
  c_scales <- c(delta + exp(theta[5:12]), 1.0)
  
  list(
    beta_scale      = beta_scale,
    C_contact_scale = c_scales,
    delta           = delta,
    alpha           = alpha,
    beta_rate       = beta_rate
  )
}

#' Build a full parameter list from unconstrained theta
#'
#' The 8 free row scaling factors are hierarchical draws on the log scale,
#' while the 9th row remains fixed at 1.
build_params_from_theta <- function(theta, base_params) {
  p  <- constrain_theta(theta)
  pm <- base_params

  # Scale beta uniformly across all age groups
  pm$beta <- base_params$beta * p$beta_scale

  # Scale rows 1:8; row 9 is the reference row and stays unchanged.
  for (j in 1:8) {
    pm$C_contact[j, ] <- base_params$C_contact[j, ] * p$C_contact_scale[j]
  }
  pm
}


# =============================================================================
# 2.  POISSON MEANS FROM FINAL ODE STATE
# =============================================================================

log_beta_binomial <- function(x, size, prob, phi) {
  if (length(x) != length(size) || length(x) != length(prob)) {
    stop("x, size, and prob must have the same length")
  }
  if (!is.finite(phi) || phi <= 0) return(-Inf)
  if (any(!is.finite(x)) || any(!is.finite(size)) || any(!is.finite(prob))) return(-Inf)
  if (any(x < 0) || any(size < 0) || any(x > size)) return(-Inf)

  prob <- pmin(pmax(prob, 1e-12), 1 - 1e-12)
  a <- prob * phi
  b <- (1 - prob) * phi

  sum(lchoose(size, x) + lbeta(x + a, size - x + b) - lbeta(a, b))
}

#' Compute model-implied age-group totals and HCV-positive counts from the
#' J-stratum (s = 1) at the final ODE state.
compute_age_quantities <- function(y_final) {
  age_total <- numeric(9L)
  age_pos <- numeric(9L)

  for (i in 0:8) {
    total_i <- 0.0
    pos_i <- 0.0
    for (k in 1:4) {
      total_i <- total_i + y_final[idx(s = 1L, k = k, h = 0L, i = i)]
      total_i <- total_i + y_final[idx(s = 1L, k = k, h = 1L, i = i)]
      total_i <- total_i + y_final[idx(s = 1L, k = k, h = 2L, i = i)]
      total_i <- total_i + y_final[idx(s = 1L, k = k, h = 3L, i = i)]

      # Preserve the current positive-state definition used in the model.
      pos_i <- pos_i + y_final[idx(s = 1L, k = k, h = 0L, i = i)]
      pos_i <- pos_i + y_final[idx(s = 1L, k = k, h = 2L, i = i)]
    }
    age_total[i + 1L] <- total_i
    age_pos[i + 1L] <- pos_i
  }

  n_model_total <- sum(age_total)
  if (!is.finite(n_model_total) || n_model_total <= 0) return(NULL)
  if (any(!is.finite(age_total)) || any(!is.finite(age_pos))) return(NULL)
  if (any(age_total <= 0)) return(NULL)
  if (any(age_pos < 0) || any(age_pos > age_total)) return(NULL)

  list(
    total_by_age = age_total,
    pos_by_age   = age_pos,
    n_model_total = n_model_total,
    p_age = age_total / n_model_total,
    q_age = pmin(pmax(age_pos / age_total, 1e-12), 1 - 1e-12)
  )
}


# =============================================================================
# 3.  LOG-POSTERIOR AND GRADIENT
# =============================================================================

# ── 3a. Log-prior ────────────────────────────────────────────────────────────
#
# Prior specification on original scale:
#   beta_scale ~ LogNormal(0, 1)                 [median 1]
#   delta      ~ LogNormal(log(0.25), 0.5)       [positive shift/floor]
#   alpha      ~ LogNormal(log(2.5), 0.3)        [Gamma shape hyperparameter]
#   beta_rate  ~ LogNormal(log(1.0), 0.3)        [Gamma rate hyperparameter]
#
# The 8 free row scalings are modeled as delta + Gamma(alpha, beta_rate),
# which places the prior mass for the resulting scales around 1–2.
#
# Change-of-variables: theta_j = log(X_j), X_j = exp(theta_j)
#   log p(theta_j) = log p_X(exp(theta_j)) + theta_j
#
log_prior <- function(theta) {
  beta_scale <- exp(theta[1L])
  delta      <- exp(theta[2L])
  alpha      <- exp(theta[3L])
  beta_rate  <- exp(theta[4L])
  excess     <- exp(theta[5:12])
  
  # LogNormal priors on each parameter (using dnorm on log scale + Jacobian)
  lp_beta_scale <- dnorm(theta[1L], mean =  0.0, sd = 1.0, log = TRUE)
  lp_delta      <- dnorm(theta[2L], mean = log(0.25), sd = 0.5, log = TRUE)
  lp_alpha      <- dnorm(theta[3L], mean = log(2.5), sd = 0.3, log = TRUE)
  lp_beta_rate  <- dnorm(theta[4L], mean = log(1.0), sd = 0.3, log = TRUE)
  lp_contact    <- sum(dgamma(excess, shape = alpha, rate = beta_rate, log = TRUE) + theta[5:12])
  
  # Sum across all parameters
  lp_beta_scale + lp_delta + lp_alpha + lp_beta_rate + lp_contact
}

# ── 3b. Analytical gradient of log-prior ─────────────────────────────────────
#
# With shifted Gamma priors on the row scalings:
#   d/d theta_j [dgamma(exp(theta_j), shape, rate, log=TRUE) + theta_j]
#     = alpha - beta_rate * exp(theta_j)
#
grad_log_prior_analytical <- function(theta) {
  grad <- numeric(12L)
  grad[1L] <- -(theta[1L] - 0.0) / 1.0^2
  grad[2L] <- -(theta[2L] - log(0.25)) / 0.5^2
  grad[3L] <- -(theta[3L] - log(2.5)) / 0.3^2
  grad[4L] <- -(theta[4L] - log(1.0)) / 0.3^2
  alpha     <- exp(theta[3L])
  beta_rate <- exp(theta[4L])
  grad[5:12] <- alpha - beta_rate * exp(theta[5:12])
  grad
}

# ── 3c. Log-likelihood ───────────────────────────────────────────────────────
log_likelihood <- function(theta, base_params, data) {
  pm  <- build_params_from_theta(theta, base_params)
  out <- run_sim(pm, data)

  if (!is.matrix(out) || nrow(out) == 0L) return(-Inf)
  if (any(!is.finite(out)))               return(-Inf)

  y_final <- as.numeric(out[nrow(out), -1L])   # drop time column
  obs     <- compute_age_quantities(y_final)
  if (is.null(obs)) return(-Inf)

  sigma_N <- if (!is.null(data$sigma_N)) data$sigma_N else 0.10
  phi     <- if (!is.null(data$phi_overdisp)) data$phi_overdisp else 50.0

  ll_total <- dlnorm(N_total_obs, meanlog = log(obs$n_model_total), sdlog = sigma_N, log = TRUE)
  ll_age   <- dmultinom(obs_tot, prob = obs$p_age, log = TRUE)
  ll_case  <- log_beta_binomial(obs_pos, size = obs_tot, prob = obs$q_age, phi = phi)

  ll_total + ll_age + ll_case
}

# ── 3d. Log-posterior ────────────────────────────────────────────────────────
log_posterior <- function(theta, base_params, data) {
  lp <- log_prior(theta)
  if (!is.finite(lp)) return(-Inf)

  ll <- tryCatch(
    log_likelihood(theta, base_params, data),
    error = function(e) -Inf
  )
  if (!is.finite(ll)) return(-Inf)
  lp + ll
}

# ── 3e. Numerical gradient (central finite differences) ───────────────────────
#
# Analytical gradients are not available because the likelihood is computed by
# a black-box C++ ODE solver. We therefore use central finite differences:
#
#   ∂ log π / ∂ theta_j  ≈
#       [log π(theta + eps·e_j) - log π(theta - eps·e_j)] / (2·eps)
#
# This requires 2 × N_PARAMS = 20 ODE runs per gradient evaluation.
# With L leapfrog steps, each HMC proposal costs ≈ (L + 1) × 20 ODE runs.
#
# Fallback: if either perturbed evaluation is -Inf, the component is set to 0
# (safe — the chain will reject or stay put, not diverge).
#
grad_log_posterior <- function(theta, base_params, data, eps = 1e-4) {
  n    <- length(theta)
  grad <- numeric(n)
  lp0  <- NULL   # lazily evaluated if needed for fallback

  for (j in seq_len(n)) {
    theta_p    <- theta; theta_p[j] <- theta[j] + eps
    theta_m    <- theta; theta_m[j] <- theta[j] - eps

    lp_p <- log_posterior(theta_p, base_params, data)
    lp_m <- log_posterior(theta_m, base_params, data)

    if (is.finite(lp_p) && is.finite(lp_m)) {
      grad[j] <- (lp_p - lp_m) / (2 * eps)
    } else if (is.finite(lp_p)) {
      if (is.null(lp0)) lp0 <- log_posterior(theta, base_params, data)
      grad[j] <- (lp_p - lp0) / eps
    } else if (is.finite(lp_m)) {
      if (is.null(lp0)) lp0 <- log_posterior(theta, base_params, data)
      grad[j] <- (lp0 - lp_m) / eps
    } else {
      grad[j] <- 0.0   # both sides infeasible — zero gradient (safe)
    }

  }
  grad
}

# =============================================================================
# 4.  LEAPFROG INTEGRATOR
# =============================================================================
#
# Störmer–Verlet / leapfrog scheme for Hamiltonian H(θ, r) = U(θ) + K(r):
#   U(θ) = -log π(θ)   (potential energy)
#   K(r) = ½ rᵀr       (kinetic energy, unit mass matrix)
#
# Algorithm:
#   r_{½}     ← r₀ + (ε/2) ∇log π(θ₀)          half-step momentum
#   θ_l       ← θ_{l-1} + ε · r_{l-½}            full position step  } ×L
#   r_{l+½}   ← r_{l-½} + ε · ∇log π(θ_l)       full momentum step  }
#   r_L       ← r_{L-½} + (ε/2) ∇log π(θ_L)     final half-step
#
# Total gradient evaluations: L + 1
#
leapfrog <- function(theta, r, grad_fn, eps, L) {
  g <- grad_fn(theta)             # ∇log π(θ₀)
  r <- r + 0.5 * eps * g         # half-step momentum

  for (l in seq_len(L)) {
    theta <- theta + eps * r      # full position update
    g     <- grad_fn(theta)       # ∇log π(θ_l)
    if (l < L) {
      r <- r + eps * g            # full momentum update (not last step)
    }
  }
  r <- r + 0.5 * eps * g         # final half-step (g = ∇log π(θ_L))

  list(theta = theta, r = r)
}


# =============================================================================
# 5.  HMC CHAIN WITH DUAL-AVERAGING STEP-SIZE ADAPTATION
# =============================================================================
#
# Step-size adaptation follows Hoffman & Gelman (2014) Algorithm 5:
#   H̄_m  ← (1 - 1/(m+t₀)) H̄_{m-1} + (δ - α_m)/(m+t₀)
#   ε_m   ← exp(μ - √m/γ · H̄_m)
#   ε̄_m  ← exp(m^{-κ} log ε_m + (1-m^{-κ}) log ε̄_{m-1})
# After warmup: fix ε = ε̄_{n_warmup}.
#
run_hmc_chain <- function(theta_init,
                           n_iter,
                           n_warmup,
                           eps_init    = 0.01,
                           L           = 10L,
                           adapt_delta = 0.65,
                           base_params,
                           data,
                           seed        = NULL,
                           chain_id    = 1L) {
  if (!is.null(seed)) set.seed(seed)

  n_p      <- length(theta_init)
  samples  <- matrix(NA_real_, n_iter, n_p,
                     dimnames = list(NULL, param_names_log))
  accepted  <- logical(n_iter)
  lp_trace  <- numeric(n_iter)
  eps_trace <- numeric(n_iter)

  theta <- theta_init
  lp    <- log_posterior(theta, base_params, data)

  # Dual-averaging state
  eps     <- eps_init
  eps_bar <- eps_init
  H_bar   <- 0.0
  mu      <- log(10 * eps_init)   # μ = log(10ε₀)
  gamma   <- 0.05
  t0      <- 10
  kappa   <- 0.75

  cat(sprintf(
    "[Chain %d] Starting  n_iter=%d | n_warmup=%d | eps0=%.4f | L=%d\n",
    chain_id, n_iter, n_warmup, eps_init, L))

  for (iter in seq_len(n_iter)) {
    # Sample momenta from N(0, I)
    r0 <- rnorm(n_p)
    K0 <- 0.5 * sum(r0^2)

    # Leapfrog proposal
    grad_fn <- function(t) grad_log_posterior(t, base_params, data)
    prop    <- leapfrog(theta, r0, grad_fn, eps, L)

    lp_prop <- log_posterior(prop$theta, base_params, data)
    K_prop  <- 0.5 * sum(prop$r^2)

    # Metropolis-Hastings acceptance
    log_alpha <- (lp_prop - K_prop) - (lp - K0)
    alpha     <- if (is.finite(log_alpha)) min(1.0, exp(log_alpha)) else 0.0

    if (runif(1L) < alpha) {
      theta <- prop$theta
      lp    <- lp_prop
      accepted[iter] <- TRUE
    }

    samples[iter, ]  <- theta
    lp_trace[iter]   <- lp
    eps_trace[iter]  <- eps

    # ── Dual-averaging adaptation (warmup only) ────────────────────────────
    if (iter <= n_warmup) {
      H_bar   <- (1 - 1 / (iter + t0)) * H_bar +
                 (adapt_delta - alpha) / (iter + t0)
      log_eps <- mu - sqrt(iter) / gamma * H_bar
      eps     <- exp(log_eps)
      eps_bar <- exp(iter^(-kappa) * log_eps +
                     (1 - iter^(-kappa)) * log(eps_bar))

      if (iter == n_warmup) {
        eps <- eps_bar
        cat(sprintf("[Chain %d] Warmup done. Adapted eps = %.5f\n",
                    chain_id, eps))
      }
    }

    # Progress report every 100 iterations
    if (iter %% 100L == 0L) {
      phase   <- if (iter <= n_warmup) "warmup" else "sample"
      win_acc <- mean(accepted[max(1L, iter - 99L):iter])
      cat(sprintf("[Chain %d] iter=%4d (%s) | lp=%8.2f | acc=%.2f | eps=%.5f\n",
                  chain_id, iter, phase, lp, win_acc, eps))
    }
  }
  
  idx_post <- seq(n_warmup + 1L, n_iter)
  list(
    samples     = samples[idx_post, , drop = FALSE],  # post-warmup only
    samples_all = samples,                             # full trace (for plots)
    lp_trace    = lp_trace,
    eps_trace   = eps_trace,
    accepted    = accepted,
    acc_rate    = mean(accepted[idx_post]),
    eps_final   = eps,
    n_warmup    = n_warmup
  )
}


# =============================================================================
# 6.  CONVERGENCE DIAGNOSTICS: R-HAT AND ESS
# =============================================================================

# ── 6a. Gelman-Rubin R-hat ───────────────────────────────────────────────────
#
# R-hat measures mixing across chains.  R-hat < 1.01 indicates convergence
# (Vehtari et al. 2021); values 1.01–1.05 are borderline.
#
# Formula (split-chain pooled variance):
#   B = (N/(M-1)) Σ_m (ψ̄_m - ψ̄)²    between-chain variance
#   W = (1/M)     Σ_m s²_m             within-chain variance
#   var̂(ψ) = (N-1)/N · W + B/N
#   R-hat   = sqrt(var̂(ψ) / W)
#
compute_rhat <- function(chains) {
  M <- length(chains)
  N <- nrow(chains[[1L]])
  J <- ncol(chains[[1L]])

  rhat <- setNames(numeric(J), param_names_log)

  for (j in seq_len(J)) {
    chain_j     <- lapply(chains, function(ch) ch[, j])
    chain_means <- vapply(chain_j, mean, numeric(1L))
    chain_vars  <- vapply(chain_j, var,  numeric(1L))

    grand_mean <- mean(chain_means)
    B          <- N / (M - 1) * sum((chain_means - grand_mean)^2)
    W          <- mean(chain_vars)
    var_hat    <- (N - 1) / N * W + B / N

    rhat[j]    <- sqrt(var_hat / W)
  }
  rhat
}

# ── 6b. Effective Sample Size (ESS) ──────────────────────────────────────────
#
# ESS based on Geyer's initial monotone sequence estimator applied to the
# pooled autocorrelation function across all chains.
#
# ESS = (M × N) / τ,   where τ = 1 + 2 Σ_{t=1}^{T*} ρ_t
# T* = last lag where ρ_{2t} + ρ_{2t+1} ≥ 0 (monotone pairs truncation).
#
# Target: ESS > 400, or ESS/N_total > 0.1 per parameter.
#
compute_ess <- function(chains) {
  M <- length(chains)
  N <- nrow(chains[[1L]])
  J <- ncol(chains[[1L]])

  ess <- setNames(numeric(J), param_names_log)

  for (j in seq_len(J)) {
    all_samps <- unlist(lapply(chains, function(ch) ch[, j]))
    n_total   <- length(all_samps)
    max_lag   <- min(n_total - 1L, 500L)

    acf_vals  <- acf(all_samps, lag.max = max_lag,
                     plot = FALSE, type = "correlation")$acf[-1L]

    # Geyer monotone pairs:  Γ_t = ρ_{2t} + ρ_{2t+1}
    n_pairs <- floor(length(acf_vals) / 2L)
    if (n_pairs < 1L) { ess[j] <- n_total; next }

    pairs   <- acf_vals[2L * seq_len(n_pairs) - 1L] +
               acf_vals[2L * seq_len(n_pairs)]
    cutoff  <- which(pairs <= 0)[1L]
    if (is.na(cutoff)) cutoff <- n_pairs

    # Sum autocorrelations up to cut-off pair
    rho_sum <- 1 + 2 * sum(acf_vals[seq_len(2L * (cutoff - 1L) + 1L)])
    ess[j]  <- n_total / max(rho_sum, 1)
  }
  ess
}

# ── 6c. Divergence check ─────────────────────────────────────────────────────
#
# A "divergence" occurs when the leapfrog trajectory explodes (H changes by
# > threshold).  Here we flag iterations where |ΔH| > 1000 as divergent.
# In practice, this is visible as sudden jumps in the lp trace.
#
count_divergences <- function(chains_raw, threshold = 1000) {
  lapply(chains_raw, function(ch) {
    lp  <- ch$lp_trace[ch$n_warmup + seq_len(nrow(ch$samples))]
    sum(abs(diff(lp)) > threshold, na.rm = TRUE)
  })
}

# ── 6d. Print formatted diagnostics table ───────────────────────────────────
print_diagnostics <- function(post_warmup_list, chains_raw = NULL) {
  rhat_vals  <- compute_rhat(post_warmup_list)
  ess_vals   <- compute_ess(post_warmup_list)
  all_samps  <- do.call(rbind, post_warmup_list)
  orig_samps <- exp(all_samps)

  diag_df <- data.frame(
    Parameter  = param_names_log,
    Mean_log   = round(colMeans(all_samps),               4),
    SD_log     = round(apply(all_samps,  2, sd),          4),
    Mean_orig  = round(colMeans(orig_samps),              4),
    Q2.5_orig  = round(apply(orig_samps, 2, quantile, 0.025), 4),
    Q97.5_orig = round(apply(orig_samps, 2, quantile, 0.975), 4),
    Rhat       = round(rhat_vals, 4),
    ESS        = round(ess_vals,  1),
    Conv_OK    = ifelse(rhat_vals < 1.05 & ess_vals > 100, "YES", "NO"),
    stringsAsFactors = FALSE
  )

  cat("\n======================================================\n")
  cat("          Convergence Diagnostics Summary\n")
  cat("  Target: R-hat < 1.05   |   ESS > 100 per parameter\n")
  cat("======================================================\n")
  print(diag_df, row.names = FALSE)

  if (!is.null(chains_raw)) {
    div_counts <- count_divergences(chains_raw)
    cat("\nDivergences per chain (post-warmup):\n")
    print(setNames(unlist(div_counts), paste0("Chain ", seq_along(div_counts))))
  }

  invisible(diag_df)
}


# =============================================================================
# 7.  POSTERIOR PREDICTIVE CHECKS (PPC)
# =============================================================================

# ── 7a. Generate PPC replicates ──────────────────────────────────────────────
#
# For each sampled theta_s, run the ODE model and draw:
#   y_rep_pos[i] ~ Poisson(lambda_pos[i] | theta_s)
#   y_rep_tot[i] ~ Poisson(lambda_tot[i] | theta_s)
#
# Returns replicated counts (integer) and the Poisson means (real),
# plus posterior predictive p-values (ppp) per age group.
#
generate_ppc_samples <- function(post_samples, base_params, data,
                                  n_ppc = 200L) {
  n_avail  <- nrow(post_samples)
  draw_idx <- sort(sample(n_avail, min(n_ppc, n_avail)))
  n_draws  <- length(draw_idx)

  ppc_pos  <- matrix(NA_integer_, n_draws, 9L)
  ppc_tot  <- matrix(NA_integer_, n_draws, 9L)
  lam_pos  <- matrix(NA_real_,    n_draws, 9L)
  lam_tot  <- matrix(NA_real_,    n_draws, 9L)

  cat(sprintf("Generating %d PPC draws from %d posterior samples...\n",
              n_draws, n_avail))

  for (ii in seq_len(n_draws)) {
    theta <- post_samples[draw_idx[ii], ]

    result <- tryCatch({
      pm      <- build_params_from_theta(theta, base_params)
      out     <- run_sim(pm, data)
      y_fin   <- as.numeric(out[nrow(out), -1L])
      compute_age_quantities(y_fin)
    }, error = function(e) NULL)

    if (!is.null(result)) {
      lam_tot[ii, ] <- N_total_obs * result$p_age
      lam_pos[ii, ] <- lam_tot[ii, ] * result$q_age

      ppc_tot[ii, ] <- as.integer(rmultinom(1L, size = N_total_obs, prob = result$p_age))
      p_draw <- rbeta(9L, shape1 = result$q_age * phi_overdisp, shape2 = (1 - result$q_age) * phi_overdisp)
      ppc_pos[ii, ] <- rbinom(9L, size = ppc_tot[ii, ], prob = p_draw)
    }

    if (ii %% 50L == 0L)
      cat(sprintf("  PPC draw %d / %d\n", ii, n_draws))
  }

  # Posterior predictive p-values: Pr(y_rep >= y_obs | y_obs)
  # Calibrated model: ppp near 0.5; poor fit: near 0 or 1.
  ppp_pos <- vapply(1:9, function(i)
    mean(ppc_pos[, i] >= obs_pos[i], na.rm = TRUE), numeric(1L))
  ppp_tot <- vapply(1:9, function(i)
    mean(ppc_tot[, i] >= obs_tot[i], na.rm = TRUE), numeric(1L))

  list(ppc_pos = ppc_pos, ppc_tot = ppc_tot,
       lam_pos = lam_pos, lam_tot = lam_tot,
       ppp_pos = ppp_pos, ppp_tot = ppp_tot)
}

# ── 7b. PPC histogram plot ───────────────────────────────────────────────────
#
# One panel per age group: histogram of y_rep (9 × n_ppc), vertical line at
# observed count, and posterior predictive p-value annotated per panel.
#
plot_ppc_histograms <- function(ppc, type = c("pos", "tot")) {
  type <- match.arg(type)

  ppc_mat  <- if (type == "pos") ppc$ppc_pos else ppc$ppc_tot
  obs_vals <- if (type == "pos") obs_pos     else obs_tot
  ppp_vals <- if (type == "pos") ppc$ppp_pos else ppc$ppp_tot

  fill_col <- if (type == "pos") "#2166AC" else "#D6604D"
  main_ttl <- if (type == "pos")
    "PPC: HCV-Positive Counts (J-stratum)" else
    "PPC: Total PWID Counts (J-stratum)"

  age_lbl <- paste0("Age group ", 1:9)

  # Drop failed draws
  ppc_mat  <- ppc_mat[complete.cases(ppc_mat), , drop = FALSE]

  ppc_long <- as.data.frame(ppc_mat) %>%
    setNames(age_lbl) %>%
    pivot_longer(everything(), names_to = "Age", values_to = "y_rep") %>%
    mutate(Age = factor(Age, levels = age_lbl))

  obs_df <- data.frame(
    Age = factor(age_lbl, levels = age_lbl),
    obs = obs_vals
  )
  ppp_df <- data.frame(
    Age = factor(age_lbl, levels = age_lbl),
    ppp = round(ppp_vals, 2)
  )

  ggplot(ppc_long, aes(x = y_rep)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins    = 30,
                   fill    = fill_col,
                   alpha   = 0.65,
                   colour  = "white",
                   linewidth = 0.2) +
    geom_vline(data = obs_df, aes(xintercept = obs),
               colour = "black", linewidth = 1.1) +
    geom_text(data = ppp_df,
              aes(label = paste0("ppp=", ppp)),
              x = Inf, y = Inf, hjust = 1.15, vjust = 1.6,
              size = 3.0, colour = "grey25") +
    facet_wrap(~ Age, scales = "free", ncol = 3) +
    labs(
      title    = main_ttl,
      subtitle = "Histogram: posterior predictive  |  Vertical line: observed  |  ppp near 0.5 = good fit",
      x        = "Count",
      y        = "Density"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background  = element_rect(fill = "grey92", colour = "grey70"),
      panel.grid.minor  = element_blank(),
      plot.subtitle     = element_text(size = 9, colour = "grey40")
    )
}

# ── 7c. PPC interval plot ─────────────────────────────────────────────────────
#
# Dot-and-whisker: posterior predictive median ± 50% and 90% intervals per age
# group, with observed count overlaid as a red point.
#
plot_ppc_intervals <- function(ppc) {
  age_lbl <- paste0("Age ", 1:9)

  summarise_ppc_mat <- function(ppc_mat, obs, type_label) {
    mat <- ppc_mat[complete.cases(ppc_mat), , drop = FALSE]
    data.frame(
      Age    = factor(age_lbl, levels = age_lbl),
      obs    = obs,
      med    = apply(mat, 2, median),
      lo50   = apply(mat, 2, quantile, 0.25),
      hi50   = apply(mat, 2, quantile, 0.75),
      lo90   = apply(mat, 2, quantile, 0.05),
      hi90   = apply(mat, 2, quantile, 0.95),
      type   = type_label,
      stringsAsFactors = FALSE
    )
  }

  plot_df <- bind_rows(
    summarise_ppc_mat(ppc$ppc_pos, obs_pos, "HCV-positive"),
    summarise_ppc_mat(ppc$ppc_tot, obs_tot, "Total PWID")
  )

  ggplot(plot_df, aes(x = Age)) +
    geom_linerange(aes(ymin = lo90, ymax = hi90),
                   linewidth = 1.0, colour = "#4292C6", alpha = 0.45) +
    geom_linerange(aes(ymin = lo50, ymax = hi50),
                   linewidth = 3.0, colour = "#08519C", alpha = 0.75) +
    geom_point(aes(y = med),
               shape = 18, size = 3.5, colour = "#08306B") +
    geom_point(aes(y = obs),
               shape = 21, size = 3.5,
               fill = "#D7191C", colour = "black", stroke = 1.2) +
    facet_wrap(~ type, scales = "free_y", ncol = 1) +
    labs(
      title    = "PPC: Posterior Predictive Intervals vs Observed",
      subtitle = "Diamond: predictive median  |  Thick bar: 50% PI  |  Thin bar: 90% PI  |  Red dot: observed",
      x        = "Age group",
      y        = "Count"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x       = element_text(angle = 35, hjust = 1),
      strip.background  = element_rect(fill = "grey92", colour = "grey70"),
      panel.grid.minor  = element_blank(),
      plot.subtitle     = element_text(size = 9, colour = "grey40")
    )
}

# ── 7d. Trace plots ───────────────────────────────────────────────────────────
#
# Full-chain trace for selected parameters (all iterations, all chains),
# with a dashed vertical line marking the end of warmup.
#
plot_traces <- function(chains_raw, param_idx = 1:min(5L, N_PARAMS)) {

  trace_list <- lapply(seq_along(chains_raw), function(ch) {
    df           <- as.data.frame(chains_raw[[ch]]$samples_all)
    colnames(df) <- param_names_log
    df$iter      <- seq_len(nrow(df))
    df$chain     <- factor(ch)
    df
  })
  trace_df  <- bind_rows(trace_list)
  n_warmup  <- chains_raw[[1L]]$n_warmup
  params_to <- param_names_log[param_idx]

  trace_df %>%
    select(iter, chain, all_of(params_to)) %>%
    pivot_longer(all_of(params_to),
                 names_to  = "param",
                 values_to = "value") %>%
    mutate(param = factor(param, levels = params_to)) %>%
    ggplot(aes(x = iter, y = value, colour = chain)) +
    geom_line(alpha = 0.70, linewidth = 0.35) +
    geom_vline(xintercept = n_warmup,
               linetype = "dashed", colour = "grey35", linewidth = 0.8) +
    facet_wrap(~ param, scales = "free_y", ncol = 1) +
    scale_colour_brewer(palette = "Set1") +
    labs(
      title    = "Trace plots — all chains (log scale)",
      subtitle = "Dashed line = end of warmup",
      x        = "Iteration",
      y        = "Value (log-transformed parameter)",
      colour   = "Chain"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background  = element_rect(fill = "grey92", colour = "grey70"),
      panel.grid.minor  = element_blank(),
      legend.position   = "right",
      plot.subtitle     = element_text(size = 9, colour = "grey40")
    )
}

# ── 7e. Posterior density plot (original scale) ────────────────────────────
plot_posterior_densities <- function(post_warmup_list) {
  all_samps <- do.call(rbind, lapply(seq_along(post_warmup_list), function(ch) {
    df           <- as.data.frame(exp(post_warmup_list[[ch]]))
    colnames(df) <- param_names_orig
    df$chain     <- factor(ch)
    df
  }))

  all_samps %>%
    pivot_longer(-chain, names_to = "param", values_to = "value") %>%
    mutate(param = factor(param, levels = param_names_orig)) %>%
    ggplot(aes(x = value, fill = chain, colour = chain)) +
    geom_density(alpha = 0.30, linewidth = 0.6) +
    facet_wrap(~ param, scales = "free", ncol = 2) +
    scale_fill_brewer(palette  = "Set1") +
    scale_colour_brewer(palette = "Set1") +
    labs(
      title    = "Posterior densities (original scale, all chains)",
      x        = "Parameter value",
      y        = "Density",
      fill     = "Chain",
      colour   = "Chain"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "grey92", colour = "grey70"),
      panel.grid.minor = element_blank(),
      legend.position  = "right"
    )
}


# =============================================================================
# 8.  MAIN: RUN HMC
# =============================================================================

# ── Sampler settings ──────────────────────────────────────────────────────────
N_CHAINS    <- 4L      # parallel chains for R-hat / ESS
N_CORES     <- 10L      # cores for parallel gradient batches (set to 1 for debugging)
N_WARMUP    <- 200L    # adaptation (discarded)
N_ITER      <- 1000L    # total iterations per chain  (post-warmup = N_ITER - N_WARMUP)
EPS_INIT    <- 0.01    # initial step size (dual averaging will adapt)
L_STEPS     <- 10L     # leapfrog steps per proposal
ADAPT_DELTA <- 0.65    # target acceptance rate

# ── Initial points ────────────────────────────────────────────────────────────
set.seed(114514)
inits <- lapply(seq_len(N_CHAINS), function(ch) {
  c(
    log(runif(1L, 0.01, 1.0)),      # log(beta_scale):    mean=0 prior, start in (0,1)
    log(runif(1L, 0.05, 0.5)),      # log(delta):         shift/floor, start small and positive
    log(runif(1L, 1.5, 4.0)),       # log(alpha):         Gamma shape hyperparameter
    log(runif(1L, 0.5, 2.0)),       # log(beta_rate):     Gamma rate hyperparameter
    log(runif(8L, 0.5, 3.0))        # log(C_contact_scale[1:8]): free row scalings
  )
})

# ── Run chains (sequential; replace with parallel::mclapply for speed) ────────
cat("=== HMC Calibration: HCV PWID Model ===\n")
cat(sprintf("Chains: %d | Warmup: %d | Sampling: %d | L: %d\n",
            N_CHAINS, N_WARMUP, N_ITER - N_WARMUP, L_STEPS))
cat(sprintf("Approx ODE runs per HMC step: %d\n", (L_STEPS + 1L) * 2L * N_PARAMS))

hmc_chains <- parallel::mclapply(seq_len(N_CHAINS), function(ch) {
  run_hmc_chain(
    theta_init  = inits[[ch]],
    n_iter      = N_ITER,
    n_warmup    = N_WARMUP,
    eps_init    = EPS_INIT,
    L           = L_STEPS,
    adapt_delta = ADAPT_DELTA,
    base_params = params,
    data        = data,
    seed        = ch * 314L,
    chain_id    = ch
  )
}, mc.cores = N_CORES)

# ── Pool post-warmup samples ──────────────────────────────────────────────────
post_warmup_list <- lapply(hmc_chains, function(ch) ch$samples)
all_samples      <- do.call(rbind, post_warmup_list)

cat(sprintf("\nTotal post-warmup samples: %d (%d × %d chains)\n",
            nrow(all_samples), nrow(post_warmup_list[[1L]]), N_CHAINS))

# ── Diagnostics ───────────────────────────────────────────────────────────────
diag_table <- print_diagnostics(post_warmup_list, chains_raw = hmc_chains)

# ── Acceptance rates ──────────────────────────────────────────────────────────
cat("\nAcceptance rates per chain (post-warmup):\n")
acc_df <- data.frame(
  Chain       = seq_len(N_CHAINS),
  Acc_rate    = round(vapply(hmc_chains, function(ch) ch$acc_rate, numeric(1L)), 3),
  Eps_final   = round(vapply(hmc_chains, function(ch) ch$eps_final, numeric(1L)), 5)
)
print(acc_df, row.names = FALSE)

# ── Posterior summary (original scale) ───────────────────────────────────────
orig_samples <- exp(all_samples)
post_summary <- data.frame(
  Parameter  = param_names_orig,
  Mean       = round(colMeans(orig_samples), 4),
  SD         = round(apply(orig_samples, 2, sd), 4),
  Q2.5       = round(apply(orig_samples, 2, quantile, 0.025), 4),
  Q50        = round(apply(orig_samples, 2, quantile, 0.500), 4),
  Q97.5      = round(apply(orig_samples, 2, quantile, 0.975), 4)
)
cat("\n=== Posterior Summary (original scale) ===\n")
print(post_summary, row.names = FALSE)

# ── PPC ───────────────────────────────────────────────────────────────────────
cat("\n=== Generating Posterior Predictive Samples ===\n")
ppc_out <- generate_ppc_samples(all_samples, params, data, n_ppc = 200L)

cat("\nPosterior predictive p-values (near 0.5 = good fit):\n")
ppp_df <- data.frame(
  Age     = paste0("Age ", 1:9),
  ppp_pos = round(ppc_out$ppp_pos, 3),
  ppp_tot = round(ppc_out$ppp_tot, 3)
)
print(ppp_df, row.names = FALSE)

# ── Plots ─────────────────────────────────────────────────────────────────────
p_trace    <- plot_traces(hmc_chains, param_idx = 1:N_PARAMS)
p_density  <- plot_posterior_densities(post_warmup_list)
p_hist_pos <- plot_ppc_histograms(ppc_out, type = "pos")
p_hist_tot <- plot_ppc_histograms(ppc_out, type = "tot")
p_interval <- plot_ppc_intervals(ppc_out)

print(p_trace)
print(p_density)
print(p_hist_pos)
print(p_hist_tot)
print(p_interval)

# ── Save ─────────────────────────────────────────────────────────────────────
saveRDS(
  list(
    hmc_chains       = hmc_chains,
    post_warmup_list = post_warmup_list,
    all_samples      = all_samples,
    diag_table       = diag_table,
    post_summary     = post_summary,
    ppc_out          = ppc_out
  ),
  file = "hmc_output.rds"
)
cat("\nAll results saved to hmc_output.rds\n")