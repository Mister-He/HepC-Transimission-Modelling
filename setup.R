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
  q      = 0.0057,          # transmission prob per needle-sharing event
  kappa  = 0.26,            # prob spontaneous clearance after acute infection
  iota1  = 180 / 365,       # mean acute duration (years); 180 days
  iota2  =  84 / 365,       # mean treatment duration (years); 84 days = 12 weeks

  # ── Genotype ──────────────────────────────────────────────────────────────
  rho    = 0.78,            # proportion with genotype 3

  # ── SVR12 rates (effective, genotype-weighted) ─────────────────────────────
  # Stages 1-2 (NC & CC), RBV-free
  # alpha_NC = rho * alpha3_NC + (1-rho) * alphao_NC
  #          = 0.78*0.976 + 0.22*0.988 = 0.9787
  alpha_NC     = 0.9787,

  # Stage 3 (DC) with RBV+
  # alpha_DC_pos = 0.78*0.87 + 0.22*0.92 = 0.8818
  alpha_DC_pos = 0.8818,

  # Stage 3 (DC) without RBV (used in sensitivity)
  # alpha_DC_neg = 0.78*0.45 + 0.22*0.90 = 0.5490
  alpha_DC_neg = 0.5490,

  # Stage 4 (HCC): assumed equal to DC+ (no dedicated trial data; see Section X)
  alpha_HCC    = 0.8818,

  # ── Treatment initiation rates (per year, per stage) ──────────────────────
  # GUESS: set to 0 (no treatment) as base scenario; override in scenarios
  tau = c(
    0.0,   # tau_1: NC   — SCENARIO-DEFINED; 0 = no treatment
    0.0,   # tau_2: CC   — SCENARIO-DEFINED
    0.0,   # tau_3: DC   — SCENARIO-DEFINED
    0.0    # tau_4: HCC  — SCENARIO-DEFINED
  ),

  # ── Baseline progression rates (per year, other genotype) ─────────────────
  p_NC_CC   = 0.027,   # NC  → CC  (Thein et al. 2008)
  p_CC_DC   = 0.039,   # CC  → DC  (Lim et al. 2018)
  p_CC_HCC  = 0.014,   # CC  → HCC (Lim et al. 2018)
  p_DC_HCC  = 0.014,   # DC  → HCC (Lim et al. 2018)

  # ── GT3 relative risks ────────────────────────────────────────────────────
  r3_NC_CC   = 1.36,   # r3^{NC→CC}  derived; CI: 1.27–1.46
  r3_CC_DC   = 1.36,   # r3^{CC→DC}  derived; CI: 1.27–1.46
  r3_CC_HCC  = 1.93,   # r3^{CC→HCC} derived; CI: 1.71–2.17
  r3_DC_HCC  = 1.93,   # r3^{DC→HCC} derived; CI: 1.71–2.17

  # ── Post-SVR progression modifiers ────────────────────────────────────────
  phi_CC_DC   = 0.07,  # relative risk CC→DC after SVR (van Meer 2012)
  phi_CC_HCC  = 0.23,  # relative risk CC→HCC after SVR (Morgan et al. 2013)
  phi_DC_HCC  = 1.00,  # assumed no reduction DC→HCC after SVR

  # ── Background mortality (per year, age-varying) ──────────────────────────
  # CALIBRATED: placeholder values from SingStat life table
  # Age groups (example boundaries): 15-19, 20-24, 25-29, 30-34,
  #                                   35-39, 40-44, 45-49, 50-54, 55+
  mu = c(
    0.001267,  # age group 1 (15-19)  
    0.000300,  # age group 2 (20-24) 
    0.000300,  # age group 3 (25-29) 
    0.000400,  # age group 4 (30-34) 
    0.000500,  # age group 5 (35-39) 
    0.000700,  # age group 6 (40-44) 
    0.001400,  # age group 7 (45-49) 
    0.002300,  # age group 8 (50-54) 
    0.016100   # age group 9 (55+)   
  ),

  # ── Standardized Mortality rate of ever-PWIDs  ─────────────────────────────
  omega = 14.68,  # SMR for ever-PWIDs (Degenhardt et al. 2011)

  # ── Disease-specific excess mortality ──────────────────────────────────────
  mu_DC   = 0.130,  # additional mortality in decompensated cirrhosis (Lim 2018)
  mu_HCC  = 0.430,  # additional mortality in HCC                     (Lim 2018)

  # ── SVR mortality modifiers ────────────────────────────────────────────────
  psi_DC  = 0.45,   # relative mortality risk DC after SVR  (Premkumar 2024)
  psi_HCC = 0.37,   # relative mortality risk HCC after SVR (Dang 2020)

  # ── Incarceration rates (per year, age-varying) ────────────────────────────
  # CALIBRATED: placeholder values — replace with SPS-fitted rates
  lambda1 = c(0.9769631, 0.7962602, 0.7842585, 
              0.8721055, 0.8675996, 0.9522804, 
              1.0377386, 0.9348018, 0.9021672),  # first-arrest rate   lambda_i^(1) — GUESS
  lambda2 = c(0.489, 0.620, 0.663,
              0.628, 0.533, 0.475, 
              0.472, 0.441, 0.451),  # release rate        lambda_i^(2) — GUESS (0.5yr avg)
  lambda3 = c(0.9769631, 0.7962602, 0.7842585, 
              0.8721055, 0.8675996, 0.9522804, 
              1.0377386, 0.9348018, 0.9021672),  # re-arrest rate      lambda_i^(3) — GUESS
  pi_recid = 0.61610608,          # recidivism probability (fitted to SPS; Assumption)

  # ── Needle-sharing contact rate ────────────────────────────────────────────
  # CALIBRATED: scalar homogeneous mixing — replace with 9×9 matrix post-calib.
  C_contact = rbind(c(7, 4, 1, 1, 0, 1, 1, 0, 0),
                    c(11, 34, 21, 11, 6, 2, 1, 0.5, 0.5),
                    c(7, 30, 80, 62, 30, 10, 2, 1.5, 1.5), 
                    c(2, 10, 60, 121, 65, 38, 15, 2, 2), 
                    c(1, 11, 22, 67, 107, 41, 18, 2.5, 2.5),
                    c(0, 4, 6, 22, 32, 31, 10, 2, 2),
                    c(0, 1, 1, 8, 10, 15, 11, 0.5, 0.5),
                    c(0, 5, 5, 5.5, 6.5, 6.5, 3.5, 1.5, 1.5),
                    c(0, 5, 5, 5.5, 6.5, 6.5, 3.5, 1.5, 1.5)) * 4,

  # ── Population entry rates (per year, age-varying) ────────────────────────
  # CALIBRATED: constant-in-time placeholder — replace with beta_i(t) from calib.
  beta = c(
    235,        # age group 1
    565/2,      # age group 2 
    565/2,      # age group 3 
    301/2,      # age group 4 
    301/2,      # age group 5 
    111/2,      # age group 6 
    111/2,      # age group 7 
    33/2,       # age group 8 
    33/2 + 4    # age group 9 
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
pos = c(55, 145, 183, 164, 212, 299, 222, 190, 133)
tot = c(307, 797, 829, 633, 598, 642, 481, 439, 366)
for (i in 0:8) {
  y0[idx(s=0, k=1, h=0, i=i)] <- tot[i+1] - pos[i+1]  # D_{u,1,i} 
  y0[idx(s=0, k=1, h=1, i=i)] <- pos[i+1]             # D_{a,1,i}
}

# =============================================================================
# DATA / SIMULATION SETTINGS
# =============================================================================
data <- list(
  t_start = 0.0,    # start year (0 = model year 0; map to calendar year in R)
  t_end   = 100.0,   # simulate 100 years
  dt      = 1/365,   # daily time steps (1/365 of a year)
  y0      = y0      # initial conditions (length-576 vector)
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
stage_names  <- c("NC", "CC", "DC", "HCC")
state_names  <- c("u", "a", "c", "t")
age_names    <- paste0("age", 1:9)

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
print(paste("Simulation time:", round(difftime(end, start, units="secs"), 2), "seconds"))
colnames(out) <- col_names
print(paste("Total population:", round(tail(rowSums(out[, -1]),1), 2)))


plot(out[, "time"], rowSums(out[, grep("J_u_", colnames(out))]),
  type = "l", xlab = "Year", ylab = "Total susceptible",
  main = "Status quo — no treatment")

# plot(out[,"time"], rowSums(out[, grep("_a_age9", colnames(out))]),
#      type = "l", xlab = "Year", ylab = "Total acute HCV",
#      main = "Status quo — no treatment")

# plot(out[,"time"], rowSums(out[, grep("_c_", colnames(out))]),
#      type = "l", xlab = "Year", ylab = "Total chronic HCV",
#      main = "Status quo — no treatment")
# plot(out[,"time"], rowSums(out[, grep("_t_", colnames(out))]),
#      type = "l", xlab = "Year", ylab = "Total treated",
#      main = "Status quo — no treatment")
    
# =============================================================================
# MODEL CALIBRATION
# =============================================================================
# TODO: implement calibration procedure to fit beta_i and C_contact to SPS data
# Firstly, we need to ensure the model has reached equalibrum
# scaling factor of inflow beta
# initial population scaling factor


# Secondly, simulation results should be compared to observed data in 2017 with CIs