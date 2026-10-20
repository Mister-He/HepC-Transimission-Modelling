# =============================================================================
# tests/unit/test_likelihood.R — unit tests for the 12-parameter likelihood and parameterisation
# (likelihood.R). Run: Rscript tests/unit/test_likelihood.R
# =============================================================================

.args0 <- commandArgs(FALSE)
.file_arg <- sub("^--file=", "", .args0[grepl("^--file=", .args0)])
root <- normalizePath(file.path(dirname(normalizePath(.file_arg)), "..", ".."))
source(file.path(root, "tests", "helper.R"))
source(file.path(root, "scripts", "calibration", "targets.R"))
targets <- read_targets(file.path(root, "data", "calibration_targets.json"))
source(file.path(root, "scripts", "calibration", "likelihood.R"))

# Synthetic base parameters: the unit test exercises the parameterisation
# (exp transforms, bounds, likelihoods), not the fitted values.
base_params <- list(
  C_contact = matrix(1, nrow = 6, ncol = 6),
  beta = rep(1, 6), beta_base = rep(1, 6)
)

check_equal(N_THETA, 12, "exactly 12 fitted parameters")
check_true(all(CONTACT_LO > 0 & CONTACT_LO < CONTACT_HI),
           "contact bounds are valid")
check_true(all(BETA_LO > 0 & BETA_LO < BETA_HI),
           "beta bounds are valid")

# theta = 0  =>  all scales = 1  =>  params unchanged
pm0 <- build_params(rep(0, N_THETA), base_params)
check_equal(pm0$C_contact, base_params$C_contact,
            "theta = 0 leaves C_contact unchanged")
check_equal(pm0$beta, base_params$beta,
            "theta = 0 leaves beta unchanged")

# exp transform: a single unit shift scales one row / one inflow by e
th1 <- rep(0, N_THETA); th1[1] <- 1
pm1 <- build_params(th1, base_params)
check_equal(pm1$C_contact[1, ], exp(1) * base_params$C_contact[1, ],
            "theta[1] = 1 scales contact row 1 by e")
check_equal(pm1$C_contact[2, ], base_params$C_contact[2, ],
            "other contact rows unchanged")
th7 <- rep(0, N_THETA); th7[7] <- 1
pm7 <- build_params(th7, base_params)
check_equal(pm7$beta[1], exp(1) * base_params$beta[1],
            "theta[7] = 1 scales beta group 1 by e")

# bounds
check_true(within_bounds(rep(0, N_THETA)), "zero theta is within bounds")
check_true(!within_bounds(rep(50, N_THETA)), "large theta is out of bounds")

# Binomial likelihood: perfect prediction beats a poor one
check_true(nll_prev(targets$p_prev_binom) <
             nll_prev(rep(0.01, 6)),
           "perfect prevalence prediction has lower NLL")
check_true(is.finite(nll_prev(targets$p_prev_binom)),
           "NLL is finite at the target prevalence")

# Population likelihood: exact match gives zero contribution
check_equal(nll_pop(targets$population), 0,
            "nll_pop is 0 when N_hat == N_obs")

# Plausibility penalty: a diagonal peaking at 30-39 with small 50+/60+
# entries and non-decreasing beta scales scores 0.
diag_ok <- base_params$C_contact
diag_ok[1, 1] <- 0.3; diag_ok[2, 2] <- 0.8; diag_ok[3, 3] <- 1.0
diag_ok[4, 4] <- 0.6; diag_ok[5, 5] <- 0.3; diag_ok[6, 6] <- 0.2
pm_ok <- base_params; pm_ok$C_contact <- diag_ok; pm_ok$beta <- seq(1, 3, length.out = 6)
th_ok <- log(c(rep(1, 6), seq(1, 3, length.out = 6)))
check_true(pattern_penalty(pm_ok, th_ok) == 0,
           "pattern penalty is 0 for a realistic contact/beta pattern")

# Mis-specified shape must attract a penalty
diag_bad <- diag_ok
diag_bad[3, 3] <- 0.2   # 30-39 below the 40-49 value
pm_bad <- base_params; pm_bad$C_contact <- diag_bad; pm_bad$beta <- seq(3, 1, length.out = 6)
th_bad <- log(c(rep(1, 6), seq(3, 1, length.out = 6)))  # decreasing beta
check_true(pattern_penalty(pm_bad, th_bad) > 0,
           "pattern penalty is positive for an implausible pattern")

# Entry scale: theta[7:12] multiplies every entry group by exp(theta) and nothing else.
fp <- modifyList(base_params, list(beta = c(10, 300, 100, 50, 20, 1)))
th_rt <- c(rep(0, 6), log(c(0.5, 2, 1, 3, 1, 4)))
check_equal(build_params(th_rt, fp)$beta, fp$beta * c(0.5, 2, 1, 3, 1, 4),
            "build_params sets beta = beta_base * exp(theta[7:12])")
check_true(!within_bounds(c(rep(0, 6), rep(log(0.005), 6))), "entry scales below 0.01 are out of bounds")

# build_params rejects the wrong number of parameters
check_error(build_params(rep(0, 11), base_params),
            "build_params rejects length-11 theta")

# Prior centre: interior fits are unchanged; a fit on a bound moves 5 prior SDs inside it.
th_in <- log(c(rep(1, 6), rep(2, 6)))
check_equal(prior_centre(th_in), th_in, "prior centre equals an interior best fit")
th_pin <- th_in; th_pin[6] <- log(CONTACT_LO[6])
pc <- prior_centre(th_pin)
check_equal(pc[6], log(CONTACT_LO[6]) + 1, "a bound-pinned fit is moved 5 prior SDs inside the bound")
check_equal(pc[-6], th_pin[-6], "only the pinned coordinate moves")
check_true(pnorm(log(CONTACT_LO[6]), pc[6], 0.20) < 1e-6, "almost no prior mass falls outside the bound")

# Nelder-Mead restarts never end worse and count every evaluation in the history.
source(file.path(root, "scripts", "calibration", "calibrate_nm.R"))
rosen <- function(x) sum(100 * (x[-1] - x[-12]^2)^2 + (1 - x[-12])^2)
nm0 <- run_nm_start(rep(-1, 12), rosen, maxit = 500)
nm3 <- run_nm_start(rep(-1, 12), rosen, maxit = 500, restarts = 3)
check_true(nm3$objective <= nm0$objective && nm3$fn_evaluations == nrow(nm3$history),
           "Nelder-Mead restarts do not worsen the fit and log every evaluation")

n_fail <- test_report("tests/unit/test_likelihood.R")
if (n_fail > 0) quit(status = 1)
