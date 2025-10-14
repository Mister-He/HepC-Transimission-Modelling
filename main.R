#######################################################
# Oct 14, 2025
# 1. Estimate age-specific scaling factors for Contact 
#   Matrix using MCMC 
# 2. Convergence Diagnostics using G-R and ESS
######################################################
library(data.table)
library(tictoc)
library(dplyr)
library(ggplot2)

# =========================================#
# 1. Preparation
# =========================================#
cat("======Loading libraries and data...======\n")
# Load required functions and data
source("functions/utils.r")
strategies <- read.csv("settings/Final_strategies.csv")

# Change annual proportions to monthly proportions
strategies$pD <- strategies$pD / 12
strategies$pJ1 <- strategies$pJ1 / 12
strategies$pJ2 <- strategies$pJ2 / 12
strategies$pJ3 <- strategies$pJ3 / 12
strategies$pX <- strategies$pX / 12

# Load observed HCV prevalence by age group
obs_data <- list(
    pos = c(55, 145, 183, 164, 212, 299, 222, 190, 133),
    tot = c(307, 797, 829, 633, 598, 642, 481, 439, 366)
)

# Set up to calculate logLL, which can be used for MCMC
nsteps_permonth <- 1 # 4 # 30
intervention <- strategies[1, ]
configuration <- create_configuration(ny = 600, nsteps_permonth = nsteps_permonth) # 250
parameters <- set_parameters(nsteps_permonth = configuration$nsteps_permonth)
initialstates <- make_initial_states(configuration, equilibrium = FALSE, parameters)

# Get this from some personal guess
xparam_init_scaling <- log(c(9.0, 0.2, 0.18, 0.13, 0.44, 0.3, 0.08, 0.02, 0.01))
xparam_init_arrest <- log(c(0.067, 0.21, 0.21, 0.15, 0.2, 0.31, 0.14, 0.4, 0.025) / nsteps_permonth / 12)
xparam_init <- c(xparam_init_scaling, xparam_init_arrest)

# Set bounds for parameters
LB_scaling <- log(c(8.5, 0.1, 0.1, 0.1, 0.35, 0.1, 0.06, 0.018, 0.008))
UB_scaling <- log(c(9.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.1, 0.022, 0.012))

LB_arrest <- xparam_init_arrest - abs(xparam_init_arrest) * 0.1
UB_arrest <- xparam_init_arrest + abs(xparam_init_arrest) * 0.1

LB <- c(LB_scaling, LB_arrest)
UB <- c(UB_scaling, UB_arrest)

labs <- c(paste("<", 20, sep = ""), paste(seq(20, 50, by = 5), seq(20 + 5 - 1, 55 - 1, by = 5), sep = "-"), paste(55, "+", sep = ""))

# =========================================#
# 2. Markov Chain Monte Carlo (MCMC)
# =========================================#
library(mcmc)
library(parallel)
library(coda)
library(future)
library(future.apply)

# Parameters
set.seed(114514)
n_chains  <- 4
n_burnin  <- 3e+4
n_samples <- 2e+4
n_batch_length <- 100
n_batches <- (n_burnin + n_samples) / n_batch_length

# Create some sensible inits for each parallel chain
inits <- replicate(n_chains, runif(length(LB), min = LB, max = UB), simplify = FALSE)

# 2.1 Use grid search to find a good scale for pilot
target_accept <- 0.25
scales_to_try <- c(0.01, 0.02, 0.03, 0.04, 0.05)
accept_rates  <- numeric(length(scales_to_try))

cat("======Start evaluating scalar scales...======\n")
for (i in seq_along(scales_to_try)) {
    test_pilot <- metrop(fun_logLL_MCMC, initial = xparam_init, nbatch = 200, blen = 1, scale = scales_to_try[i])
    accept_rates[i] <- test_pilot$accept
    cat("Scale =", scales_to_try[i], ", Accept rate =", accept_rates[i], "\n")
}

# Choose the best scale
best_scale_idx <- which.min(abs(accept_rates - target_accept))
best_scale <- scales_to_try[best_scale_idx]
cat("Best scale:", best_scale, "with accept rate:", accept_rates[best_scale_idx], "\n")

# Run pilot with the best scale
cat("======Run pilot with the best scale...======\n")
pilot <- metrop(fun_logLL_MCMC, initial = xparam_init, nbatch = 1e3, blen = 1, scale = best_scale)
cat("Final pilot acceptance rate:", pilot$accept, "\n")

# Double check if the accept rate is still sensible
if (pilot$accept > 0.5) {
    cat("Warning: Acceptance rate is still too high\n")
    # Lower scale
    pilot <- metrop(fun_logLL_MCMC, initial = xparam_init, nbatch = 1e3, blen = 1, scale = best_scale * 0.1)
    cat("Reduced scale pilot acceptance rate:", pilot$accept, "\n")
} else if (pilot$accept < 0.2) {
    cat("Warning: Acceptance rate is too low\n")
    # Increase scale
    pilot <- metrop(fun_logLL_MCMC, initial = xparam_init, nbatch = 1e3, blen = 1, scale = best_scale * 2)
    cat("Increased scale pilot acceptance rate:", pilot$accept, "\n")
}

# 2.2 Get covariance matrix to further scaling the previous best scale
cat("======Start calculating covariance matrix...======\n")
# Calculate covariance matrix, using a larger regularization
S <- cov(pilot$batch[501:1e3, ]) + 1e-6 * diag(length(xparam_init))
L <- chol(S)   
d <- length(xparam_init)
c.opt <- 2.38 / sqrt(d)
scale_mat <- c.opt * t(L)

# Test final scale_vec
final_test <- metrop(fun_logLL_MCMC, initial = xparam_init, nbatch = 1000, blen = 1, scale = scale_mat)
cat("Final test acceptance rate:", final_test$accept, "\n")

cat("======Start parallel MCMC...======\n")
time1 <- Sys.time()
run_chain <- function(seed, init) {
  set.seed(seed)
  fit <- metrop(fun_logLL_MCMC, 
                initial = init, 
                nbatch = n_batches,
                blen = n_batch_length,
                scale = scale_mat)
  return(fit)
}

plan(multisession, workers = n_chains)
fits <- future_lapply(1:n_chains, function(k) run_chain(k, inits[[k]]), future.seed = TRUE)
saveRDS(fits, file = "output/mcmc_raw_results.rds")

# Print acceptance rates for each chain
for (k in 1:n_chains) {
    cat("Chain", k, "acceptance rate:", fits[[k]]$accept, "\n")
}

mcmc_list <- mcmc.list(lapply(fits, function(f) mcmc(f$batch[-c(1:n_burnin/n_batch_length),])))
cat("======Parallel MCMC completes!...======\n")

# =========================================#
# 3. Convergence Diagnostics
# =========================================#
rhat <- gelman.diag(mcmc_list)  # Gelman-Rubin Diagnostics
ess <- effectiveSize(mcmc_list) # Effective Sample Size

cat("======Gelman-Rubin Diagnostics======\n")
print(rhat$psrf)

cat("======Effective Sample Size======\n")
print(ess)

cat("======Bundle Results======\n")
results <- list(
    samples = mcmc_list,
    Rhat    = rhat,
    ESS     = ess
)

# Save to RDS
saveRDS(results, file = "output/mcmc_results.rds")

time2 <- Sys.time()
cat("======Done======\n")
cat("It takes", difftime(time2, time1, units = "mins"), "minutes in total.\n")

# =========================================#
# 4. Post-result analysis
# =========================================#
results <- readRDS('output/mcmc_results.rds')
scaling_params <- do.call(rbind, results$samples)

# 4.1 Save the equilibrium status
nsteps_permonth <- 1
intervention <- strategies[1, ]
configuration <- create_configuration(ny = 600, nsteps_permonth = nsteps_permonth)
parameters <- set_parameters(nsteps_permonth = configuration$nsteps_permonth)
initialstates <- make_initial_states(configuration, equilibrium = FALSE, parameters)
X <- epimodel_fast(initialstates, parameters,
    configuration, intervention,
    adjContactMat = colMeans(scaling_params)
)
final_prevalence(X)
store_states(X, outfile = "output/equilibrium_241410_low_SVR.rdata")

# 4.2 Assess intervention impact using posterior mean
intervention_file_1 <- "settings/tableHIGH_1_wFINALstrategies_yichen_low_SVR.csv"
tabulator(file = intervention_file_1, header = TRUE)
for (str in 1:dim(strategies)[1]) {
    cat("Running strategy", str, "of", dim(strategies)[1], "\n")
    intervention <- strategies[str, ]
    configuration <- create_configuration(50)
    parameters <- set_parameters(nsteps_permonth = configuration$nsteps_permonth)
    parameters$theta_1 <- parameters$theta_1 * 0.8
    parameters$theta_2 <- parameters$theta_2 * 0.8
    initialstates <- make_initial_states(configuration, equilibrium = TRUE, parameters)
    X <- epimodel_fast(initialstates, parameters,
        configuration, intervention,
        adjContactMat = colMeans(scaling_params)
    )
    tabulator(X, intervention, file = intervention_file_1, resolution = "HIGH")
    print(final_prevalence(X))
}

# 4.3 Interventions with Bayesian Credible Intervals
for (iteration in 1:100) {
    # Sample once from the posterior distribution
    pos_sample <- apply(scaling_params, 2, function(x) sample(x, 1))

    # Create a list to save all results in this iteration
    result <- list()

    for (str in 1:2) {
        cat("Running strategy", str, "of", dim(strategies)[1], "in iteration", iteration, "\n")
        intervention <- strategies[str, ]
        configuration <- create_configuration(50)
        parameters <- set_parameters(nsteps_permonth = configuration$nsteps_permonth)
        initialstates <- make_initial_states(configuration, equilibrium = TRUE, parameters)
        X <- epimodel_fast(initialstates, parameters,
            configuration, intervention,
            adjContactMat = pos_sample
        )
        print(final_prevalence(X))
        result[[str]] <- rdstabulator(X, intervention)
    }
    sim_data <- do.call(cbind, result)
    colnames(simdata) <- c("Strategy,t,Du1,Du2,Du3,Du4,Du5,Du6,Du7,Du8,Du9,D01,D02,D03,D04,D05,D06,D07,D08,D09,D11,D12,D13,D14,D15,D16,D17,D18,D19,D21,D22,D23,D24,D25,D26,D27,D28,D29,D31,D32,D33,D34,D35,D36,D37,D38,D39,D41,D42,D43,D44,D45,D46,D47,D48,D49,D51,D52,D53,D54,D55,D56,D57,D58,D59,D61,D62,D63,D64,D65,D66,D67,D68,D69,Ju1,Ju2,Ju3,Ju4,Ju5,Ju6,Ju7,Ju8,Ju9,J01,J02,J03,J04,J05,J06,J07,J08,J09,J11,J12,J13,J14,J15,J16,J17,J18,J19,J21,J22,J23,J24,J25,J26,J27,J28,J29,J31,J32,J33,J34,J35,J36,J37,J38,J39,J41,J42,J43,J44,J45,J46,J47,J48,J49,J51,J52,J53,J54,J55,J56,J57,J58,J59,J61,J62,J63,J64,J65,J66,J67,J68,J69,Xu1,Xu2,Xu3,Xu4,Xu5,Xu6,Xu7,Xu8,Xu9,X01,X02,X03,X04,X05,X06,X07,X08,X09,X11,X12,X13,X14,X15,X16,X17,X18,X19,X21,X22,X23,X24,X25,X26,X27,X28,X29,X31,X32,X33,X34,X35,X36,X37,X38,X39,X41,X42,X43,X44,X45,X46,X47,X48,X49,X51,X52,X53,X54,X55,X56,X57,X58,X59,X61,X62,X63,X64,X65,X66,X67,X68,X69,ctx,C41,C42,C43,C44,C45,C46,C47,C48,C49,C51,C52,C53,C54,C55,C56,C57,C58,C59,C61,C62,C63,C64,C65,C66,C67,C68,C69")
}


# =========================================#
# 5. Plotting
# =========================================#
Bo <- F
source("functions/step_1.R")
source("functions/step_2.R")
