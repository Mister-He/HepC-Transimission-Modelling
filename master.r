library(data.table)
library(tictoc)
library(dplyr)
library(ggplot2)

source("functions/utilis.r")
strategies = read.csv("settings/Final_strategies.csv" )

############################################################################

# Change annual proportions to monthly proportions
strategies$pD = strategies$pD/12  
strategies$pJ1 = strategies$pJ1/12
strategies$pJ2 = strategies$pJ2/12
strategies$pJ3 = strategies$pJ3/12
strategies$pX = strategies$pX/12

FIRSTTIME = F

#---- initialise ----#
# (only do this once)
if (FIRSTTIME){
  intervention = strategies[1,]
  configuration = create_configuration(ny=600) #250
  parameters = set_parameters(nsteps_permonth=configuration$nsteps_permonth)
  initialstates = make_initial_states(configuration,equilibrium = F,parameters)
  tic(msg = "Single epimodel round")
  Y = epimodel(initialstates,parameters,configuration,intervention)
  X = epimodel_fast(initialstates,parameters,configuration,intervention)
  compare(X,Y)
  #[1] "The lists are exactly the same."
  # So without intervention, the two models are the same
  toc()
  
  print(final_prevalence(X)) 
  store_states(X)
}

#prison data (age-specific for 9 age groups)
obs_data = list(pos=c(55, 145, 183, 164, 212, 299, 222, 190, 133),
                tot=c(307, 797, 829, 633, 598, 642, 481, 439, 366))

# set up to calculate logLL, which can be used for MCMC
nsteps_permonth = 1 # 4 # 30
intervention = strategies[1,]
configuration = create_configuration(ny=600, nsteps_permonth=nsteps_permonth) #250
parameters = set_parameters(nsteps_permonth=configuration$nsteps_permonth)
initialstates = make_initial_states(configuration,equilibrium = FALSE,parameters)

# Get this from some personal guess
xparam_init = c(9.0,0.2,0.18,0.13,0.44,0.3,0.08,0.02,0.01, parameters$alpha)

## MCMC
# logLL fun for adaptive MCMC
LB = c(8.5, 0,   0,   0,   0,   0.1, 0,    0,   0,    rep(0.001,9))
UB = c(9.5, 0.5, 0.5, 0.3, 0.5, 0.5, 0.2,  0.1, 0.02, rep(0.04,9))  #10*rep(1, n_xparam)
labs = c(paste("<",20,  sep = ""), paste(seq(20, 50, by = 5), seq(20 + 5 - 1, 55 - 1, by = 5), sep = "-"), paste(55, "+", sep = ""))

#=========================================#
# MCMC procedure
library(mcmc)
library(parallel)
library(coda)

getmode = function(v) {
  uniqv = unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

n_burnin = 2e+3
n_MCMC_base = 2e+4
n_MCMC = n_MCMC_base + n_burnin # MCMC iterations, include burn-in
n_thinning = 50

# # New MCMC method
# # grid search scale = 0.01, 0.005, 0.001
# for (scale in c(0.01, 0.005, 0.001)) {
#   time1 = Sys.time()
#   out = metrop(obj     = fun_logLL_MCMC,
#                initial = xparam_init,
#                nbatch  = 50,
#                blen    = 1000,
#                scale   = scale)
#   time2 = Sys.time()
#   
#   # save output
#   print(paste0('For scale = ', scale))
#   print(list('Accept runs' = out$accept, 
#              'Time took'   = time2 - time1))
# }

# According to the grid search result, smaller the scale,
# the better the acceptance rate, but the longer the time it takes,
# so we choose 0.005 as the scale, whose acceptance is about 0.14
time1 = Sys.time()
out = metrop(obj     = fun_logLL_MCMC,
             initial = xparam_init,
             nbatch  = 20,
             blen    = 1000,
             scale   = 0.05 * xparam_init)
             
run_chain <- function(seed, init) {
  set.seed(seed)
  fit <- metrop(
    obj     = fun_logLL_MCMC,
    initial = init,
    nbatch  = 50,
    blen    = 1000,
    scale   = 0.05 * init
  )
  return(fit)
}

inits <- replicate(4, xparam_init * runif(length(xparam_init), 0.8, 1.2), simplify = FALSE)
seeds <- 1:4
fits <- mclapply(seq_along(seeds), function(k) run_chain(seeds[k], inits[[k]]), mc.cores = 4)

# Convert to coda
mcmc_list <- mcmc.list(lapply(fits, function(f) mcmc(f$batch)))
gelman.diag(mcmc_list) # R-hat
effectiveSize(mcmc_list) # ESS
time2 = Sys.time()
print(time2 - time1)

#out = readRDS('output/MCMC_out_normal_SVR.rds')
out_pos_mean = apply(out$batch, 2, mean)
out_pos_sd   = apply(out$batch, 2, sd)

# Generate Bayesian credible intervals
samples_matrix <- do.call(rbind, mcmc_list)
n_samples <- 100
pred_prev <- matrix(NA, nrow = n_samples, ncol = 9)
for (i in 1:n_samples) {
  pred_prev[i, ] = fun_final_prev(pmax(samples_matrix[i, ],0))$prevalence_byage
  cat('Done for sample', i, 'of', n_samples, '\n')
}

# Calculate the median, 0.025 and 0.975 quantiles for each age group
pred_prev_median <- apply(pred_prev, 2, median)
pred_prev_lower  <- apply(pred_prev, 2, quantile, probs = 0.025)
pred_prev_upper  <- apply(pred_prev, 2, quantile, probs = 0.975)
fun_plotprev(pred_prev_median, obs_data)
fun_plotprev(pred_prev_lower, obs_data)
fun_plotprev(pred_prev_upper, obs_data)

#out_pos_mean = readRDS('output/MCMC_out_normal_SVR.rds')
out_final_prev = fun_final_prev(out_pos_mean)
fun_plotprev(out_final_prev$prevalence_byage, obs_data)

# Save the equilibrium status
nsteps_permonth = 1
intervention = strategies[1,]
configuration = create_configuration(ny=600, nsteps_permonth=nsteps_permonth)
parameters = set_parameters(nsteps_permonth=configuration$nsteps_permonth)
initialstates = make_initial_states(configuration,equilibrium = FALSE,parameters)
X = epimodel_fast(initialstates,parameters,
                  configuration,intervention, 
                  adjContactMat = out_pos_mean[1:9], 
                  alpha         = out_pos_mean[10:18])
final_prevalence(X)
store_states(X, outfile='output/equilibrium_241410_low_SVR.rdata')
#=========================================#

################################## for LATER (after fitting ONLY) ############################
#---- interventions ----#
outputfile1 = 'settings/tableHIGH_1_wFINALstrategies_yichen_low_SVR.csv'
tabulator(file=outputfile1,header=TRUE)
for(str in 1:dim(strategies)[1]){
  cat('Running strategy',str,'of',dim(strategies)[1],'\n')
  intervention = strategies[str,]
  configuration = create_configuration(50)
  parameters = set_parameters(nsteps_permonth=configuration$nsteps_permonth)
  parameters$alpha = out_pos_mean[10:18] / configuration$nsteps_permonth
  parameters$theta_1 = parameters$theta_1 * 0.8
  parameters$theta_2 = parameters$theta_2 * 0.8
  initialstates = make_initial_states(configuration,equilibrium = TRUE,parameters)
  X = epimodel(initialstates,parameters,
               configuration,intervention,
               adjContactMat = out_pos_mean[1:9])
  tabulator(X,intervention,file=outputfile1,resolution='HIGH')
  print(final_prevalence(X))
}

#---- interventions with Bayesian Credible Intervals ----#
# Sample parameters from MCMC posterior parameters 50 times
out = readRDS('output/MCMC_out_normal_SVR.rds')

for (iteration in 1) {
  # Sample once from the 50*18 posterior distribution matrix
  pos_sample = apply(out$batch,2,function(x) sample(x,1))
  
  # Create a list to save all results in this iteration
  result = list()
  
  for(str in 1:2){
    cat('Running strategy',str,'of',dim(strategies)[1],'in iteration',iteration,'\n')
    intervention = strategies[str,]
    configuration = create_configuration(50)
    parameters = set_parameters(nsteps_permonth=configuration$nsteps_permonth)
    parameters$alpha = pos_sample[10:18] / configuration$nsteps_permonth
    # parameters$theta_1 = parameters$theta_1 * 0.8
    # parameters$theta_2 = parameters$theta_2 * 0.8
    initialstates = make_initial_states(configuration,equilibrium = TRUE,parameters)
    X = epimodel(initialstates,parameters,
                 configuration,intervention,
                 adjContactMat = pos_sample[1:9])
    print(final_prevalence(X))
    result[[str]] = rdstabulator(X,intervention)
  }
  sim_data = do.call(cbind,result)
  colnames(simdata) = c('Strategy,t,Du1,Du2,Du3,Du4,Du5,Du6,Du7,Du8,Du9,D01,D02,D03,D04,D05,D06,D07,D08,D09,D11,D12,D13,D14,D15,D16,D17,D18,D19,D21,D22,D23,D24,D25,D26,D27,D28,D29,D31,D32,D33,D34,D35,D36,D37,D38,D39,D41,D42,D43,D44,D45,D46,D47,D48,D49,D51,D52,D53,D54,D55,D56,D57,D58,D59,D61,D62,D63,D64,D65,D66,D67,D68,D69,Ju1,Ju2,Ju3,Ju4,Ju5,Ju6,Ju7,Ju8,Ju9,J01,J02,J03,J04,J05,J06,J07,J08,J09,J11,J12,J13,J14,J15,J16,J17,J18,J19,J21,J22,J23,J24,J25,J26,J27,J28,J29,J31,J32,J33,J34,J35,J36,J37,J38,J39,J41,J42,J43,J44,J45,J46,J47,J48,J49,J51,J52,J53,J54,J55,J56,J57,J58,J59,J61,J62,J63,J64,J65,J66,J67,J68,J69,Xu1,Xu2,Xu3,Xu4,Xu5,Xu6,Xu7,Xu8,Xu9,X01,X02,X03,X04,X05,X06,X07,X08,X09,X11,X12,X13,X14,X15,X16,X17,X18,X19,X21,X22,X23,X24,X25,X26,X27,X28,X29,X31,X32,X33,X34,X35,X36,X37,X38,X39,X41,X42,X43,X44,X45,X46,X47,X48,X49,X51,X52,X53,X54,X55,X56,X57,X58,X59,X61,X62,X63,X64,X65,X66,X67,X68,X69,ctx,C41,C42,C43,C44,C45,C46,C47,C48,C49,C51,C52,C53,C54,C55,C56,C57,C58,C59,C61,C62,C63,C64,C65,C66,C67,C68,C69')
}

#=========================================#
# Plotting
Bo = F
source('functions/step_4.R')
source('functions/step_5.R')
#=========================================#