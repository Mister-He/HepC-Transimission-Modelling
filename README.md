# Hepatitis C Transmission and Progression Among People Who Have Injected Drugs in Singapore: Modelling Treatment for Eradication

## Overview

This repository contains a comprehensive mathematical model for simulating Hepatitis C virus (HCV) transmission and progression among people who inject drugs (PWID) in Singapore. The model uses a compartmental approach with age-stratified populations to evaluate various intervention strategies for HCV eradication.

## Summary

The study provides a detailed analysis of HCV infection trends among different populations:

- **Modeling Approach:** We used a combination of compartmental and stochastic models to simulate infection dynamics across age groups.
- **Population Groups:** Analysis covered current PWID, prisoners, and former PWID, both individually and in combined scenarios.
- **Intervention Scenarios:** Various intervention strategies were modeled, showing the impact of different annual reduction rates on infection outcomes.
- **Statistical Methods:** Bayesian Markov Chain Monte Carlo (MCMC) methods were employed to calibrate age-specific contact matrices and quantify uncertainty.
- **Results Interpretation:** The results demonstrate significant reductions in infections with increased intervention efforts, highlighting the importance of targeted strategies for different groups.
- **Policy Implications:** The findings suggest that comprehensive, combined interventions are the most effective in reducing overall infection rates.

## Publication

[[Preprint]](https://doi.org/10.1101/2025.10.24.25338708)

## Project Structure

```
HepC-Transmission-Modelling/
├── main.R                      # Main analysis pipeline
├── README.md                   # This file
├── functions/                  # Core functions
│   ├── utils.r                # Model implementation and utilities
│   ├── step_1.R               # Visualization: infections by strategy
│   ├── step_2.R               # Visualization: complications by strategy
│   └── step_3.R               # Additional analyses and validation
├── settings/                   # Configuration files
│   ├── Final_strategies.csv   # Treatment intervention strategies
│   ├── graph_config.csv       # Plotting configurations
│   └── table*.csv             # Scenario-specific parameters
└── output/                     # Generated results
    ├── equilibrium_*.rdata    # Equilibrium model states
    ├── mcmc/                  # MCMC calibration results
    ├── fig/                   # Generated figures
    └── BI_normal_SVR/         # Bayesian credible intervals
```

## Requirements

### Software Requirements
- **R version:** ≥ 4.5.1 (tested on 4.5.1)
- **Operating System:** macOS, Linux, or Windows

### Required R Packages

The following packages will be automatically loaded by `main.R`:

```r
# Core packages
data.table
tictoc
dplyr
ggplot2

# MCMC and Bayesian analysis
mcmc
coda
posterior
bayesplot

# Parallel computing
parallel
future
future.apply

# Visualization
grid
```

### Installing Dependencies

Install all required packages using:

```r
install.packages(c("data.table", "tictoc", "dplyr", "ggplot2", 
                   "mcmc", "coda", "posterior", "bayesplot",
                   "parallel", "future", "future.apply"))
```

## How to Run

### Complete Analysis Pipeline

To reproduce all outputs from the model, run the main analysis script:

**Option 1: Using RStudio or R console**
```r
source("main.R")
```

**Option 2: Command line (foreground)**
```bash
Rscript main.R
```

**Option 3: Command line (background, with logging)**
```bash
nohup Rscript main.R >output.log 2>&1 &
```

### What `main.R` Does

The main analysis script performs the following steps:

1. **Data Preparation**
   - Loads observed HCV prevalence data by age group
   - Reads intervention strategy configurations
   - Initializes model parameters and initial states

2. **MCMC Calibration** (Steps 1-3)
   - Calibrates age-specific contact scaling factors
   - Runs 4 parallel MCMC chains with 100,000 samples each
   - Performs convergence diagnostics (Gelman-Rubin R̂ and ESS)
   - Saves results to `output/mcmc/`

3. **MCMC Visualization** (Step 4)
   - Generates trace plots, density overlays, ACF plots, and rank histograms
   - Saves diagnostic plots to `output/fig/mcmc_bayesplot/`

4. **Post-MCMC Analysis** (Step 5)
   - Computes equilibrium states with calibrated parameters
   - Evaluates all intervention strategies with posterior means
   - Generates Bayesian credible intervals (100 posterior samples)
   - Performs posterior predictive checks

5. **Visualization** (Step 6)
   - Generates publication-ready figures
   - Creates comparison plots for different strategies

### Running Individual Analysis Steps

You can also run specific visualization or analysis functions:

```r
# Step 1: Generate infection trajectory plots
source("functions/step_1.R")

# Step 2: Generate complication trajectory plots
source("functions/step_2.R")

# Step 3: Additional validation analyses
source("functions/step_3.R")
```

## Key Model Features

### Age-Stratified Compartmental Model
- **9 age groups:** <20, 20-24, 25-29, 30-34, 35-39, 40-44, 45-49, 50-54, 55+
- **Population compartments:** Current PWID (D), prisoners (J), former PWID (X)
- **Disease stages:** Uninfected, acute (stages 0-3), chronic (stages 4-6)

### Intervention Strategies
The model evaluates multiple treatment scenarios defined in `settings/Final_strategies.csv`:
- Varying treatment coverage rates for different populations
- Different sustained virologic response (SVR) rates
- Combined intervention strategies

### Bayesian Calibration
- Uses MCMC to estimate age-specific contact matrix scaling factors
- Fits model to observed age-stratified HCV prevalence data
- Quantifies uncertainty through posterior distributions

## Expected Runtime

- **Full MCMC analysis:** 30-60 minutes (depending on CPU cores)
- **Post-MCMC simulations:** 10-20 minutes
- **Bayesian credible intervals (100 samples):** 1-2 hours
- **Total pipeline:** ~2-3 hours on a modern multi-core system

## Output Files

### MCMC Results
- `output/mcmc/mcmc_raw_results.rds` - Raw MCMC chains
- `output/mcmc/mcmc_results.rds` - Processed results with diagnostics
- `output/mcmc/mcmc_bci.rds` - Bayesian credible intervals for prevalence

### Model States
- `output/equilibrium_*.rdata` - Equilibrium model states for different scenarios

### Intervention Results
- `settings/tableHIGH_1_*.csv` - Detailed simulation results by strategy
- `output/BI_normal_SVR/sim_data_*.rds` - Posterior predictive simulations

### Figures
- `output/fig/mcmc_bayesplot/` - MCMC diagnostic plots
- Additional publication figures generated by step_1.R and step_2.R

## Configuration

### Modifying Intervention Strategies

Edit `settings/Final_strategies.csv` to define custom intervention scenarios. Key parameters:
- `pD` - Annual treatment proportion for current PWID
- `pJ1`, `pJ2`, `pJ3` - Treatment proportions for different prisoner stages
- `pX` - Annual treatment proportion for former PWID
- `strategy_class` - Strategy grouping for visualization

### Model Parameters

Core model parameters are defined in `functions/utils.r` within the `set_parameters()` function:
- Contact rates and mixing patterns
- Disease progression rates
- Treatment efficacy (SVR rates)
- Recidivism probabilities

## Troubleshooting

### Memory Issues
If you encounter memory errors, reduce the number of posterior samples:
```r
# In main.R, section 5.3, change:
for (iteration in 1:100) 
# to a smaller number like:
for (iteration in 1:50)
```

### Convergence Issues
If MCMC chains don't converge (R̂ > 1.1), try:
- Increasing burn-in period: `n_burnin <- 2e+4`
- Running longer chains: `n_samples <- 2e+5`
- Adjusting the proposal scale in the pilot run

### Parallel Computing
The number of parallel chains can be adjusted:
```r
n_chains <- 4  # Reduce if you have fewer CPU cores
```