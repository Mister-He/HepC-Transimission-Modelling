# Training

## Deterministic calibration

```bash
Rscript scripts/run_calibration.R \
  --run-id r_20260929_agepop --maxit 1000 --n-starts 5 --seed 20261001
```

Targets are read from `data/calibration_targets.json`. Every objective
evaluation runs model time 0--150, checks the final five-year equilibrium, and
uses the final state as the 2015 snapshot. The likelihood contains only six
binomial seroprevalence terms and six log-normal population terms.

The six fitted contact scales and six age-specific inflow scales multiply the
2015 SingStat/CNB baseline in `data/singstat_new_drug_abusers.csv`. Values from
2016 onward use the same fitted scales; the latest observed year is held
constant through 2030.

## Neural posterior estimation

```bash
Rscript scripts/run_npe.R --step all \
  --fit output/calibration/r_20260929_agepop/fit.rds \
  --out-dir output/calibration/npe_agepop_recal \
  --n-sims 3000 --n-draws 5000 --n-cores 4 --n-rounds 1 \
  --max-epochs 500 --seed 20261001 \
  --npe-contact-sd 0.20 --npe-beta-sd 0.20 \
  --npe-hidden-features 128 --npe-num-transforms 10 \
  --python /path/to/python
```

NPE uses binomial prevalence noise and log-normal population noise
(`sigma = 0.10`). SBC coverage and posterior spread are quality gates. No
Laplace approximation is run or stored. The corrected-target accepted run is
in `output/calibration/npe_agepop_recal/`; the R prior simulator and Python
prior now use the same SD. The current machine uses
`/opt/homebrew/Caskroom/miniconda/base/envs/bayes-venv/bin/python`.

## MCMC validation and 2015 predictions

```bash
Rscript scripts/run_mcmc_validation.R \
  --fit output/calibration/r_20260929_agepop/fit.rds \
  --npe-dir output/calibration/npe_agepop_recal \
  --out-dir output/calibration/mcmc_recal \
  --n-iter 120000 --burnin 20000 --thin 1 \
  --n-chains 4 --parallel-chains 4 --seed 20261002 \
  --contact-sd 0.20 --beta-sd 0.20 --dt 0.0833 \
  --ess-min 4000 --rhat-max 1.001 --same-start true
Rscript scripts/predict_2015_posterior.R \
  --fit output/calibration/r_20260929_agepop/fit.rds \
  --out-dir output/calibration/npe_agepop_recal --year 2015 --cores 4
Rscript scripts/plot_calibration.R \
  --fit output/calibration/r_20260929_agepop/fit.rds \
  --npe-dir output/calibration/npe_agepop_recal \
  --mcmc-dir output/calibration/mcmc_recal --cores 4 --max-draws 2000
```

The MCMC chain uses a monthly integration step over the full 150-year
equilibrium horizon. Its strict gate requires every parameter to have pooled
ESS above 4,000 and R-hat below 1.001. Run the calibration plot after MCMC
finishes so it uses the new samples. The 2015 stratum prediction uses both NPE
seed posteriors; the calibration comparison randomly selects 2,000 draws per
method.

## 2026 treatment scenarios

```bash
Rscript scripts/scenarios_2026.R \
  --fit output/calibration/r_20260929_agepop/fit.rds \
  --posterior-dir output/calibration/npe_agepop_recal
```

The script retains the accepted NPE posterior median and 2015 equilibrium, runs
J background treatment from 2019, and branches all policies from the same 2026
state. D has no treatment; F/X interventions cover all age groups. Annual
screening quotas are 600, 1,000, 2,000, 5,000 or 8,000, with F receiving 25%, 50% or 75%,
plus J-only and infinite-screening comparators. Quotas generate chronic
start capacities from a fixed year-start non-treatment population snapshot;
the C++ solver enforces them at the actual c-to-t event. See
`data/scenario_policy_config.json` for externally editable assumptions and
`output/scenarios_2026/README.md` for definitions, cost evidence and results.
Costs are illustrative SGD public-price proxies, not confirmed Singapore
healthcare-payer net costs. The outputs are expected course events under one
NPE median parameter vector, not unique patients or posterior policy
probabilities.
