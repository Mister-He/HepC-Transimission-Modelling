# Model scripts guide

This folder contains `sim.cpp`, model setup/data preprocessing, and rate
inputs. The C++ simulator implements the five HCV states and 150-year
equilibrium convention described in `docs/AGENTS.md`.

- `sim.cpp`: the ODE right-hand side and RK4 solver. The first-arrest rate is
  derived as `lambda1 = lambda3 * c_composite / beta_scale`.
- `setup.R`: parameters, initial state and the inflow schedule (source).
- `preprocess.R`, `release_rate.stan`, `rearrest_rate.stan`: clean the prison
  records (`data/sps.csv`, not in the repository) and fit release and re-arrest
  rates; they write `output/release_rate_fit.rds` and `output/rearrest_rate_fit.rds`.
- `update_rate_params.R`: copies those fits into the `lambda2`, `lambda3`,
  `pi_recid` and `c_composite` entries of `setup.R` after validating them.
