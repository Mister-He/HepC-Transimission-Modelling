# Model scripts guide

This folder contains `sim.cpp`, model setup/data preprocessing, and rate
inputs. The C++ simulator implements the five HCV states and 150-year
equilibrium convention described in `docs/AGENTS.md`.

- `sim.cpp`: the ODE right-hand side and RK4 solver.
  - **First-arrest rate:** `lambda1` is derived from the Stan first-arrest flow `f_first`
    through the HCV-free steady state of D:
    `D_i = (beta_i + 0.1*D_{i-1} - f_i) / (mu_i*omega + a_i)`, `lambda1_i = f_i / D_i`.
    `first_arrest_rate()` exposes it to R. An entry flow too small to supply `f_first`
    stops the run (infeasible parameter set).
  - **Viraemic excess mortality:** chronic viraemic people (state c) die at `hr_viraemic`
    times the background rate. `omega` is rescaled by `1 + (hr - 1) * viraemic_share_ref`.
- `setup.R`: parameters, initial state and the inflow schedule (source).
- `preprocess.R`, `release_rate.stan`, `rearrest_rate.stan`: clean the prison
  records (`data/sps.csv`, not tracked by git; people ever held for drug consumption)
  and fit release and re-arrest rates.
  - The re-arrest model's 2015 prison steady state uses the ODE's ageing rate (0.1)
    and mortality (`mu*omega`), and gives the first-arrest flow `f_pop`.
  - Run with `Rscript scripts/model/preprocess.R` from the repository root. It writes
    `output/release_rate_fit.rds` and `output/rearrest_rate_fit.rds`.
- `update_rate_params.R`: copies those fits into the `lambda2`, `lambda3`,
  `pi_recid` and `f_first` entries of `setup.R` after validating them.
