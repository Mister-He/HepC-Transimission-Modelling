# Model scripts guide

This folder contains `sim.cpp`, model setup/data preprocessing, and rate
inputs. The C++ simulator implements the five HCV states and 150-year
equilibrium convention described in `docs/AGENTS.md`.

- `sim.cpp`: the ODE right-hand side and RK4 solver.
  - **First-arrest rate:** `lambda1` (D -> J) is read from `setup.R`. The prison records never
    observe the never-arrested, so it is assumed equal to the re-arrest rate `lambda3`.
  - **Viraemic excess mortality:** chronic viraemic people (state c) die at `hr_viraemic`
    times the background rate. `omega` is rescaled by `1 + (hr - 1) * viraemic_share_ref`.
- `setup.R`: parameters, initial state and the inflow schedule (source).
- `preprocess.R`, `release_rate.stan`, `rearrest_rate.stan`: clean the prison
  records (`data/sps.csv`, not tracked by git; people ever held for drug consumption)
  and fit release and re-arrest rates.
  - The re-arrest model's 2015 prison steady state uses the ODE's ageing rate (0.1)
    and mortality (`mu*omega`). Its first-arrest flow `f_pop` informs `lambda3` and `pi`
    there only; it is not passed to the ODE.
  - Run with `Rscript scripts/model/preprocess.R` from the repository root. It writes
    `output/release_rate_fit.rds` and `output/rearrest_rate_fit.rds`.
- `update_rate_params.R`: copies those fits into the `lambda2`, `lambda3` and `pi_recid`
  entries of `setup.R` after validating them, and writes `lambda1` equal to `lambda3`.
