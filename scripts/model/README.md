# Model scripts guide

This folder contains `sim.cpp`, model setup/data preprocessing, and rate
inputs. The C++ simulator implements the five HCV states and 150-year
equilibrium convention described in `docs/AGENTS.md`.

- `sim.cpp`: the ODE right-hand side and RK4 solver.
  - **First-arrest rate:** D -> J at `c_first * lambda3`. The prison records never observe
    the never-arrested, so `c_first` is calibrated (θ7 = log c).
  - **Mortality:** background deaths are the general-population rate `mu` times a PWID SMR:
    `smr_viraemic` for states a, c, t and `smr_nonviraemic` for u, r (van Santen 2014). DC and
    HCC add liver deaths on top. The viraemic share follows the model state over time.
- `setup.R`: parameters, initial state and the inflow schedule (source). Entry into D is the
  CNB new drug abuser count, unscaled.
- `preprocess.R`, `release_rate.stan`, `rearrest_rate.stan`: clean the prison
  records (`data/sps.csv`, not tracked by git; people ever held for drug consumption)
  and fit release and re-arrest rates (cure-type survival model for `pi` and `lambda3`).
  - Run with `Rscript scripts/model/preprocess.R` from the repository root. It writes
    `output/release_rate_fit.rds` and `output/rearrest_rate_fit.rds`.
- `update_rate_params.R`: copies those fits into the `lambda2`, `lambda3` and `pi_recid`
  entries of `setup.R` after validating them.
