# Model scripts guide

This folder contains `sim.cpp`, model setup/data preprocessing, and rate
inputs. The C++ simulator implements the five HCV states and 150-year
equilibrium convention described in `docs/AGENTS.md`.

- `sim.cpp`: the ODE right-hand side and RK4 solver.
  - **First-arrest rate:** D -> J at `c_first * lambda3`. The prison records never observe
    the never-arrested, so `c_first` is calibrated (θ7 = log c).
  - **Mortality (one route per cause, 2026-11-03):** non-liver deaths are the general-population
    rate `mu` (60+ = ages 60–74) times `smr_nonviraemic` = 4.73 (HCV-uninfected PWID, van Santen
    2014) for u, r, and `smr_viraemic` = 4.73 × 1.47 (extrahepatic HR, Lee 2012) for a, c, t.
    Liver deaths come only from `mu_DC`/`mu_HCC` in DC and HCC. The WHO outcome counts liver deaths.
- `setup.R`: parameters and initial state (source). Entry into D is the CNB 2015 new drug
  abuser count, unscaled, in every model year (no year-by-year schedule).
- `preprocess.R`, `release_rate.stan`, `rearrest_rate.stan`: clean the prison
  records (`data/sps.csv`, not tracked by git; people ever held for drug consumption)
  and fit release and re-arrest rates (cure-type survival model for `pi` and `lambda3`).
  - Run with `Rscript scripts/model/preprocess.R` from the repository root. It writes
    `output/release_rate_fit.rds` and `output/rearrest_rate_fit.rds`.
- `update_rate_params.R`: copies those fits into the `lambda2`, `lambda3` and `pi_recid`
  entries of `setup.R` after validating them.
