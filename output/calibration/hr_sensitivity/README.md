# HR sensitivity of the Nelder–Mead fit

These are sensitivity fits only; nothing here feeds the production path.
- **Command:** `scripts/run_calibration.R --hr-viraemic <HR> --n-starts 2 --maxit 3000 --nm-restarts 5
  --sd-perturb 0.9 --seed 20261001`, which gives 3 starts.
- **Model:** `lambda1 = lambda3`, free entry scales (no entry anchor), `viraemic_share_ref = 0.205`.
- **Files kept:** `predictions.csv`, `solutions.csv`, `run_config.json` and `equilibrium.csv`.

| HR | Source | 60+ seroprevalence (target 0.429, SE 0.041) | Objective |
| --- | --- | --- | --- |
| 1 | no excess mortality | 0.440 | 18.91 |
| 1.11 | Piselli et al. 2021, non-liver, RNA+ | 0.436 | 18.79 |
| 1.23 | Lee et al. 2012, lower 95% limit | 0.431 | 18.74 |
| 1.47 | Lee et al. 2012 (production, `../r_20261020_agepop/`, 8 starts) | 0.429 | 18.73 |
| 1.77 | Lee et al. 2012, upper 95% limit | 0.429 | 18.73 |

**Reading:** with free entry scales, 60+ fits within 0.3 SE of its target at every HR. The
viraemic excess mortality improves the fit only slightly, so the HR is not identified by the
prevalence data; its effect on the scenarios comes through the deaths it adds, not the fit.
