# HR sensitivity of the Nelder–Mead fit

These are sensitivity fits only; nothing here feeds the production path.
- **Command:** `scripts/run_calibration.R --hr-viraemic <HR> --n-starts 2 --maxit 3000 --nm-restarts 5
  --sd-perturb 0.9 --seed 20261001`, which gives 3 starts.
- **Viraemic share:** `viraemic_share_ref = 0.242` (an earlier fixed-point round; the production
  fit uses 0.205).
- **Files kept:** `predictions.csv`, `solutions.csv`, `run_config.json` and `equilibrium.csv`.

| HR | Source | 60+ seroprevalence (target 0.429, SE 0.041) | Objective |
| --- | --- | --- | --- |
| 1 | no excess mortality | 0.494 | 22.70 |
| 1.11 | Piselli et al. 2021, non-liver, RNA+ | 0.493 | 22.16 |
| 1.23 | Lee et al. 2012, lower 95% limit | 0.481 | 21.75 |
| 1.47 | Lee et al. 2012 (production, `../r_20261014_agepop/`, 8 starts) | 0.455 | 21.44 |
| 1.77 | Lee et al. 2012, upper 95% limit | 0.439 | 21.19 |

**Reading:** with the entry anchor, the viraemic excess mortality is what brings 60+ close to
its target. The objective improves as the HR rises.
