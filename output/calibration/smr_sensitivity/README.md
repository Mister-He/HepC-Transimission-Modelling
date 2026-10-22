# SMR sensitivity fits

Nelder–Mead refits (3 uninformative starts each, maxit 3000, `--nm-restarts 5`, seed
20261001) with the PWID standardised mortality ratios at the bounds of their 95% CIs
(van Santen et al. 2014). The production fit (`../r_20261021_agepop`) uses 9.52 (viraemic)
and 4.73 (non-viraemic). Each folder keeps `predictions.csv`, `solutions.csv`,
`run_config.json` and `targets.json`.

| Run | SMR viraemic | SMR non-viraemic | Objective (−log posterior) | c | Seroprevalence RMSE |
| --- | --- | --- | --- | --- | --- |
| `smr_low` | 8.06 | 3.47 | −7.847 | 1.100 | 0.0119 |
| production | 9.52 | 4.73 | −7.829 | 1.099 | 0.0130 |
| `smr_high` | 11.26 | 6.45 | −7.787 | 1.097 | 0.0148 |

The SMRs barely change the fit or `c`; all starts converged and passed the equilibrium check.

Rerun:

```bash
Rscript scripts/run_calibration.R --out-dir output/calibration/smr_sensitivity --run-id smr_low \
  --smr-viraemic 8.06 --smr-nonviraemic 3.47 --maxit 3000 --nm-restarts 5 --n-starts 3 \
  --sd-perturb 0.9 --seed 20261001 --cores 3
# smr_high: --smr-viraemic 11.26 --smr-nonviraemic 6.45
```
