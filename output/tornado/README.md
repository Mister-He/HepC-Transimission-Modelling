# One-way sensitivity analysis (tornado)

Written by `Rscript scripts/tornado.R --root . --cores 14` (2026-11-03 review, `REVIEW_REPORT.md` R-05).
The 13 calibrated parameters are held at the Nelder–Mead best fit (`output/calibration/r_20261103_agepop`);
each fixed input is set to the low and then the high end of its range, without refitting, and
`scripts/scenarios_2026.R --param-file` runs all 17 policies for that single parameter set.

- `ranges.csv`: every input varied, its base, low and high values and the source of the range (the
  Notation table of `Hepc_modelling_technical_details.pdf`, otherwise the cited literature or the Stan
  rate fits).
- `best_theta.csv`: the 13 calibrated parameters used in every run.
- `tornado_results.csv`: 2030 reductions (incidence, liver deaths, liver + extrahepatic deaths) and
  pass/fail for every run and policy; `run = base` is the result at the best fit.
- `pass_flips.csv`: runs in which a policy's pass/fail differs from the base run.
- `failed_runs.csv`: runs whose 2015 snapshot was not at equilibrium after the 150-year burn-in.
- `runs/<input>_<low|high>/`: the full scenario output, the override (`override.rds`) and the log.
- `tornado.png/.svg`: drawn by `python scripts/figures.py . --only tornado` (S1 Text Fig G).
