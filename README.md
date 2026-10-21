# Repository guide

This repository contains the Singapore PWID HCV transmission model. Source
model code is in `scripts/model/`, calibration code in `scripts/calibration/`,
tests in `tests/`, inputs in `data/`, documentation in `docs/`, and generated
results in `output/` and `figures/`. Read `AGENTS.md` and the relevant child
directory guide before traversing files. Current 2026–2030 policy methods,
results and cost evidence are in `output/scenarios_2026/README.md`. The calibration path
(Nelder–Mead posterior mode → reference MCMC → scenarios; NPE + SBC as cross-validation only)
is in `docs/training.md`;
`回答.md` preserves interpretation of an earlier scenario design.
