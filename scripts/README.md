# Scripts guide

Executable model and analysis code. Each script starts with a short description.

- **`model/`**: the simulator (`sim.cpp`), parameters and initial state (`setup.R`), prison-record
  preprocessing (`preprocess.R`) and its Stan models, and the rate-fit handoff into `setup.R`
  (`update_rate_params.R`).
- **`calibration/`**: shared functions:
  - targets, J-stratum summaries and the equilibrium check;
  - the likelihood;
  - Nelder–Mead (`calibrate_nm.R`);
  - prior sampling, by LHS or at random (`prior_sampling.R`; LHS for NPE training, random for SBC);
  - the exact log posterior and adaptive Metropolis sampler (`mcmc.R`).

Workflow, in order:

1. **Posterior mode:** `run_calibration.R`, a Nelder–Mead multi-start fit run in parallel with `--cores`.
2. **Production posterior:** `run_mcmc_reference.R`, the reference MCMC.
3. **Cross-validation only** (nothing here feeds the production path):
   - `run_npe.R` and `npe_train.py`: single-round NPE (5-network ensemble) trained on Latin hypercube
     prior simulations, with SBC and its coverage gate.
   - `compare_npe_mcmc.R`: NPE versus reference MCMC, for parameters and policy outcomes.
   - `sbc_mcmc.R`: a small direct SBC of the reference MCMC.
4. **Checks:**
   - `predict_2015_posterior.R`: 2015 predictions by stratum.
   - `plot_calibration.R`: the calibration figure.
5. **Scenarios:** `scenarios_2026.R` runs every policy for each posterior draw.
   `post_calibration_analysis.R` is a thin entry point that runs it.
6. **Manuscript:**
   - `figures.py`: redraws every figure (manuscript, SI, calibration diagnostics) from the saved
     outputs with the `scientific-visualization` skill, as PNG (300 dpi) plus editable SVG. Run it
     last; the older ggplot drawings in `paper_figures.R` and the other R scripts write the same
     file names.
   - `build_paper.py` (with `paper/facts.py`, `paper/main_text.py`, `paper/si_text.py`,
     `paper/registry.py`, `paper/fmt.py`): builds the PLOS Medicine manuscript, S1 Text, S1 Checklist,
     TIFF figures and the evidence registries in `docs/paper/` from the saved outputs.
7. **Tests:** `run_tests.R` runs everything under `tests/`.

