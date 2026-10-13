# Scripts guide

Executable model and analysis code. Each script starts with a short description.

- **`model/`**: the simulator (`sim.cpp`), parameters and initial state (`setup.R`), prison-record
  preprocessing (`preprocess.R`) and its Stan models. The top-level `preprocess.R` and `*.stan` files
  are an older copy that uses `../` paths.
- **`calibration/`**: shared functions:
  - targets, J-stratum summaries and the equilibrium check;
  - the likelihood;
  - Nelder–Mead (`calibrate_nm.R`);
  - prior sampling, by LHS or at random (`prior_sampling.R`);
  - the exact log posterior and adaptive Metropolis sampler (`mcmc.R`).

Workflow, in order:

1. **Prior centre:** `run_calibration.R`, a Nelder–Mead multi-start fit run in parallel with `--cores`.
2. **Inference:**
   - `run_mcmc_reference.R`: reference MCMC, the production posterior.
   - `run_npe.R` and `npe_train.py`: single-round NPE trained on Latin hypercube prior simulations.
3. **Checks:**
   - `compare_npe_mcmc.R`: NPE versus reference MCMC, for parameters and policy outcomes.
   - `predict_2015_posterior.R`: 2015 predictions by stratum.
   - `plot_calibration.R`: the calibration figure.
4. **Scenarios:** `scenarios_2026.R` runs every policy for each posterior draw.
   `post_calibration_analysis.R` is a thin entry point that runs it.
5. **Manuscript:**
   - `paper_figures.R`: the manuscript figures.
   - `build_paper.py`: fills every number from the result files and runs pandoc to write
     `docs/main.docx` and `docs/Supplementary Information.docx`.
6. **Tests:** `run_tests.R` runs everything under `tests/`.

