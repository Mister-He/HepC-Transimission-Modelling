# Scripts guide

Executable model and analysis code lives here. `model/` contains the simulator,
`calibration/` contains reusable calibration functions, and top-level scripts
run calibration, NPE, MCMC validation, scenarios, tests, and plotting.

Manuscript pipeline: `scenarios_2026.R` (posterior-draw scenarios, also writes the
no-treatment counterfactual) → `paper_figures.R` (writes `paper_fig*.png` into
`output/scenarios_2026/`) → `build_paper.py` (fills all numbers from the result CSVs
and runs pandoc to write `docs/main.docx` and `docs/Supplementary Information.docx`).
