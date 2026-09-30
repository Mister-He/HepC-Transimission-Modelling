# Structured execution prompt

## Objective

Continue work on the Singapore PWID HCV transmission model on a dedicated
`feature/AddRecovery` branch and worktree, while preserving the current model
changes and keeping result paths stable.

## Repository and Git

1. Add persistent agent rules requiring minimal creation of files/folders.
2. Add a concise Markdown directory guide to every substantive existing
   directory; document generated timestamped log directories from their parent.
3. Move all current working-tree changes onto `feature/AddRecovery`.
4. Create a worktree for that branch in the parent directory of this repo.

## Calibration

1. Re-run the current NPE in the existing NPE result directory; do not create a
   new result directory.
2. Reduce the excessive NPE 95% predictive BCI width for ages 30--49 using
   defensible prior, likelihood, simulation-budget, or estimator changes.
3. Retain SBC/coverage diagnostics and report any remaining structural
   non-identifiability rather than hiding it.
4. Relax the MCMC convergence gate to ESS > 4,000 and Rhat < 1.001. Initialize
   adaptive MCMC from the NPE posterior covariance and overwrite the existing
   MCMC result files.

## Literature review

Review current and primary/authoritative evidence relevant to Singapore HCV:

- HCV treatment cost and cost-effectiveness;
- treatment access, linkage, and barriers among PWID and people in the
  Singapore care system;
- treatment capacity, service delivery, and feasible treatment volume.

Write the synthesis, evidence links, uncertainty, and model implications to
`docs/background.md`.

## Policy optimization

1. Expand treatment levels beyond the current grid where computationally and
   substantively justified.
2. Consider treatment in additional strata, subject to the model's stratum
   definitions and the user's policy scope.
3. Define a transparent objective balancing HCV benefit against treatment
   burden, access constraints, and treatment-capacity constraints.
4. Implement or prepare Bayesian optimization to identify the optimal treatment
   policy, retaining reproducible policy definitions and diagnostics in the
   existing output paths.

## Guardrails

- Do not add unrelated features or run post-calibration analyses that are not
  needed for this objective.
- Preserve the 150-year equilibrium convention and the recovered state model.
- Clearly distinguish exact-model results, approximations, literature facts,
  and assumptions.
