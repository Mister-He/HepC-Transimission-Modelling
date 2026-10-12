Please clean up the project’s inference workflows, historical outputs, and script descriptions, and update the single-round NPE sampling and training pipeline.

### 1. Retain only single-round NPE and reference MCMC

Keep only:

- Single-round NPE.
- The existing long-chain reference MCMC.

Remove code, configurations, execution entry points, results, figures, and documentation references specific to multi-round APT/SNPE-C and short-chain MCMC. Update any analysis or reporting scripts that depend on them.

Important distinctions:

- Single-round NPE may still use the `NPE_C` class from `sbi`. Do not delete it based on its name. Remove sequential proposal and APT-loss logic while retaining the implementation needed for single-round NPE.
- Short-chain and reference MCMC may share a sampler. Preserve all shared code and configurations required by reference MCMC.
- Preserve the existing reference MCMC results. Do not rerun reference MCMC merely because directories have been reorganized.

### 2. Clean up and reorganize `output/`

Inspect existing outputs and their dependencies. Remove deprecated, superseded, or obsolete results so that `output/` is easier to navigate.

Requirements:

- Clearly separate single-round NPE outputs, reference MCMC outputs, and their comparison results.
- Preserve valid calibration inputs, the reference posterior, essential run metadata, and files still required by downstream analyses.
- Determine whether results are deprecated from their generating scripts, configurations, and current dependencies—not merely their age.
- Update affected paths and documentation to prevent broken references.
- Before deletion, identify the files to be removed and briefly explain why. Retain files whose purpose remains uncertain and flag them in the final report.

### 3. Update single-round NPE to use 50,000 LHS prior simulations

Update the training pipeline to:

> Generate 50,000 parameter sets from the existing prior using Latin hypercube sampling (LHS), generate the corresponding simulated data, and train single-round NPE using the existing network and training configuration.

Requirements:

- Preserve the existing prior, parameter transformations, model, observation-noise model, constraints, and numerical solver settings.
- LHS must target the existing prior distribution. Do not replace it with uniform sampling between parameter bounds.
- Inspect the prior structure first. For independent marginal priors, generate an LHS design in the unit hypercube and transform each coordinate using the corresponding prior inverse CDF. Apply normal-prior transformations in the log-parameter space where those priors are defined.
- If the prior contains dependencies or truncation, use an appropriate transformation that preserves the intended joint prior. Do not silently replace it with independent priors.
- Set and record random seeds for reproducible sampling and training.
- Retain the existing network architecture and training hyperparameters. Do not introduce additional hyperparameter tuning in this task.
- Use the standard single-round NPE training objective. Do not introduce sequential proposals or APT loss.

Define 50,000 as the number of initially generated prior parameter sets. For failed simulations or equilibrium violations, retain the existing handling rules after verifying their correctness. Record attempted simulations, valid simulations, and exclusion reasons separately. Do not silently generate replacement samples to obtain exactly 50,000 valid training examples.

After training, update the parameter-posterior and policy-outcome comparisons against reference MCMC. Clearly label the new results as LHS-based to distinguish them from previous random-sampling results.

### 4. Add a brief description to every project code script

Read the description at the top of `helper.R` and follow its format when adding or standardizing descriptions across all project-owned code scripts.

Requirements:

- Use approximately one or two lines to explain each script’s purpose.
- Follow the language and style used in `helper.R`.
- Use the appropriate comment syntax for each programming language, preserving the required placement of shebangs and other special headers.
- Describe the script’s purpose, not its modification history. Avoid comments such as “added,” “removed,” or “updated in this change.”
- Exclude third-party code, generated code, and non-code result files.
- Do not duplicate existing descriptions that are already concise and accurate.

### 5. Validation and deliverables

Check the retained entry points, configurations, paths, and downstream dependencies. Confirm that:

- Single-round NPE and reference MCMC remain usable.
- Active workflows no longer reference removed methods or outputs.
- The LHS parameter dimensions, support, and prior transformations are correct.
- The updated NPE training and comparison outputs have been generated. Clearly report any failed or incomplete steps.
- The output structure matches the documentation.

In the final report, briefly summarize:

1. The main code and results removed or retained.
2. The reorganized `output/` structure.
3. The LHS implementation, random seeds, and valid simulation count.
4. The NPE training configuration and updated comparison results.
5. Completed validation and any unresolved issues.
