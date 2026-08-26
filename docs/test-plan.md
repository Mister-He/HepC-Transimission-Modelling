# Test plan

## 1. Objectives

Verify, on every change, that:

1. the model compiles and produces valid dynamics (no negative values, no
   NaN/Inf, finite positive population, prevalence in [0, 1]);
2. the calibration parameterisation and likelihoods behave correctly;
3. the equilibrium gate correctly separates converged from unconverged
   trajectories;
4. the scenario machinery correctly targets strata and age groups
   (regression protection for the `tau_stratum` / `tau_min_age` wiring);
5. the end-to-end pipeline runs without external secrets or heavy compute.

## 2. Test levels

### Unit tests — `tests/unit/`

Fast, dependency-light, no C++ compilation required.

| File | Covers |
|---|---|
| `test_targets.R` | age-group structure; prevalence inside (0,1); integer prison totals; Binomial count reconstruction within rounding tolerance |
| `test_likelihood.R` | 12-parameter exp parameterisation; bounds; Binomial/log-Normal likelihoods; plausibility pattern penalty |
| `test_metrics.R` | 384-compartment indexing; J summary extraction on synthetic outputs; fit metrics at observed values |
| `test_equilibrium.R` | equilibrium gate passes on stable trajectories and fails on perturbed ones; default criteria |

### Integration tests — `tests/integration/`

Compile `src/sim.cpp` and exercise the real pipeline.

| File | Covers |
|---|---|
| `test_sim_and_pipeline.R` | compile + output shape; no negatives/NaN; population invariants; determinism; objective finite at theta=0; scenario CSV loaded; prison-only/age-restricted/community treatment effects (regression for the strata/age wiring) |

## 3. Running tests

```bash
# all tests (unit then integration)
Rscript scripts/run_tests.R

# a single scope
Rscript scripts/run_tests.R --scope unit
Rscript scripts/run_tests.R --scope integration

# a single file
Rscript tests/unit/test_likelihood.R
Rscript tests/integration/test_sim_and_pipeline.R
```

Each test file runs in its own R session; the runner exits non-zero if any
file fails. Integration tests require `Rcpp` and `RcppArmadillo`.

## 4. CI mapping

`.github/workflows/ci.yml` runs on every push/PR:

1. `r-tests` (ubuntu): set up R, install R packages, compile `src/sim.cpp`,
   run unit tests, run integration tests;
2. `python-npe-syntax`: syntax-check `scripts/npe_train.py` (heavy NPE
   training is deliberately excluded from CI).

## 5. Acceptance / regression gates

- Unit + integration suites pass with zero failures.
- `Rcpp::sourceCpp("src/sim.cpp")` compiles cleanly.
- Known regression: an all-age prison-only scenario must reduce total HCV
  more than a 40+-only scenario over the same horizon (guards the
  `tau_min_age` wiring that previously made all scenarios equivalent).
- The calibration acceptance criteria (see `requirements.md` section 3)
  are verified by `run_calibration.R` output; the Bayesian criteria
  (R-hat in [0.995, 1.005], pooled ESS > 1000) by `run_npe.R` output.

## 6. Out of CI scope

- Full multi-start calibration and full NPE training (hours of compute);
  these are run on demand via `scripts/run_calibration.R` and
  `scripts/run_npe.R` and reviewed from `output/calibration/`.
- Platform matrix testing (macOS/Windows) — currently ubuntu-only in CI.
