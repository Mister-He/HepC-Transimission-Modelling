# HCV transmission model among PWID in Singapore — NumPyro/JAX edition

An age-structured compartmental model of hepatitis C (HCV) transmission among
people who inject drugs (PWID) in Singapore, migrated from the legacy R/C++
implementation to **NumPyro / JAX**, with Bayesian inference via **NUTS**.

Chinese version: [README.zh-CN.md](README.zh-CN.md). Full experiment record:
[memo.md](memo.md). Design philosophy (Chinese): [docs/design.md](docs/design.md).

## Model

- 4 strata (D/J/F/X) × 4 liver stages (NC/CC/DC/HCC) × 4 HCV states (u/a/c/t)
  × 6 age groups = **384 compartments**;
- continuous-time ODE, fixed-step RK4 (`dt = 1/365`), 10-year age-band
  progression, simulated to t = 150 (targets at t = 45 ≈ 2015);
- frequency-dependent force of infection among active injectors (D, F);
- **12 fitted log-parameters**: 6 contact-matrix row scales + 6 beta inflow
  scales, calibrated to prison (J-stratum) HCV seroprevalence and population.

## Layout

```text
config/    fixed parameters as JSON (model, simulation, targets, priors)
src/       JAX simulator (sim.cpp port) + NumPyro model + inference tools
scripts/   run_nuts.py — NUTS inference pipeline
tests/     unit (incl. R-reference equivalence) + integration (smoke)
legacy/    read-only archive of the R/C++ results and source
docs/      Chinese documentation (requirements/architecture/test-plan/design)
results/   inference outputs (posterior, diagnostics, PPC, figures)
```

## Quickstart

```bash
python3.13 -m venv .venv && source .venv/bin/activate
pip install -r requirements-dev.txt

pytest tests/unit -v          # incl. numerical equivalence with R reference
pytest tests/integration -v   # short-horizon NUTS/PPC smoke

python scripts/run_nuts.py --num-warmup 500 --num-samples 1000 \
  --num-chains 4 --dt-days 7 --out-dir results/numpyro_nuts
```

For parallel chains set `XLA_FLAGS="--xla_force_host_platform_device_count=4"`.

## Key results (see memo.md and results/numpyro_nuts/)

- JAX port matches the R/C++ reference to ~1e-12 at the legacy NM optimum;
- NUTS (whitened parameterization, NM point-estimate init) produces a posterior
  with R-hat/ESS meeting the same strict criteria as the legacy run;
- posterior predictive checks cover the observed targets; per-parameter
  comparison with the legacy posterior (medians, 95% intervals, KS tests) is
  reported in `results/numpyro_nuts/compare_legacy.csv`.

## License / provenance

Derived from the legacy branch `dev_10yrs_age_interval`; all legacy results and
source are archived under `legacy/`.
