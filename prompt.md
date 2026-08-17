# prompt.md

You are the calibration agent for the Singapore HCV PWID model. Redo the
project from the clean baseline and follow `AGENTS.md`.

## Baseline

- Original progression rates: `0.027 / 0.039 / 0.014 / 0.014`.
- `m_min = m_max = 1`.
- 2015 Singapore mortality, scalar PWID SMR `omega = 14.68`.
- **No excess-mortality parameter**: `eta_s = 1`.
- 12 fitted parameters: 6 contact scaling factors, 6 beta scaling factors.
- **Simulation horizon**: at least 150 model years so that all states,
  including DC/HCC totals, converge before targets are extracted.
- **Beta scaling factors have no upper bound**; only a small positive lower
  bound is used.

## Contact-matrix preference (soft)

After scaling, the **diagonal** of `C_contact` should be:

- unimodal;
- peak as close as possible to 30-39 (or 20-29);
- 50+ and 60+ diagonal entries relatively small.

This is a soft plausibility penalty, not a hard constraint. A fit that is
statistically good and passes convergence is acceptable even if the diagonal
is not perfectly unimodal.

## Beta-scale preference (soft)

- Beta scales should be monotone non-decreasing with age.
- Values `>1` are preferred.
- Both are soft; never force them at the cost of a passing fit.

## Figure style

Every ggplot figure must:

- use a pure white background;
- have **no grid lines**;
- use 95% credible intervals in analysis figures;
- order analysis scenarios as in `analysis_report.md`;
- place the analysis legend on the right.

## Workflow

1. Delete previous results.
2. Run deterministic multi-start Nelder-Mead with the soft plausibility
   penalties.
3. Report equilibrium, acceptance criteria, Laplace intervals, fitted
   diagonal values, and beta scales.
4. Run NPE + strict MCMC with wide thinning and strict convergence.
5. After a passing MCMC result, run the 50-year sensitivity analysis in
   `prompt_analysis.md`.
6. Keep generated files minimal.

## Acceptance criteria

```text
prevalence RMSE                    <= 0.02
maximum absolute prevalence error <= 0.03
population MAPE                    <= 0.10
maximum population APE            <= 0.20
equilibrium                       pass
```
