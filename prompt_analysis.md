# prompt_analysis.md

Run this phase only after NM and MCMC pass all checks.

## Goal

Project the calibrated model 50 years from 2017 under a broad set of
treatment and screening strategies, propagate MCMC posterior uncertainty,
and recommend a practical strategy.

## Scenario configuration

Scenarios are stored in a single machine-readable table:

```text
src/calibration/scenarios.csv
```

Columns:

```text
scenario, description, tau_NC, tau_CC, tau_DC, tau_HCC, rho,
alpha_DC_mode, phi_mode, elig_D, elig_J, elig_F, elig_X, min_age_group
```

- `elig_D / elig_J / elig_F / elig_X`: 0/1 flags for the four model strata
  (D = never-arrested active, J = currently detained, F = ever-arrested
  active, X = former). Any combination is allowed, e.g. prison-only
  `(0,1,0,0)`, community `(1,0,1,0)`, post-release linkage `(0,0,1,0)`,
  prison plus post-release continuation `(0,1,1,0)`, all strata
  `(1,1,1,1)`.
- `min_age_group`: 1-based minimum eligible age group
  (e.g. `4` = treat only 40+; `5` = 50+; `6` = 60+).
- `alpha_DC_mode`: `pos` (with RBV) or `neg` (without RBV).
- `phi_mode`: `baseline` or `none` (disable post-SVR progression modifiers).
- `tau_*`: annual DAA treatment-initiation rate among eligible diagnosed
  chronic patients. Screening coverage is folded into `tau_*` (diagnosis is
  the gateway to treatment), so an opt-out screening programme is
  represented by a higher `tau_*` in the targeted strata/ages.

Adding/modifying a strategy only requires editing one CSV row; the analysis
script reads the table and applies it automatically. The scenario key
(`scenario` column) is the canonical lookup id used in all outputs, figures
and reports.

## Scenario inventory

The canonical inventory lives in `src/calibration/scenarios.csv`. It covers:

1. status quo no treatment;
2. community treatment (NC/CC only; broad all-stage);
3. community treatment age-targeted 40+;
4. prison-only treatment (all ages);
5. prison-only treatment 30+;
6. prison-only treatment 40+;
7. prison-only treatment 50+;
8. prison-only treatment 60+;
9. prison-only aggressive treatment;
10. prison universal opt-out screening + immediate DAA;
11. aggressive treatment in all strata;
12. post-release linkage (treat F only);
13. prison plus post-release continuation (J + F);
14. GT3 proportion low/high without treatment;
15. DC SVR without RBV;
16. post-SVR modifiers disabled;
17. combined broad community treatment + GT3 0.6;
18. prison-only 50+ with GT3 0.9.

Scenario keys:

```text
no_treatment, early_treatment_community, broad_treatment_community,
community_40plus, prison_only_treatment, prison_only_treatment_30plus,
prison_only_treatment_40plus, prison_only_treatment_50plus,
prison_only_treatment_60plus, prison_only_aggressive,
prison_universal_optout, aggressive_treatment_all,
post_release_linkage, prison_plus_postrelease,
rho_low_no_tx, rho_high_no_tx, dc_without_rbv,
no_postsvr_modifiers, broad_community_rho60, prison_only_50plus_rho90
```

## Screening interpretation

The model has no separate "screening" compartment: a screened-and-treated
patient moves from chronic (`c`) to treatment (`t`) at rate `tau` in the
eligible strata and age groups. Scenario tau values therefore represent the
annual probability that an eligible chronic patient is diagnosed **and**
initiates DAA. A screening-only programme without treatment has no effect in
this model and is reported as `no_treatment`; the practical analysis
discusses screening coverage, opt-out policy, and linkage as the drivers of
the effective `tau`.

## Required analyses

For each scenario report:

- total HCV (a+c+t), DC, and HCC counts per year, 2017-2067;
- median and 95% equal-tailed CrI at key years (2017, 2027, 2037, 2047,
  2057, 2067);
- comparison of age-targeted prison strategies to the all-age prison
  strategy in the medium term (2027-2037), where age effects are largest;
- scenario key-year table and CSV output;
- status quo always included as reference.

## Projection

- Run each posterior draw from `t_start = 0` to `t_end = 150` to full
  equilibrium.
- Use the equilibrated state as the 2017 initial condition.
- Project 2017-2067 under each scenario.

## Figure requirements

- Pure white background, no grid lines.
- 95% credible intervals.
- Status quo always included.
- Legend on the right, ordered as in `scenarios.csv`.

## Practical recommendation

After the projections, recommend which treatment/screening strategy to
implement, considering:

- prison population coverage and feasibility of prison-based DAA;
- age-targeted prioritisation (40+ / 50+ higher liver-disease burden);
- community vs prison-only reach and post-release linkage;
- cost, linkage-to-care, retention, and real-world DAA uptake;
- equity and public-health impact (HCV incidence, DC/HCC);
- order-of-magnitude robustness: strategies that still eliminate by 2067
  under pessimistic GT3 / SVR assumptions.

Cite literature/websites for each practical claim (e.g. WHO HCV guidelines,
prison-based HCV treatment reviews, PWID DAA programmes).

## Outputs

```text
output/analysis/scenario_summary.csv   (per-scenario yearly median/CrIs)
output/analysis/scenario_key_years.csv (2017-2067 key-year table)
output/analysis/fig_hcv_trajectories.png
output/analysis/fig_dc_hcc_trajectories.png
figures/fig_hcv_trajectories.png       (stable copies)
figures/fig_dc_hcc_trajectories.png
docs/calibration/analysis_report.md    (method + results + recommendation)
```

Keep the scenario inventory, results, and recommendation consistent across
these files.
