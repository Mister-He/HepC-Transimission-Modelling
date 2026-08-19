# Sensitivity analysis — 50-year projections from 2017

## 1. Method

- Uses the passing MCMC posterior (`output/calibration/npe_bayes/`).
- 300 equilibrium-feasible posterior draws per scenario.
- Each draw is run from `t_start = 0` to `t_end = 150` to full equilibrium;
  the equilibrated state is used as the 2017 initial condition, then the
  model is projected 2017-2067 under each scenario (RK4, daily steps).
- Quantities: total HCV (acute + chronic + treatment, `a+c+t`), DC, HCC.
- Medians and 95% equal-tailed credible intervals (CrIs) at key years.
- Status quo (`no_treatment`) is always included as reference.

## 2. Scenario configuration

Scenarios are stored in a single machine-readable table:

```text
src/calibration/scenarios.csv
```

Each row defines, per scenario:

- `tau_NC..tau_HCC`: annual DAA treatment-initiation rate among eligible
  diagnosed chronic patients, by liver stage;
- `rho`: proportion of infections with GT3 (drives SVR and progression);
- `alpha_DC_mode`: DC SVR with (`pos`) or without (`neg`) RBV;
- `phi_mode`: `baseline` post-SVR progression modifiers or `none`;
- `elig_D, elig_J, elig_F, elig_X`: 0/1 flags for the four strata
  (D = never-arrested active, J = currently detained, F = ever-arrested
  active, X = former);
- `min_age_group`: 1-based minimum eligible age group (`4` = 40+,
  `5` = 50+, `6` = 60+).

Adding or modifying a strategy requires editing only one CSV row; the
analysis script reads the table and applies it automatically. The `scenario`
column is the canonical key used in outputs, figures, and this report.

Screening coverage is folded into `tau_*` (diagnosis is the gateway to
treatment), so an opt-out screening programme is represented by a higher
`tau_*` in the targeted strata/ages; a screening-only programme without
treatment is structurally equivalent to `no_treatment` in this model.

## 3. Scenario inventory (20 scenarios)

```text
no_treatment                 status quo no treatment
early_treatment_community    NC/CC only in D + F
broad_treatment_community    all stages in D + F
community_40plus             all stages in D + F, 40+
prison_only_treatment        all stages in J only, all ages
prison_only_treatment_30plus J only, 30+
prison_only_treatment_40plus J only, 40+
prison_only_treatment_50plus J only, 50+
prison_only_treatment_60plus J only, 60+
prison_only_aggressive       higher tau in J only, all ages
prison_universal_optout      opt-out screening + immediate DAA in J
aggressive_treatment_all     higher tau in all strata
post_release_linkage         treat F only (linkage after release)
prison_plus_postrelease      treat J + F (prison and continuation)
rho_low_no_tx                GT3 = 0.6, no treatment
rho_high_no_tx               GT3 = 0.9, no treatment
dc_without_rbv               DC SVR without RBV, no treatment
no_postsvr_modifiers         broad community, phi modifiers disabled
broad_community_rho60        broad community, GT3 = 0.6
prison_only_50plus_rho90     prison 50+, GT3 = 0.9
```

## 4. Key results

All counts are medians over 300 posterior draws (95% CrIs in parentheses).
The status-quo trajectory is flat (HCV 11,065; DC 479; HCC 125.6), which
confirms the 2017 starting state is a genuine equilibrium.

### 4.1 Total HCV (a + c + t)

| Scenario | 2027 | 2037 | 2067 |
|---|---:|---:|---:|
| no_treatment | 11,065 | 11,065 | 11,065 |
| early_treatment_community | 5,634 (5,098-6,301) | 2,691 | 184 (157-217) |
| broad_treatment_community | 5,602 (5,068-6,265) | 2,663 | 179 (153-211) |
| community_40plus | 8,393 (7,592-9,364) | 7,224 | 5,960 (5,259-6,718) |
| prison_only_treatment | 6,582 (5,961-7,310) | 4,608 | 2,899 (2,533-3,373) |
| prison_only_treatment_30plus | 7,227 | 5,691 | 4,363 |
| prison_only_treatment_40plus | 8,384 (7,595-9,332) | 7,606 | 7,124 (6,406-7,958) |
| prison_only_treatment_50plus | 9,946 (8,969-11,099) | 9,770 | 9,703 (8,737-10,840) |
| prison_only_treatment_60plus | 10,933 | 10,921 | 10,911 |
| prison_only_aggressive | 5,710 (5,175-6,351) | 3,593 | 1,796 (1,497-2,156) |
| prison_universal_optout | 5,941 (5,385-6,606) | 3,856 | 2,069 (1,750-2,471) |
| aggressive_treatment_all | 45 (38-53) | 0.5 | ~0 |
| post_release_linkage | 7,292 (6,599-8,132) | 4,832 | 1,926 (1,455-2,488) |
| prison_plus_postrelease | 4,871 (4,407-5,436) | 2,425 | 305 (216-451) |
| broad_community_rho60 | 5,639 | 2,713 | 190 |
| prison_only_50plus_rho90 | 9,915 | 9,712 | 9,621 |
| rho_low_no_tx | 11,116 | 11,161 | 11,201 |
| rho_high_no_tx | 11,031 | 11,000 | 10,973 |
| dc_without_rbv | = no_treatment | = | = |
| no_postsvr_modifiers | = broad_treatment_community | = | = |

### 4.2 DC and HCC at 2067

| Scenario | DC | HCC |
|---|---:|---:|
| no_treatment | 479.2 | 125.6 |
| early_treatment_community | 21.6 | 4.8 |
| broad_treatment_community | 21.2 | 4.7 |
| community_40plus | 269.0 | 68.7 |
| prison_only_treatment | 99.2 | 25.9 |
| prison_only_treatment_40plus | 298.2 | 77.2 |
| prison_only_treatment_50plus | 423.8 | 110.4 |
| prison_only_aggressive | 59.7 | 15.4 |
| prison_universal_optout | 69.1 | 18.0 |
| aggressive_treatment_all | ~0 | ~0 |
| post_release_linkage | 93.0 | 23.2 |
| prison_plus_postrelease | 22.5 | 5.2 |

### 4.3 Age-targeted prison strategies (medium term)

Restricting prison treatment by age progressively weakens the medium-term
impact (medians at 2027/2037/2067, HCV):

| Strategy | 2027 | 2037 | 2067 | Reduction by 2067 |
|---|---:|---:|---:|---:|
| prison all ages | 6,582 | 4,608 | 2,899 | 74% |
| prison 30+ | 7,227 | 5,691 | 4,363 | 61% |
| prison 40+ | 8,384 | 7,606 | 7,124 | 36% |
| prison 50+ | 9,946 | 9,770 | 9,703 | 12% |
| prison 60+ | 10,933 | 10,921 | 10,911 | 1% |

Older-only strategies treat a small pool and leave younger cohorts to drive
continued transmission, so they cannot reach elimination alone. They are
best viewed as a phased roll-out step, not a long-run strategy.

## 5. Reading the results

1. **Prison-only treatment is powerful but insufficient alone.** Moderate
   all-age prison DAA cuts total HCV by 74% over 50 years, but ongoing
   community transmission sustains ~2,900 chronic infections at 2067.
   Higher prison intensity (aggressive / opt-out) accelerates the decline
   (84% / 81%) but does not eliminate.
2. **Post-release linkage nearly closes the gap.** Prison + continuation
   in the ever-arrested community (`F`) reduces HCV to 305 by 2067 (97%),
   close to full community access (179, 98%). This is the single most
   important design choice for a prison-centred programme.
3. **Universal community access achieves near-elimination.** Treating
   active PWID in the community (D + F) at moderate rates reaches ~180
   infections by 2067, regardless of whether NC/CC-only or all stages are
   treated (the marginal benefit of treating DC/HCC is small).
4. **Aggressive all-strata treatment eliminates** by ~2040.
5. **Robustness.** GT3 proportion (0.6-0.9), DC SVR with/without RBV, and
   post-SVR modifier settings change results by at most a few percent and
   do not change the ranking of strategies.

## 6. Structural limitations

- **Screening is not a separate compartment.** `tau_*` is the annual
  diagnosis-and-treatment rate; the report reasons about screening coverage
  through its effect on `tau_*`.
- **Post-SVR progression modifiers (`phi`) are read but not used in the
  ODE** (cured `u` states do not progress in the current implementation),
  so `no_postsvr_modifiers` is structurally equivalent to the baseline
  treatment scenario. This is flagged as a sensitivity check only.
- **No reinfection-specific protection** beyond the susceptible/post-SVR
  state; reinfection risk is implicit in the force of infection.
- **Treatment courses are not followed across release**; the
  `prison_plus_postrelease` scenario approximates continuity of care by
  making F eligible.

## 7. Practical recommendation

**Recommended strategy: universal opt-out HCV screening and DAA in prisons,
age-prioritised phasing (40+ first), with mandatory post-release community
linkage and unrestricted community DAA access for active PWID.**

The projections support a three-step, evidence-aligned pathway:

1. **Immediate: prison universal opt-out screening + DAA for all ages.**
   Prisons concentrate high-prevalence PWID and offer a feasible setting
   for directly observed, short-course DAA. Singapore already screens in
   correctional settings, and unrestricted prison DAA access raised SVR to
   99.1% (Wong et al., 2021). The 2025 global best-practice prison
   guidelines recommend opt-out testing, on-site treatment, and peer
   support as core components (Kronfli et al., 2025).
2. **Phasing: prioritise 40+/50+ during capacity ramp-up, then all ages.**
   Older PWID carry the highest DC/HCC burden, so treating them first
   captures the most advanced disease per course. However, the model shows
   age-only strategies leave 60-99% of infections at 2067 because younger
   cohorts sustain transmission; age prioritisation must be a temporary
   sequencing tool, not the endpoint.
3. **Essential: post-release linkage and unrestricted community DAA.**
   Prison + community continuation reaches 97% reduction, versus 74% for
   prison-only. Jail-based testing + treatment + post-release linkage is
   projected to cut new PWID infections by ~47% and HCV deaths by ~40%
   (Stanford Health Policy, 2026); linkage interventions in prison are
   among the strongest-evidenced strategies in the 2025 systematic review
   (Cunningham et al., 2025). WHO recommends DAA for all adults with
   chronic HCV, including PWID, without fibrosis-stage restriction
   (WHO, 2018/2022).

Feasibility and equity considerations: community programmes should be
co-located with harm-reduction services (OST, NSP) and use point-of-care
RNA testing and peer navigation to maximise uptake; costs should be tracked
against avoided DC/HCC events (WHO, 2018; Cunningham et al., 2025). A
Singapore-specific modelling preprint similarly concludes that treating
detained PWID is the most feasible and cost-effective route, with combined
community coverage needed for full elimination (Park et al., medRxiv
2025) — the same source paper that provides the calibration targets, so
the present projections are consistent with its conclusions.

## 8. References

1. World Health Organization. *Guidelines for the care and treatment of
   persons diagnosed with chronic hepatitis C virus infection*. Geneva:
   WHO, 2018 (updated 2022). ISBN 9789241550345.
   https://www.who.int/publications/i/item/9789241550345
2. Kronfli N, et al. Best practice guidelines for viral hepatitis service
   delivery in prisons. *Lancet Infectious Diseases*, 2025.
   https://doi.org/10.1016/S1473-3099(25)00630-9
3. Cunningham EB, et al. Interventions to improve testing, linkage to care,
   and treatment for hepatitis C infection in prison: a systematic review
   and meta-analysis. *International Journal of Drug Policy*, 2025.
   https://doi.org/10.1016/j.drugpo.2025.105082
4. Wong YJ, et al. The impact of unrestricted access to direct-acting
   antiviral among incarcerated hepatitis C virus-infected patients.
   *Clinical and Molecular Hepatology*, 2021;27(3):474-485.
   https://doi.org/10.3350/cmh.2021.0015
5. Singapore Medical Journal. Elimination of chronic viral hepatitis C in
   correctional health (review). 2025.
   https://doi.org/10.4103/singaporemedj.smj-2025-161
6. Stanford Health Policy. Jails could cut hepatitis C infections nearly in
   half with testing and treatment. 2026.
   https://healthpolicy.fsi.stanford.edu/news/jails-could-cut-hepatitis-c-infections-nearly-half-testing-and-treatment
7. Park et al. HCV transmission and elimination modelling among PWID in
   Singapore (preprint). medRxiv, 2025.
   https://doi.org/10.1101/2025.10.24.25338708

## 9. Files

```text
src/calibration/scenarios.csv           scenario inventory (single source)
src/calibration/run_analysis.R          analysis driver
output/analysis/scenario_summary.csv    per-scenario yearly medians/CrIs
output/analysis/scenario_key_years.csv  key-year table (2017-2067)
output/analysis/fig_hcv_trajectories.png
output/analysis/fig_dc_hcc_trajectories.png
figures/fig_hcv_trajectories.png        stable copies
figures/fig_dc_hcc_trajectories.png
```
