# Model card

## Purpose

This deterministic compartmental model represents HCV transmission and liver
disease progression among Singaporean PWID moving through four strata: never
incarcerated active PWID (`D`), incarcerated (`J`), formerly incarcerated
active PWID (`F`), and formerly incarcerated ex-PWID (`X`). It has six age
groups, four liver stages, and five HCV states (`u`, `a`, `c`, `t`, `r`).

Recovered people remain antibody-reactive and can be reinfected. Calibration
seroprevalence is `(r+a+c+t)/(u+r+a+c+t)`.

## Intended use

- estimate the 2015 PWID/prison equilibrium from the supplied age targets;
- project age-specific drug-abuser inflow through 2030; and
- compare prison, community, mixed, age-targeted, and stage-targeted treatment
  strategies beginning in 2026.
- evaluate modelled-PWID incidence and HCV-related mortality reductions against
  the WHO 2030 relative targets.

## Limitations

This is not a national HCV model. It does not represent the general population
or every non-PWID transmission route, and SingStat new-drug-abuser counts are
an enforcement-based inflow proxy. The WHO incidence and mortality indicators
are model-specific PWID-network proxies, not national estimates: incidence is
new infection flow into acute infection, while HCV-related mortality is the
DC/HCC disease-specific mortality component. Results are policy comparisons
within the modelled PWID network.

Calibration limits:
- **Targets:** only the 2015 age-specific prison seroprevalence is fitted (Beta likelihood, φ =
  100); all six are within 0.022. Population size is not fitted: the modelled 2015 prison
  population (8,451) is 2.5 times the observed count of screened PWID (3,367).
- **Rate fits:** release and re-arrest rates come from SPS survival models. First arrest is
  `c` × re-arrest with `c` ≥ 0.5; seroprevalence barely informs `c`, whose posterior (1.24,
  0.57–3.13) follows its prior.
- **Entry:** the CNB new-drug-abuser counts enter D unscaled.
- **Mortality:** non-liver deaths use all-cause PWID SMRs against the general population from an
  Amsterdam cohort grouped by antibody status (9.52 HCV-positive, 4.73 negative; van Santen 2014). HCV deaths include the viraemic excess,
  which cure removes; any confounding in that excess overstates the benefit of treatment.
- **Prior-dominated parameters:** the 50–59 and 60+ contact scales and `c`. The
production posterior is the reference MCMC; single-round NPE is used only as an
independent cross-check.
