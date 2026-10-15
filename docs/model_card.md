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
- **2015 targets:** 60+ seroprevalence is 0.455 against 0.429 (0.65 SE); the other targets are
  within 0.01. The CNB counts anchor the 2015 entry into D (log-normal, SD 1.0).
- **Rate fits:** the first-arrest flow and the release and re-arrest rates come from SPS
  records, under a steady-state prison balance. The modelled PWID population is the people with
  a known screening result.
- **Excess mortality:** non-liver deaths of chronic viraemic people use an HR of 1.47 from a
  Taiwanese general-population cohort (Lee et al. 2012), not from PWID.
- **Prior-dominated parameters:** the 60+ contact scale and the older entry scales. The
production posterior is the reference MCMC; single-round NPE is used only as an
independent cross-check.
