# Liver-disease progression rates — evidence review (revision v1)

Date: 2026-08-08. Status: adopted into `src/setup.R` (highest defensible
literature values, per user directive). See DECISIONS.md
"Progression rates: highest defensible literature values".

## Adopted rates

| Parameter | Value (/yr) | Mean sojourn | Source |
|---|---:|---:|---|
| `p_NC_CC` | 0.027 | 37 yr | Thein et al. 2008 (retained) |
| `p_CC_DC` | 0.0788 | 12.7 yr | Alazawi et al. 2010, Aliment Pharmacol Ther 32(3):344-55, PMID 20497143 — untreated-only decompensation in compensated HCV cirrhosis, 7.88%/yr (pooled 6.37%/yr; treated-included 5.34%/yr) |
| `p_CC_HCC` | 0.0479 | 20.9 yr | Same review — untreated-only HCC, 4.79%/yr (pooled 3.36%/yr; treated-included 2.52%/yr) |
| `p_DC_HCC` | 0.0464 | 21.6 yr | Rivera-Irigoin et al. 2006, AIDS Res Hum Retroviruses 22(12):1236-41, PMID 17209765 — HCC in HCV-monoinfected decompensated cirrhosis 3.31/100 py, 95% CI 2.70-4.64 (upper bound adopted) |

Why the untreated-only values: the model has no treatment (`tau = 0`), so
the comparator arm that excludes treated patients is the relevant one. The
Alazawi review explicitly separates studies including treated patients
(lower rates) from untreated-only studies (higher rates).

## Genotype-3 weighting

`src/sim.cpp` applies genotype-weighted effective rates:

```cpp
ptc_CC_DC  = (rho*r3_CC_DC  + (1-rho)) * p_CC_DC;   // rho = 0.78, r3 = 1.36
ptc_CC_HCC = (rho*r3_CC_HCC + (1-rho)) * p_CC_HCC;  // r3 = 1.93
ptc_DC_HCC = (rho*r3_DC_HCC + (1-rho)) * p_DC_HCC;  // r3 = 1.93
```

Effective annual rates with the adopted baselines:

| Transition | Baseline | Effective (GT-weighted) | Literature comparison |
|---|---:|---:|---|
| CC -> DC | 0.0788 | 0.101 | above pooled 6.37%/yr and untreated 7.88%/yr (GT3 RR 1.36 applied on top of a mixed-genotype baseline) |
| CC -> HCC | 0.0479 | 0.083 | within Asian cirrhosis-cohort ranges (4-10%/yr) but above Western pooled 3.36-4.79%/yr |
| DC -> HCC | 0.0464 | 0.080 | above the 3.31/100 py (CI 2.70-4.64) mixed-genotype estimate |

Interpretation: the adopted baselines are the highest directly sourced
values, and after GT3 weighting they sit at or above the upper end of the
international literature. This is intentional per the user directive
("as high as possible"); the consequence is that the model errs toward fast
liver-disease progression, which mainly affects the older chronic
population.

## Supporting context

- Fattovich et al. 1997, Gastroenterology 112(2):463-72, PMID 9024300
  (384 compensated HCV cirrhosis patients): 18% decompensation and 7% HCC
  at 5 years (~3.6-3.9%/yr and ~1.4%/yr) — older, interferon-era, and
  lower than the adopted values.
- Post-DAA meta-analyses report HCC ~1.5-3.0 per 100 person-years after SVR
  in advanced fibrosis/cirrhosis — lower than untreated, consistent with
  the adopted untreated rates being an upper bound.
- Japanese decompensated-HCV-cirrhosis cohorts report even higher
  post-decompensation HCC (cumulative ~14% at 1 yr in one cohort of 1412
  patients), but these are older, GT1b-dominated populations and are used
  only as context, not adopted.

## Caveats

1. Rates vs probabilities: the model applies these as annual hazard-style
   rates; the one-year transition probability is `1 - exp(-rate)` (~7.6%,
   4.7%, 4.5% respectively), i.e. close to the stated annual percentages.
2. Age dependence: cirrhosis complications are age-dependent; constant
   rates applied to all ages slightly overstate progression in younger
   PWID. A fully age-stratified specification was not adopted (transparent
   constant rates preferred; flagged as a limitation).
3. Double-counting risk with GT3 RRs: the literature rates are mixed
   genotype; multiplying by GT3 RRs on top inflates effective rates. The
   chosen values are therefore deliberately conservative-upper.
4. These are progression (natural-history) parameters, not transmission
   probabilities (transmission is governed by `q`, `C_contact`, `m(t)=1`).
