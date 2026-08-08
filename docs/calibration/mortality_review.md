# Mortality evidence review — 2015 Singapore baseline

Date: 2026-08-08. Status: primary specification.

## Singapore 2015 age-specific death rates (SingStat)

Source: Singapore Department of Statistics, "Age-Specific Death Rates,
Annual" (data.gov.sg), reference year 2015, total residents. Mapped to the
six model age groups:

| Model age group | SingStat 2015 bands used | mu (per year) |
|---|---|---:|
| <20 | 15-19 | 0.0002 |
| 20-29 | mean(20-24, 25-29) | 0.00025 |
| 30-39 | mean(30-34, 35-39) | 0.00045 |
| 40-49 | mean(40-44, 45-49) | 0.0012 |
| 50-59 | mean(50-54, 55-59) | 0.00335 |
| 60+ | population-weighted mean over 60-64 .. 90+ (deaths 60+/population 60+ = 15,322/700,208) | 0.0218818 |

The 60+ open-ended rate is population-weighted (0.02188/yr); using only
60-69 would understate mortality of the open-ended group (older
alternative 0.00845/yr) and was rejected (previous-round finding, retained
in this review).

## PWID excess mortality (SMR)

Source: Mathers et al. 2013, "Mortality among people who inject drugs: a
systematic review and meta-analysis", Bull World Health Organ 91:102-123,
PMID 23554523. Pooled all-cause SMR for PWID = **14.68** (95% CI
13.3-16.1). Applied as scalar `omega = 14.68` in
`effective_mu[i] = singapore_mu[i] * omega` (model code: `mu[i] * omega`).

The SMR is a pooled (not age-stratified) estimate; age-stratified SMR
scenarios remain a documented sensitivity for future work (previous-round
M2/Glasgow scenario showed an age-declining SMR cannot reproduce the 60+
prison share under this model).

## Disease-specific mortality (retained)

| Parameter | Value | Source |
|---|---:|---|
| `mu_DC` | 0.130 /yr | Lim et al. 2018 (additional mortality in decompensated cirrhosis) |
| `mu_HCC` | 0.430 /yr | Lim et al. 2018 (additional mortality in HCC) |
| `psi_DC` / `psi_HCC` | 0.45 / 0.37 | Premkumar 2024; Dang et al. 2020 (post-SVR mortality modifiers) |

These are disease-specific and unchanged in this revision.
