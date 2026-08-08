# Model audit (revision v1)

Date: 2026-08-08. Read-only audit of the current `src/sim.cpp` /
`src/setup.R` implementation.

## Compartment structure

```text
3 strata   D (community PWID, inject), J (detained), X (former PWID, non-injecting)
4 liver stages  NC, CC, DC, HCC
5 HCV states    u (susceptible), a (acute), c (chronic), t (treatment), s (seropositive cleared/post-SVR)
6 age groups    <20, 20-29, 30-39, 40-49, 50-59, 60+
360 compartments
```

Ageing: `y_change = y[base + i] / 10.0 * dt;` — the corrected 10-year-band
implementation (verified present; 1/10 per year flows to the next band,
with the 60+ band open-ended).

## Transmission

- Force of infection exists only in D (detained PWID do not inject; X are
  ex-PWID): `gamma_i(t) = q * sum_j C_contact(i,j) * infectious_D_j / active_D_j`
  with infectious = acute + chronic. Proportionate (frequency-dependent)
  mixing within D.
- This revision: constant transmission (`m(t) = 1`); the historical level
  is absorbed by fitted contact row scales.

## Flows

- Inflow `beta_i` into `D_{u,NC,i}` (susceptible entrants).
- Arrest D -> J at `lambda1[i]`; release J -> D with probability
  `pi_recid` or J -> X with `1 - pi_recid` at `lambda2[i]`.
- HCV: acute -> chronic/clearance; chronic -> CC/DC/HCC progression with
  genotype-weighted rates; treatment `tau = 0` (no treatment scenario);
  post-SVR progression modifiers `phi_*`.
- Mortality: `mu[i] * omega` background (omega = 14.68); DC adds
  `mu_DC = 0.130`, HCC adds `mu_HCC = 0.430`; the cleared (s) state uses
  `mu[i] * omega * eta_s[i]` with `eta_s = 1` unless the excess-mortality
  contingency is activated.

## Targets

The Changi Prison universal screening (Dec 2014 - Feb 2016) is anti-HCV
serology (non-/borderline-/reactive), so the target is seroprevalence
`(a+c+t+s)/N`; the `s` state makes seroprevalence exactly representable
(`u` = never infected only). Viremic prevalence `(a+c+t)/N` is reported as
a sensitivity. Integer positives reconstructed as
`round(target_prev * prison_total)` (recorded in `targets.csv`).

## Equilibrium

The system is simulated over t = -10..55 (calendar 1960-2025, t = 0 <->
1970), targets at t = 45 (2015), and the equilibrium gate is evaluated at
t = 55 vs t = 50 (2025 vs 2020). With constant transmission the age
structure is a genuine quasi-steady state; the gate must still pass
(criteria 0.01 log-pop, 0.005 prevalence).

## Diagnosis for the 12-parameter fit (a priori)

1. Age-pattern: with constant transmission, prevalence rises with age
   (cumulative exposure): <20 low (recent entrants), 20-49 rising,
   50-59 plateau ~0.48. This matches the observed pattern up to 50-59.
2. The 60+ decline (0.48 -> 0.355) must come from either (a) uninfected
   elderly inflow (beta6) diluting J 60+, or (b) excess mortality of
   seropositives. With `eta_s = 1`, only (a) is available; the
   contact-6/beta6 trade-off determines whether (a) alone can satisfy both
   N6 = 409 and p6 = 0.3545. If not, the excess-mortality contingency
   (dedicated rationale file) is invoked.
3. <20 saturation: contact row 1 (young PWID mixing) plus inflow beta1
   control the young plateau (p1 ~ 0.111, N1 = 99). The steady-state
   balance implies a modest beta1 inflow; whether beta_scale1 > 1 is
   feasible is tested explicitly (soft rule).
4. Fast liver progression (adopted rates) removes chronics in older
   groups; this slightly lowers middle/older prevalence and is expected to
   be largely absorbed by re-optimisation of contact scales.
