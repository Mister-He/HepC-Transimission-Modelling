# Results summary (state of 2026-11-12)

Medians with 95% uncertainty intervals across 2,000 posterior draws (posterior parameter uncertainty only).
Sources: `output/scenarios_2026/`, `output/scenarios_2026_jcov/`, `output/scenarios_2026_linkage/`,
`output/decision_analysis/`, `output/tornado/`, `output/calibration/mcmc_reference/`.
Costs are the **programme budget only** (testing, consultations, DAA; 2026 SGD, discounted 3%); they exclude
care for advanced liver disease and savings from disease averted.

## Calibration and the 2015 population
- MCMC converged (largest R-hat 1.005; smallest bulk ESS 1,641); all 13 targets overlapped.
- 2015: 59,599 people who had ever injected (53,968–69,269): D 1,000, J 8,768, F 6,303, X 43,510.
  Community injectors (D + F) 7,332 (6,513–8,610) vs survey 11,000. First-arrest factor c 2.71 (1.54–5.23).
- 2015: 529 (454–628) new infections and 157 (130–191) liver deaths a year; flat without treatment.

## WHO 2030 targets (main analysis, treatment from 2026)
- **J only (status quo):** new infections −51.5% (40.6–57.3), liver deaths −22.1% (21.2–23.0). Fails both.
- **Finite screening:** 5 of 15 policies meet the incidence target (8,000/yr any F share; 5,000/yr with
  50% or 75% to F). **None meets the mortality target**; best 8,000/yr 25% to F −44.0% (41.4–45.9).
  No finite policy meets both in any draw.
- **Uncapped screening:** −99.8% infections, −71.4% (71.3–71.6) liver deaths; the only policy meeting both.
  Fails only if ψ_DC (liver death after cure in DC) is at the unfavourable end of its range (−59.4%).
- More tests to F → larger fall in infections, smaller fall in deaths. In 2026, 73% of chronic infections in
  F + X are in X, and F + X account for 90% of liver deaths.
- **Prison coverage:** every reduction rises with coverage; without prison treatment only 8,000/yr 75% to F
  meets incidence; uncapped screening just meets mortality (65.4%).
- **Prison treatment from 2019 (sensitivity):** J only −85.0% infections (meets), no capped policy meets
  mortality (best −56.1%).
- Genotype 3 share 10% / 50%: same pass sets. Observed inflow: see S1 Text Section 11.
- Budget 2026–2030: J only S$211.6M (180.3–252.4); cheapest incidence-passing policy (5,000/yr 75% F)
  S$244.3M, +S$32.5M.

## Long horizon, 2026–2050
- Median liver deaths reach 35% of 2015 in **2028** (uncapped), **2035** (8,000/yr 25% F), **2041**
  (8,000/yr 75% F), **2048** (J only). J-only infections reach 20% of 2015 in 2036.
- J only vs no treatment: 9,892 (8,279–11,843) infections and 1,651 (1,397–1,983) liver deaths averted;
  **93% of those liver deaths averted occur after 2030.**
- Extra liver deaths averted vs J only: uncapped 1,433 (1,176–1,788); 8,000/yr 25% F 947 (782–1,138);
  8,000/yr 75% F 546 (454–648).
- At every capacity, 25% to F averts more liver deaths; 75% to F averts more infections and needs fewer
  treatment courses.

## Cost per outcome averted vs J only, 2026–2050
| Policy | Extra budget, S$M | S$ per liver death averted | S$ per infection averted |
|---|---|---|---|
| 8,000/yr, 75% F | 43.7 (14.3–59.0) | ~117,000 (41,000–148,000) | ~17,500 |
| 8,000/yr, 50% F | 97.7 (67.8–120.4) | ~177,000 | ~44,100 |
| 8,000/yr, 25% F | 140.7 (111.9–170.8) | ~214,000 (186,000–225,000) | ~86,200 |
| Uncapped | 274.1 (224.9–338.8) | ~256,000 (232,000–269,000) | ~108,300 |
- 600 / 1,000 / 2,000 tests a year with 75% to F are **budget-neutral** vs J only (95% UI of the extra budget
  includes 0): the infections they prevent need no treatment later.

## Strategy selection (main decision horizon 2026–2050, 17 policies)
- **Frontier (liver deaths):** 1,000 → 2,000 → 5,000 → 8,000/yr (all 75% F) → 8,000/yr 50% F → uncapped.
  J only is dominated (1,000/yr 75% F costs no more and averts more). All 25%-to-F policies are dominated
  or extendedly dominated.
- **Rule (b), lowest cost per outcome averted vs no treatment:** 5,000/yr 75% F (both outcomes).
- **Rule (c), TOPSIS, equal weights:** 8 Pareto-efficient; first **8,000/yr 75% to F** (closeness 0.565;
  first in 99.2% of draws), then 5,000/yr 75% F (0.559), uncapped (0.551).
- **Over 2026–2030 only:** uncapped screening ranks first (only strategy meeting both targets).
- **Weights:** doubling the weight on liver deaths or on the WHO targets → uncapped first; doubling the
  budget or courses → 5,000/yr 75% F. Random weights: uncapped first in 48%, 5,000/yr 75% F 18%,
  8,000/yr 75% F 17%, 2,000/yr 75% F 15%.
- **Rule (a), budget caps (liver deaths):** S$300M → 5,000/yr 75% F; S$350M → 5,000/yr 50% F; S$400M →
  8,000/yr 50% F; S$450–500M → 8,000/yr 25% F; S$550M → uncapped.
- **68 combinations (prison coverage as a choice):** uncapped screening with full prison coverage first
  (100% of draws); top four are uncapped screening at each coverage; closeness rises with prison coverage.

## Linkage to care (sensitivity)
- 80% / 60% linkage: 5 / 3 finite policies meet incidence (main 5); none meets mortality.
- 8,000/yr 25% F: liver deaths −40.5% / −36.6% in 2030 (main −44.0%); 850 / 723 liver deaths averted vs J
  only to 2050 (main 947).
- Uncapped screening still meets both targets (−71.0% / −68.9%), because people missed are rescreened.

## Whom / where / how intensively
- **Whom:** released people — those still injecting (infections) and those who stopped (liver deaths);
  people never imprisoned (D) are not reached by any strategy.
- **Where:** prison treatment is the base; higher prison coverage is preferred in every comparison.
- **How intensively:** several thousand community tests a year; screening every released person each year
  is needed to meet the 2030 mortality target.
